rm(list=ls())
library(dplyr)
library(ggplot2)
library(MatchIt) 
library(multiwayvcov)
library(plm)
library(lmtest)
library(clubSandwich)
library(moments)
library(doParallel)
set.seed(123456789) #not needed for final version?

### this is executed in the /report subdirectory, need to ..
path <- strsplit(getwd(), "/report")[[1]]

### set this switch to TRUE if you want to produce a final report - this will save results matrices in a static directory
final_verion_swith <- TRUE
## heterogeneity analysis:
# 0 no
# 1 allow for enough time - sc level 
hetero <- 0
RI_conf_switch <- TRUE
glob_repli <- 2500
glob_sig <- c(.025,.975) ### 5 percent conf intervals


# takes raw data (baseline and endline), makes it anonymous and puts in into the data/public folder, ready to be analysed by the code chucks below
#source("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/cleaning.R")
#source("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/anonyize.R")
endline <- read.csv(paste(path,"data/public/endline.csv", sep="/"))
endline$a21 <- as.character(endline$region)
endline$region <- NULL
### EDITS SHOULD BE MADE HERE FOR FINAL VERSION###########################################################################################
### to do in final version: #
### - check skewness before taking IHS.
### - variables in endline that were transformed should also be transformed in baseline
### - number of days sick can only be determined if all data is in, as we need to establish the max number of sick household members
### there should be no duplicates in this dataset
endline <- endline[!duplicated(endline$hhid),]

endline$a21 <- as.factor(endline$a21)

########################################################### functions declarations #####################################################

trim <- function(var, dataset, trim_perc=.05) {
### function for triming a variable in a dataset - replaces with NA
 dataset[var][dataset[var] < quantile(dataset[var],c(trim_perc/2,1-(trim_perc/2)), na.rm=T)[1] | dataset[var] > quantile(dataset[var], c(trim_perc/2,1-(trim_perc/2)),na.rm=T)[2] ] <- NA
return(dataset)
}

FW_index <- function(indexer,revcols = NULL,data_orig) {
### function to make family wise index using covariance as weights (following http://cyrussamii.com/?p=2656)
### FW_index(c("baraza.B2","baraza.B3","baraza.B4.1","inputs","baraza.B5.2","baraza.B5.3"),data=endline)
##indexer <- c("baraza.B2","baraza.B3","baraza.B4.1","inputs","baraza.B5.2","baraza.B5.3")
##data_orig <- endline

data <- data_orig[complete.cases(data_orig[indexer]),]
x <- data[indexer]
  					if(length(revcols)>0){
						x[,revcols] <-  -1*x[,revcols]
					}

				for(j in 1:ncol(x)){
					x[,j] <- (x[,j] - mean(x[,j]))/sd(x[,j])
				}

					i.vec <- as.matrix(rep(1,ncol(x)))
					Sx <- cov(x)

					data$index <- t(solve(t(i.vec)%*%solve(Sx)%*%i.vec)%*%t(i.vec)%*%solve(Sx)%*%t(x))
data <- merge(data_orig,data[c("hhid","index")], by="hhid",all.x=T)

return( data )
}

### wrapper function to make summary forest plots
credplot.gg <- function(d,units, hypo, axlabs, lim){
 # d is a data frame with 4 columns
 # d$x gives variable names
 # d$y gives center point
 # d$ylo gives lower limits
 # d$yhi gives upper limits
 require(ggplot2)
 p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi, colour=as.factor(grp)))+
 geom_pointrange(position=position_dodge(-.4), size=.5)+
 geom_hline(yintercept = 0, linetype=2)+
 coord_flip(ylim = c(-lim,lim))+
 xlab('') + ylab(units)+ labs(title=hypo)  + theme_minimal()+ theme(axis.text=element_text(size=18),
        axis.title=element_text(size=14,face="bold"),legend.text=element_text(size=18), plot.title = element_text(size=22,hjust = 0.5), legend.title=element_blank())+
    geom_errorbar(aes(ymin=ylo, ymax=yhi),position=position_dodge(-.4),width=0,cex=1.5) + scale_colour_manual(values = c("#CCCCCC", "#6E8DAB", "#104E8B", "#000000")) + scale_x_discrete(labels=axlabs)
 return(p)
}
# function definitions 
RI_conf_sc <- function(i,outcomes, baseline_outcomes, dta_sim , ctrls = NULL, nr_repl = 1000, sig = c(.025,.975)) {
### a function to esimate confidence intervals using randomization inference following Gerber and Green pages 66-71 and 83.
#RI_conf_dist(2,outcomes, baseline_outcomes, subset(dta, ((information == 1 & deliberation==1) | district_baraza == 1)) , ctrls = "a21", nr_repl = 1000, sig = c(.025,.975))
#dta_sim <- dta[dta$district_baraza == 0 ,]
#ctrls <- "a21"
#nr_repl <- 1000
#sig <-  c(.025,.975)

	if (is.null(baseline_outcomes)) {
		formula1 <- as.formula(paste(outcomes[i],paste("information:deliberation",ctrls,sep="+"),sep="~"))
		formula2 <- as.formula(paste(outcomes[i],paste("information*deliberation",ctrls,sep="+"),sep="~"))
	} else {
		formula1 <- as.formula(paste(paste(outcomes[i],paste("information:deliberation",ctrls,sep="+"),sep="~"),baseline_outcomes[i],sep="+"))
		formula2 <- as.formula(paste(paste(outcomes[i],paste("information*deliberation",ctrls,sep="+"),sep="~"),baseline_outcomes[i],sep="+"))
	}


	dta_sim <- dta_sim %>%  mutate(clusterID = group_indices(., district, subcounty))
	### get ATEs for two different models
	ols_1 <- lm(formula1, data=dta_sim[dta_sim$deliberation == dta_sim$information,]) 
	ols_2 <- lm(formula2, data=dta_sim) 
	dta_sim$dep <- as.numeric(unlist(dta_sim[as.character(formula1[[2]])]))

	treat_nrs <- table(data.frame(aggregate(dta_sim[c("information","deliberation")], list(dta_sim$clusterID),mean))[,2:3])

	### calculate potential outcomes
	### for model 1 (sc effect)
	dta_sim$pot_out_0_1 <- NA
	dta_sim$pot_out_0_1[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0]
	dta_sim$pot_out_0_1[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] - coef(ols_1)["information:deliberation"]

	dta_sim$pot_out_1_1 <- NA
	dta_sim$pot_out_1_1[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1]
	dta_sim$pot_out_1_1[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] + coef(ols_1)["information:deliberation"]

	### for model 2 (factorial design with )
	### potential outcomes for I=0 D=0

	dta_sim$pot_out_00_2 <- NA
	dta_sim$pot_out_00_2[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0]
	dta_sim$pot_out_00_2[dta_sim["information"] == 1 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 0] - coef(ols_2)["information"]
	dta_sim$pot_out_00_2[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] - coef(ols_2)["information"] - coef(ols_2)["deliberation"] - coef(ols_2)["information:deliberation"]
	dta_sim$pot_out_00_2[dta_sim["information"] == 0 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 1] - coef(ols_2)["deliberation"]

	### potential outcomes for I=1 and D-1
	dta_sim$pot_out_11_2 <- NA
	dta_sim$pot_out_11_2[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1]
	dta_sim$pot_out_11_2[dta_sim["information"] == 0 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 1] +  coef(ols_2)["information"]
	dta_sim$pot_out_11_2[dta_sim["information"] == 1 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 0] + coef(ols_2)["deliberation"]
	dta_sim$pot_out_11_2[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] + coef(ols_2)["information"] + coef(ols_2)["deliberation"] + coef(ols_2)["information:deliberation"]

	### potential outcomes for I=1 and D=0
	dta_sim$pot_out_10_2 <- NA
	dta_sim$pot_out_10_2[dta_sim["information"] == 1 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 0]
	dta_sim$pot_out_10_2[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] -  coef(ols_2)["deliberation"] - coef(ols_2)["information:deliberation"]
	dta_sim$pot_out_10_2[dta_sim["information"] == 0 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 1] + coef(ols_2)["information"] -  coef(ols_2)["deliberation"]
	dta_sim$pot_out_10_2[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] +  coef(ols_2)["information"] 

	### potential outcomes for I=0 and D=1
	dta_sim$pot_out_01_2 <- NA
	dta_sim$pot_out_01_2[dta_sim["information"] == 0 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 1]
	dta_sim$pot_out_01_2[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] -  coef(ols_2)["information"] - coef(ols_2)["information:deliberation"]
	dta_sim$pot_out_01_2[dta_sim["information"] == 1 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 0] - coef(ols_2)["information"] +  coef(ols_2)["deliberation"]
	dta_sim$pot_out_01_2[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] +  coef(ols_2)["deliberation"] 

	oper <- foreach (repl = 1:nr_repl,.combine=rbind) %dopar% {
		#do the permuations	
		perm_treat <- data.frame(cbind(sample(c(rep("C",treat_nrs[1,1]), rep("D",treat_nrs[1,2]), rep("I",treat_nrs[2,1]),rep("B",treat_nrs[2,2]))),names(table(dta_sim$clusterID))))  
		names(perm_treat) <- c("perm_treat","clusterID")
		dta_perm <- merge(dta_sim, perm_treat, by.x="clusterID", by.y="clusterID")
		dta_perm$information <- ifelse(dta_perm$perm_treat %in% c("I","B"), 1, 0)
		dta_perm$deliberation <- ifelse(dta_perm$perm_treat %in% c("D","B"), 1, 0)

		dta_perm$dep_1 <- NA
		dta_perm$dep_1[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1] <- dta_perm$pot_out_1_1[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1]
		dta_perm$dep_1[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0] <- dta_perm$pot_out_0_1[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0]

		dta_perm$dep_2 <- NA
		dta_perm$dep_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1] <- dta_perm$pot_out_11_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1]
		dta_perm$dep_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0] <- dta_perm$pot_out_00_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0] 
		dta_perm$dep_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==0] <- dta_perm$pot_out_10_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==0] 
		dta_perm$dep_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==1] <- dta_perm$pot_out_01_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==1] 

### p-value
		exceed1 <- abs(coef(lm(formula1, data=dta_perm))["information:deliberation"]) > abs(coef(ols_1)["information:deliberation"])
		exceed2 <- abs(coef(lm(formula2, data=dta_perm))["information"]) > abs(coef(ols_2)["information"])
		exceed3 <- abs(coef(lm(formula2, data=dta_perm))["deliberation"]) > abs(coef(ols_2)["deliberation"])


		dta_perm[outcomes[i]] <- dta_perm$dep_1
		r1 <-coef(lm(formula1, data=dta_perm[dta_perm$deliberation == dta_perm$information,]))["information:deliberation"]

		dta_perm[outcomes[i]] <- dta_perm$dep_2
		r2 <-coef(lm(formula2, data=dta_perm))["information"]
		r3 <- coef(lm(formula2, data=dta_perm))["deliberation"]
		oper <- return(c(r1,r2,r3, exceed1, exceed2, exceed3))
	}
	return(list(conf_1 = quantile(oper[,1],sig, na.rm=T),conf_2 = quantile(oper[,2],sig, na.rm=T),conf_3 = quantile(oper[,3],sig, na.rm=T), pval_1= (sum(oper[,4])/nr_repl), pval_2= (sum(oper[,5])/nr_repl), pval_3= (sum(oper[,6])/nr_repl)))
	}


RI_conf_dist <- function(i,outcomes, baseline_outcomes, dta_sim , ctrls = NULL, nr_repl = 1000, sig = c(.025,.975)) {
#RI_conf_dist(2,outcomes, baseline_outcomes, subset(dta, ((information == 1 & deliberation==1) | district_baraza == 1)) , ctrls = "a21", nr_repl = 1000, sig = c(.025,.975))
#dta_sim <- subset(dta, ((information == 1 & deliberation==1) | district_baraza == 1))
#ctrls <- "a21"
#nr_repl <- 1000
#sig <-  c(.025,.975)

### a function to esimate confidence intervals using randomization inference following Gerber and Green pages 66-71 and 83.
	if (is.null(baseline_outcomes)) {
		formula <- as.formula(paste(outcomes[i],paste("district_baraza",ctrls,sep="+"),sep="~"))
		
	} else {
		formula <- as.formula(paste(paste(outcomes[i],paste("district_baraza",ctrls,sep="+"),sep="~"),baseline_outcomes[i],sep="+"))
	}


	dta_sim <- dta_sim %>%  mutate(clusterID = group_indices(., district))
	### get ATEs for two different models
	ols <- lm(formula, data=dta_sim) 
	
	dta_sim$dep <- as.numeric(unlist(dta_sim[as.character(formula[[2]])]))

	treat_nrs <- table(data.frame(aggregate(dta_sim["district_baraza"], list(dta_sim$clusterID),mean)[,2]))

	### calculate potential outcomes
	### for model 1 (sc effect)
	dta_sim$pot_out_0 <- NA
	dta_sim$pot_out_0[dta_sim["district_baraza"] == 0 ] <- dta_sim$dep[dta_sim["district_baraza"] == 0 ]
	dta_sim$pot_out_0[dta_sim["district_baraza"] == 1 ] <- dta_sim$dep[dta_sim["district_baraza"] == 1] - coef(ols)["district_baraza"]

	dta_sim$pot_out_1 <- NA
	dta_sim$pot_out_1[dta_sim["district_baraza"] == 0 ] <- dta_sim$dep[dta_sim["district_baraza"] == 0 ] + coef(ols)["district_baraza"]
	dta_sim$pot_out_1[dta_sim["district_baraza"] == 1 ] <- dta_sim$dep[dta_sim["district_baraza"] == 1]	

	
	oper <- foreach (repl = 1:nr_repl,.combine=rbind) %dopar% {
		#do the permuations	
		perm_treat <- data.frame(cbind(sample(c(rep("SC",treat_nrs[1]), rep("D",treat_nrs[2]))),names(table(dta_sim$clusterID))))  
		names(perm_treat) <- c("perm_treat","clusterID")
		dta_perm <- merge(dta_sim, perm_treat, by.x="clusterID", by.y="clusterID")
		dta_perm$district_baraza <- ifelse(dta_perm$perm_treat == "D", 1, 0)
	

		dta_perm$dep <- NA
		dta_perm$dep[dta_perm["district_baraza"] ==1 ] <- dta_perm$pot_out_1[dta_perm["district_baraza"] ==1 ]
		dta_perm$dep[dta_perm["district_baraza"] ==0 ] <- dta_perm$pot_out_0[dta_perm["district_baraza"] ==0 ]


		### p-value
		exceed <- abs(coef(lm(formula, data=dta_perm))["district_baraza"]) > abs(coef(ols)["district_baraza"])

		

		dta_perm[outcomes[i]] <- dta_perm$dep
		res_list <- cbind(coef(lm(formula, data=dta_perm))["district_baraza"],"exceed" = exceed)
		return(res_list) 
	}
	return(list(conf = quantile(oper[,1],sig, na.rm=T),pval= (sum(oper[,2])/nr_repl)))
}
################################################################## end of funtions declarations

#### for the mock report, I use a dummy endline - I read in a dummy endline of 3 households just to get the correct variable names
#endline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/endline.csv")[10:403]

#### I then merge with the sampling list to basically create an empty endline dataset
#### and merge in the treatments
#list <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/sampling_list_hh.csv")[c("hhid","a21","a22","a23")]
#endline <- merge(list, endline, by="hhid", all.x=T)
treats <- read.csv(paste(path,"questionnaire/final_list_5.csv", sep ="/"))
#endline <- merge(treats, endline, by.x=c("district","subcounty"), by.y=c("a22","a23"))



## baseline not needed in this first section, but used to generate fake data
baseline <- read.csv(paste(path,"data/public/baseline.csv",sep="/"))
baseline$a23 <- as.character(baseline$a23)
baseline$a23[baseline$a23 == "LUWERO TC"] <- "LUWERO_TC"
baseline$a23[baseline$a23 == "SEMBABULE TC"] <- "SEMBABULE_TC"
baseline$a23[baseline$a23 == "RAKAI TC"] <- "RAKAI_TC"
baseline$a23[baseline$a23 == "NTUSI"] <- "NTUUSI"


baseline$b21 <-  as.numeric(baseline$b21=="Yes")
baseline$b31 <-  as.numeric(baseline$b31=="Yes")
baseline$b44 <-  as.numeric(baseline$b44=="Yes")
baseline$b44[is.na(baseline$b44)] <- 0
baseline$base_inputs <- as.numeric(baseline$used_seed=="Yes" | baseline$used_fert=="Yes")
baseline$b5144 <- as.numeric(baseline$b5144=="Yes")
baseline$b5146 <- as.numeric(baseline$b5146=="Yes")
##use of unprotected water sources in dry season
###this was changed post registration to follow https://www.who.int/water_sanitation_health/monitoring/jmp2012/key_terms/en/ guidelines on what is considered improved, that also considers rainwater a protected source
baseline$base_unprotected <- as.numeric(( baseline$c11a %in%  c("Surface water","Bottled water","Cart with small tank","Unprotected dug well","Unprotected spring","Tanker truck"))    )
### is there are water committee
baseline$c10 <- as.numeric(baseline$c10=="Yes")

baseline$dist_source <- baseline$c12source 
baseline$c12source <- log(baseline$c12source + sqrt(baseline$c12source ^ 2 + 1))
baseline <- trim("c12source", baseline)
baseline <- trim("dist_source", baseline)
baseline$wait_source <- baseline$qc15
baseline$qc15 <- log(baseline$qc15 + sqrt(baseline$qc15 ^ 2 + 1))
baseline <- trim("qc15", baseline)
baseline <- trim("wait_source", baseline)
baseline$dist_road <- baseline$a6 
baseline$a6 <- log(baseline$a6 + sqrt(baseline$a6 ^ 2 + 1))
baseline <- trim("a6", baseline)
baseline <- trim("dist_road", baseline)
baseline$pub_health_access <- as.numeric((baseline$feverd21_fever %in% c("HCII","HCIII","HCIV","Regional_referral_hospital")))
baseline$maternal_health_access <- as.numeric((baseline$delivery_birthd21_delivery_birth %in% c("HCII","HCIII","HCIV","Regional_referral_hospital")))
baseline$d31 <- as.numeric(baseline$d31=="Yes")
baseline$d43 <- NA
baseline$d43[!is.na(baseline$d43a)] <- baseline$d43a[!is.na(baseline$d43a)]
baseline$d43[!is.na(baseline$d43b)] <- baseline$d43b[!is.na(baseline$d43b)] 
baseline$dist_health <- baseline$d43
baseline <- trim("dist_health", baseline)
baseline$d43 <- log(baseline$d43 + sqrt(baseline$d43 ^ 2 + 1))
baseline <- trim("d43", baseline)
## has anyone been sick
baseline$d11 <- as.numeric(baseline$d11=="Yes")

baseline$tot_sick[baseline$d11==0] <- 0 

baseline$tot_sick <- log(baseline$tot_sick + sqrt(baseline$tot_sick ^ 2 + 1))
baseline <- trim("tot_sick", baseline)

baseline$wait_time <- baseline$d410hh*60+ baseline$d410mm
baseline$wait_time_min  <- baseline$d410hh*60+ baseline$d410mm
baseline$wait_time <- log(baseline$wait_time + sqrt(baseline$wait_time ^ 2 + 1))
baseline <- trim("wait_time", baseline)
baseline <- trim("wait_time_min", baseline)

baseline$d61 <- as.numeric(baseline$d61=="Yes")

#children in public schools
baseline$base_n_children <- rowSums(cbind(baseline$e2bupe,baseline$e2aupe,baseline$e2fuse,baseline$e2muse), na.rm=T)
baseline$e5 <- rowMeans(cbind(as.numeric(baseline$e5upe) , as.numeric(baseline$e5use)), na.rm=T) 
baseline$dist_school <- rowMeans(cbind(as.numeric(baseline$e5upe) , as.numeric(baseline$e5use)), na.rm=T) 
baseline$e5[is.na(baseline$e5upe) & is.na(baseline$e5use)] <- NA
baseline$e5 <- log(baseline$e5 + sqrt(baseline$e5 ^ 2 + 1))
baseline <- trim("e5", baseline)
baseline <- trim("dist_school", baseline)


baseline$e12 <- rowSums(cbind(as.numeric(baseline$e12upe == "Yes") , as.numeric(baseline$e12use == "Yes")), na.rm=T) > 0
baseline$e12[is.na(baseline$e12upe) & is.na(baseline$e12use)] <- NA

baseline$e14 <- rowSums(cbind(as.numeric(baseline$e14upe == "Yes") , as.numeric(baseline$e14use == "Yes")), na.rm=T) > 0
baseline$e14[is.na(baseline$e14upe) & is.na(baseline$e14use)] <- NA

baseline$e22 <- rowSums(cbind(as.numeric(baseline$e22upe == "Yes") , as.numeric(baseline$e22use == "Yes")), na.rm=T) > 0
baseline$e22[is.na(baseline$e22upe) & is.na(baseline$e22use)] <- NA

baseline$e32 <- rowSums(cbind(as.numeric(baseline$e32upe == "Yes") , as.numeric(baseline$e32use == "Yes")), na.rm=T) > 0
baseline$e32[is.na(baseline$e32upe) & is.na(baseline$e32use)] <- NA

baseline$e45 <- rowSums(cbind(as.numeric(baseline$e45upe == "Yes") , as.numeric(baseline$e45use == "Yes")), na.rm=T) > 0
baseline$e45[is.na(baseline$e45upe) & is.na(baseline$e45use)] <- NA
#



##need to take logs?
baseline$log_farmsize <- log(baseline$farmsize + sqrt(baseline$farmsize ^ 2 + 1))
baseline$log_farmsize[is.infinite(baseline$log_farmsize)] <- NA
baseline <- trim("farmsize", baseline)

baseline$ironroof <- as.numeric(baseline$a512 =="Corrugated iron sheets")
baseline$improved_wall <- as.numeric(baseline$a513 %in% c("Mud_bricks_burnt_bricks","Concrete_blocks") )
baseline$head_sec <- as.numeric(baseline$a36) > 15

baseline$f241.LC1.election <- baseline$f241.LC1.election == "Yes"   
baseline$f241.LC3.election <- baseline$f241.LC3.election == "Yes"   
baseline$f241.LC5.election <- baseline$f241.LC5.election == "Yes"   
baseline$f241.Pesidential <- baseline$f241.Pesidential == "Yes"   
baseline$f241.Parliamentary <- baseline$f241.Parliamentary == "Yes"  
baseline$f241.Party.leader <- baseline$f241.Party.leader == "Yes"


baseline$f2301 <- as.numeric(baseline$f2301) %in% c(6,7)
baseline$f2303 <- as.numeric(baseline$f2303) %in% c(1,2,3,6,7)
baseline$f2307  <- as.numeric(baseline$f2307) %in% c(1,3,6,7) 
baseline$f2309 <- as.numeric(baseline$f2309) %in% c(1,2,3,6,7)
baseline$f2310 <- as.numeric(baseline$f2310) %in% c(1,2,3,6,7)

baseline$f21 <- baseline$f21 == "Yes"

baseline$roof <- baseline$a512 %in% c("Corrugated iron sheets", "Tiles")
baseline$wall <- baseline$a513 == "Mud_bricks_burnt_bricks"

baseline_desc <- baseline
baseline_matching <- merge(baseline,treats, by.x=c("a22","a23"), by.y=c("district","subcounty"))

#define endline variables -
#endline$baraza.B3 <- 0
endline$baraza.B3 <- endline$baraza.B3 ==1 |  endline$baraza.B3.3 ==1
#endline$inputs <- 0
endline$inputs <- as.numeric(endline$baraza.B1==1 | endline$baraza.B1.5==1) 
endline$baraza.B1 <- endline$baraza.B1==1
endline$baraza.B1.5 <- endline$baraza.B1.5==1

###this was changed post registration to follow https://www.who.int/water_sanitation_health/monitoring/jmp2012/key_terms/en/ guidelines on what is considered improved, that also considers rainwater a protected source
#baseline$base_unprotected <- as.numeric(( baseline$c11a %in%  c("Surface water","Bottled water","Cart with small tank","Unprotected dug well","Unprotected spring","Tanker truck"))    )
### is there are water committee
endline$unprotected <- (as.numeric(endline$baraza.C1 %in% c(5,7,9,10,11,12)) )

### here we simulate endline variables - remove if endline data is in
#endline$baraza.B2  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b21 == 1, na.rm=T))
endline$baraza.B2  <- endline$baraza.B2  == 1
###visits
### here we simulate endline variables - remove if endline data is in
#endline$baraza.B3 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b31 == 1, na.rm=T))
###naads in village
### here we simulate endline variables - remove if endline data is in
#endline$baraza.B4.1  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b44 == 1, na.rm=T))
endline$baraza.B4.1 <- endline$baraza.B4.1 == 1
###simulate an effect on this one
#endline$inputs <- rbinom(n=length(endline$inputs),size=1,prob=mean(baseline$base_inputs, na.rm=T)) 
###simulate an effect on this one
#endline$baraza.B5.2  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b5144 == 1, na.rm=T))
endline$baraza.B5.2 <- endline$baraza.B5.2 ==1 
###simulate an effect on this one
#endline$baraza.B5.3  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b5146 == 1, na.rm=T))
endline$baraza.B5.3 <- endline$baraza.B5.3 == 1
###simulate an effect on this one
#endline$unprotected <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$base_unprotected ==1, na.rm=T))
#endline$baraza.C1.2 <- sample(baseline$c12source[!is.na(baseline$c12source)],dim(endline)[1])  ### this needs to be inverse hypersine transformed and trimmed in final version
endline$baraza.C1.2 <-  as.numeric(as.character(endline$baraza.C1.2))
endline$baraza.C1.2[is.na(endline$baraza.C1.2)] <- 0 ## is na for households with piped water in compound -> distance set to 0
endline$baraza.C1.2[endline$baraza.C1.2 == 999] <- NA ## code for dont know is 999
endline$baraza.C1.2 <-  log(endline$baraza.C1.2 + sqrt(endline$baraza.C1.2 ^ 2 + 1))
endline <- trim("baraza.C1.2",endline)
#endline$baraza.C1.3 <- sample(baseline$qc15[!is.na(baseline$qc15)],dim(endline)[1]) ### this needs to be inverse hypersine transformed and trimmed in final version
endline$baraza.C1.3 <-  as.numeric(as.character(endline$baraza.C1.3))
endline$baraza.C1.3[is.na(endline$baraza.C1.3)] <- 0 ## is na for households with piped water in compound -> waiting time set to 0
endline$baraza.C1.3[endline$baraza.C1.3 == 999] <- NA ## code for dont know is 999
endline$baraza.C1.3 <-  log(endline$baraza.C1.3 + sqrt(endline$baraza.C1.3 ^ 2 + 1))
endline <- trim("baraza.C1.3",endline)

#endline$baraza.C2.3 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$c10 ==1, na.rm=T))
endline$baraza.C2.3 <- endline$baraza.C2.3 == 1
#endline$baraza.A6 <- sample(baseline$a6[!is.na(baseline$a6)],dim(endline)[1]) ### this needs to be inverse hypersine transformed and trimmed in final version
endline$baraza.A6 <-  as.numeric(as.character(endline$baraza.A6))
endline$baraza.A6[endline$baraza.A6 == 999] <- NA ## code for dont know is 999
endline$baraza.A6 <-  log(endline$baraza.A6 + sqrt(endline$baraza.A6 ^ 2 + 1))
endline <- trim("baraza.A6",endline)

### access to public health if sick
#endline$baraza.D2 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$pub_health_access, na.rm=T))
endline$baraza.D2 <- (endline$baraza.D2 %in% 2:5)
#endline$baraza.D2.4 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$maternal_health_access , na.rm=T))
endline$baraza.D2.4 <- (endline$baraza.D2.4 %in% 2:5)
#endline$baraza.D3 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$pub_health_access, na.rm=T))
endline$baraza.D3 <- endline$baraza.D3 == 1
#endline$baraza.D4.2 <- sample(baseline$d43[!is.na(baseline$d43)] ,dim(endline)[1])  ### this needs to be inverse hypersine transformed and trimmed in final version
endline$baraza.D4.2 <-  as.numeric(as.character(endline$baraza.D4.2))
endline$baraza.D4.2[endline$baraza.D4.2 == 999] <- NA ## code for dont know is 999
endline$baraza.D4.2 <-  log(endline$baraza.D4.2 + sqrt(endline$baraza.D4.2 ^ 2 + 1))
endline <- trim("baraza.D4.2",endline)

#health outcome - less people sick 
#endline$baraza.D1  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$d11, na.rm=T))
endline$baraza.D1 <- endline$baraza.D1 == 1

members <- paste(paste("baraza.labour", 1:15, sep="."),".D1.2", sep=".")

endline[members] <- lapply(endline[members], function(x) as.numeric(as.character(x)) )
endline[members] <- lapply(endline[members], function(x) replace(x, x == 999,NA) )
endline$baraza.D1.2 <- rowSums(endline[members], na.rm=T)

endline$baraza.D1.2 <- log(endline$baraza.D1.2 + sqrt(endline$baraza.D1.2 ^ 2 + 1))
endline <- trim("baraza.D1.2", endline)

#endline$baraza.D4.6 <- sample(baseline$wait_time[!is.na(baseline$wait_time)] ,dim(endline)[1]) 
endline$baraza.D4.6 <-  as.numeric(as.character(endline$baraza.D4.6))
endline$baraza.D4.6[endline$baraza.D4.6 == 999] <- NA ## code for dont know is 999
endline$baraza.D4.6 <-  log(endline$baraza.D4.6 + sqrt(endline$baraza.D4.6 ^ 2 + 1))
endline <- trim("baraza.D4.6",endline)

#endline$baraza.D6  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$d61, na.rm=T))
endline$baraza.D6 <- endline$baraza.D6 == 1

##edu
endline$n_children <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.1)),as.numeric(as.character(endline$baraza.E2.1))), na.rm=T) 
#endline$n_children <- sample(baseline$base_n_children ,dim(endline)[1]) 
endline$baraza.E1.2[endline$baraza.E1.2 == 999] <- NA
endline$baraza.E2.2[endline$baraza.E2.2 == 999] <- NA 
endline$baraza.E5  <- rowMeans(cbind(as.numeric(as.character(endline$baraza.E1.2)),as.numeric(as.character(endline$baraza.E2.2))), na.rm=T) 
endline$baraza.E5[is.na(as.numeric(as.character(endline$baraza.E1.2))) & is.na(as.numeric(as.character(endline$baraza.E2.2)))] <- NA
endline$baraza.E5 <- log(endline$baraza.E5 + sqrt(endline$baraza.E5 ^ 2 + 1))
#boundery fence
#endline$baraza.E1.4 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e12, na.rm=T))
endline$baraza.E12 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.4)) == 1 , as.numeric(as.character(endline$baraza.E2.4))==1), na.rm=T) > 0
endline$baraza.E12[is.na(as.numeric(as.character(endline$baraza.E1.4)) ) & is.na(as.numeric(as.character(endline$baraza.E2.4)) )] <- NA

endline$baraza.E14 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.6)) == 1 , as.numeric(as.character(endline$baraza.E2.6))==1), na.rm=T) > 0
endline$baraza.E14[is.na(as.numeric(as.character(endline$baraza.E1.6)) ) & is.na(as.numeric(as.character(endline$baraza.E2.6)) )] <- NA

endline$baraza.E22 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.10)) == 1 , as.numeric(as.character(endline$baraza.E2.10))==1), na.rm=T) > 0
endline$baraza.E22[is.na(as.numeric(as.character(endline$baraza.E1.10)) ) & is.na(as.numeric(as.character(endline$baraza.E2.10)) )] <- NA

endline$baraza.E32 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.13)) == 1 , as.numeric(as.character(endline$baraza.E2.13))==1), na.rm=T) > 0
endline$baraza.E32[is.na(as.numeric(as.character(endline$baraza.E1.13)) ) & is.na(as.numeric(as.character(endline$baraza.E2.13)) )] <- NA

endline$baraza.E45 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.18)) == 1 , as.numeric(as.character(endline$baraza.E2.18))==1), na.rm=T) > 0
endline$baraza.E45[is.na(as.numeric(as.character(endline$baraza.E1.18)) ) & is.na(as.numeric(as.character(endline$baraza.E2.18)) )] <- NA

### assorted outcomes
### type of roof
endline$baraza.roof <- endline$baraza.roof %in% 1:2

endline$baraza.wall <- endline$baraza.wall == 2

### type of wall



#PLEASE DESCRIBE YOUR PARTICIPATION IN DIFFERENT TYPES OF ELECTIONS:
#F1	F1. In this household, are there any members who currently hold any political/traditional positions?
endline$baraza.F1 <- (endline$baraza.F1 == 1)


#F2	F2. In the last elections did you participate in the  LCI elections? 
#F2.1	F2.1 In the last elections did you participate in the  LC3 elections? 
#F2.2	F2.2  In the last elections did you participate in the LC5  elections?   
#F2.3	F2.3 In the last elections did you participate in  the Presidential elections?  
#F2.4	F2.4  In the last elections did you participate in the Parliamentary election?
#F2.5	F2.5 In the last elections did you participate in the  Party leaders elections? 

endline$baraza.part.F2 <- endline$baraza.part.F2 == 1 
endline$baraza.part.F2.1 <- endline$baraza.part.F2.1 == 1 
endline$baraza.part.F2.2 <- endline$baraza.part.F2.2 == 1
endline$baraza.part.F2.3 <- endline$baraza.part.F2.3 == 1
endline$baraza.part.F2.4 <- endline$baraza.part.F2.4 == 1 
endline$baraza.part.F2.5 <- endline$baraza.part.F2.5 == 1


##ag

#1 


# 7 #make an pol index
endline <- FW_index(c("baraza.F1","baraza.part.F2","baraza.part.F2.1","baraza.part.F2.2","baraza.part.F2.3","baraza.part.F2.4","baraza.part.F2.5"),data=endline)
names(endline)[names(endline) == 'index'] <- 'pol_index'
baseline <- FW_index(c("f21","f241.LC1.election", "f241.LC3.election","f241.LC5.election", "f241.Pesidential", "f241.Parliamentary", "f241.Party.leader"  ),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_pol_index'
baseline_matching <- FW_index(c("f21","f241.LC1.election", "f241.LC3.election","f241.LC5.election", "f241.Pesidential", "f241.Parliamentary", "f241.Party.leader"  ),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_pol_index'

##F1.1 When was the last time that you spoke personally with the LC 1 Chairperson, for a reason relating to service provision in agriculture, health, education, water or roads?
##F1.2 When was the last time that you spoke personally with the Subcounty Chief, for a reason relating to service provision in agriculture, health, education, water or roads?
##F1.3 When was the last time that you spoke personally with the Head teacher/ SMC member, for a reason relating to service provision in agriculture, health, education, water or roads?
##F1.4  When was the last time that you spoke personally with the Health Unit Management Committee (HUMC) Member, for a reason relating to service provision in agriculture, health, education, water or roads?
##F1.5  When was the last time that you spoke personally with the Water committee member, for a reason relating to service provision in agriculture, health, education, water or roads?
endline$baraza.F1.1 <- (endline$baraza.F1.1<=2) # last month
endline$baraza.F1.2 <- (endline$baraza.F1.2<=5) # last year 
endline$baraza.F1.3 <- (endline$baraza.F1.3<=4) # last 1/2 year 
endline$baraza.F1.4 <- (endline$baraza.F1.4<=5) # last year 
endline$baraza.F1.5 <- (endline$baraza.F1.5<=5) # last year 

###make a contact index
endline <- FW_index(c("baraza.F1.1", "baraza.F1.2", "baraza.F1.3","baraza.F1.4","baraza.F1.5"),data=endline)
names(endline)[names(endline) == 'index'] <- 'contact_index'
baseline <- FW_index(c("f2301","f2303", "f2307","f2309","f2310"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_contact_index'
baseline_matching <- FW_index(c("f2301","f2303", "f2307","f2309","f2310"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_contact_index'

###priority rankings
#H1	H1.  Access to a drinking water source is a serious problem. (where 1 means you completely disagree and 10 means you completely agree) 
#H2	H2. Drinking water is usually dirty. (where 1 means you completely disagree and 10 means you completely agree)    
#H3	H3.  Access to a government health centre or hospital is a serious problem. (where 1 means you completely disagree and 10 means you completely agree)
#H4	H4  Government health centres or hospitals do not have relevant medicines. (where 1 means you completely disagree and 10 means you completely agree)
#H5	H5 Staff at government health centres or hospitals are rude to patients. (where 1 means you completely disagree and 10 means you completely agree)
#H6	H6.  Medical staff at government health centres or hospitals are often absent. (where 1 means you completely disagree and 10 means you completely agree)
#H7	H7 Access to a government primary school is a serious problem. (where 1 means you completely disagree and 10 means you completely agree)
#H8	H8 Teachers in government schools are often absent. (where 1 means you completely disagree and 10 means you completely agree)
#H9	H9 Children’s learning outcomes in government schools are poor. (where 1 means you completely disagree and 10 means you completely agree)  
#H10	H10.  Availability/ Access to all-weather roads is a serious problem (where 1 means you completely disagree and 10 means you completely agree)
#H11	H11  Agricultural inputs supplied by the government are of poor quality. (where 1 means you completely disagree and 10 means you completely agree)
#H11b	H11  Agricultural inputs supplied by the government are delivered in time. (where 1 means you completely disagree and 10 means you completely agree)
#H12	H12 There is lack of transparency in how farmers are selected to receive agricultural inputs from the government. (where 1 means you completely disagree and 10 means you completely agree)
#H13	H13  Agricultural extension agents rarely visit. (where 1 means you completely disagree and 10 means you completely agree)
#H14	H14.  Agricultural extension agents are not aware of enterprises or agricultural inputs relevant to farmers. (where 1 means you completely disagree and 10 means you completely agree)

endline$baraza.H1[endline$baraza.H1 == 11] <- NA
endline$baraza.H2[endline$baraza.H2 == 11] <- NA
endline$baraza.H3[endline$baraza.H3 == 11] <- NA
endline$baraza.H4[endline$baraza.H4 == 11] <- NA
endline$baraza.H5[endline$baraza.H5 == 11] <- NA
endline$baraza.H6[endline$baraza.H6 == 11] <- NA
endline$baraza.H7[endline$baraza.H7 == 11] <- NA
endline$baraza.H8[endline$baraza.H8 == 11] <- NA
endline$baraza.H9[endline$baraza.H9 == 11] <- NA
endline$baraza.H10[endline$baraza.H10 == 11] <- NA
endline$baraza.H11[endline$baraza.H11 == 11] <- NA
endline$baraza.H12[endline$baraza.H12 == 11] <- NA
endline$baraza.H13[endline$baraza.H13 == 11] <- NA
endline$baraza.H14[endline$baraza.H14 == 11] <- NA 

endline <- FW_index(c("baraza.H1","baraza.H2","baraza.H3","baraza.H4", "baraza.H5",  "baraza.H6", "baraza.H7", "baraza.H8", "baraza.H9", "baraza.H10", "baraza.H11", "baraza.H12", "baraza.H13","baraza.H14"),data=endline)
names(endline)[names(endline) == 'index'] <- 'priority_index'
baseline <- FW_index(c("i1","i2","i3","i4","i5","i6","i7","i8","i9","i10","i11","i12","i13","i14"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_priority_index'
baseline_matching <- FW_index(c("i1","i2","i3","i4","i5","i6","i7","i8","i9","i10","i11","i12","i13","i14"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_priority_index'

### contributions
#G1	G1: Did you ever make any contributions to the school in the last two years? 
#G1.1	G1.1 Did you ever make any contributions to the health centre in the last two years? 
#G1.2	G1.2 Did you ever make any contributions to the road/ bridge in the last two years? 
#G1.3	G1.3 Did you ever make any contributions to the drinking water facility in the last two years? 
#G1.4	G1.4 Did you ever make any contributions to the dam/ irrigation facility in the last two years? 
#G1.5	G1.5 Did you  make any contributions to any other government or community building/ structure in the last two years? 
endline$baraza.G1k <- endline$baraza.G1 %in% c(1,3)
endline$baraza.G1c <- endline$baraza.G1 %in% c(2,3)

endline$baraza.G1.1k <- endline$baraza.G1.1 %in% c(1,3)
endline$baraza.G1.1c <- endline$baraza.G1.1 %in% c(2,3)

endline$baraza.G1.2k <- endline$baraza.G1.2 %in% c(1,3)
endline$baraza.G1.2c <- endline$baraza.G1.2 %in% c(2,3)

endline$baraza.G1.3k <- endline$baraza.G1.3 %in% c(1,3)
endline$baraza.G1.3c <- endline$baraza.G1.3 %in% c(2,3)

endline$baraza.G1.4k <- endline$baraza.G1.4 %in% c(1,3)
endline$baraza.G1.4c <- endline$baraza.G1.4 %in% c(2,3)

endline$baraza.G1.5k <- endline$baraza.G1.5 %in% c(1,3)
endline$baraza.G1.5c <- endline$baraza.G1.5 %in% c(2,3)

####make a contribution  index
## in kind
endline <- FW_index(c("baraza.G1k","baraza.G1.1k","baraza.G1.2k","baraza.G1.3k","baraza.G1.4k","baraza.G1.5k"),data=endline)
names(endline)[names(endline) == 'index'] <- 'in_kind_index'
baseline <- FW_index(c("cschoolk", "chealthk","croadk","cwaterk",  "cdamk",    "cbuildk"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_in_kind_index'
baseline_matching <- FW_index(c("cschoolk", "chealthk","croadk","cwaterk",  "cdamk",    "cbuildk"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_in_kind_index'
## in cash
endline <- FW_index(c("baraza.G1c","baraza.G1.1c","baraza.G1.2c","baraza.G1.3c","baraza.G1.4c","baraza.G1.5c"),data=endline)
names(endline)[names(endline) == 'index'] <- 'in_cash_index'
baseline <- FW_index(c("cschoolc", "chealthc","croadc","cwaterc",  "cdamc",    "cbuildc"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_in_cash_index'
baseline_matching <- FW_index(c("cschoolc", "chealthc","croadc","cwaterc",  "cdamc",    "cbuildc"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_in_cash_index'


outcomes <- c("baraza.F1","baraza.part.F2","baraza.part.F2.1","baraza.part.F2.2","baraza.part.F2.3","baraza.part.F2.4","baraza.part.F2.5", "pol_index", "baraza.F1.1", "baraza.F1.2", "baraza.F1.3","baraza.F1.4","baraza.F1.5", "contact_index","baraza.H1","baraza.H2","baraza.H3","baraza.H4", "baraza.H5",  "baraza.H6", "baraza.H7", "baraza.H8", "baraza.H9", "baraza.H10", "baraza.H11", "baraza.H12", "baraza.H13","baraza.H14","priority_index","baraza.G1k","baraza.G1.1k","baraza.G1.2k","baraza.G1.3k","baraza.G1.4k","baraza.G1.5k","in_kind_index","baraza.G1c", "baraza.G1.1c","baraza.G1.2c","baraza.G1.3c","baraza.G1.4c","baraza.G1.5c", "in_cash_index" )
baseline_outcomes <- c("f21","f241.LC1.election", "f241.LC3.election","f241.LC5.election", "f241.Pesidential", "f241.Parliamentary", "f241.Party.leader" ,"base_pol_index","f2301","f2303", "f2307","f2309","f2310","base_contact_index","i1","i2","i3","i4","i5","i6","i7","i8","i9","i10","i11","i12","i13","i14","base_priority_index","cschoolk", "chealthk","croadk","cwaterk",  "cdamk",    "cbuildk","base_in_kind_index","cschoolc", "chealthc","croadc","cwaterc",  "cdamc",  "cbuildc","base_in_cash_index")

#create unique ID for clustering based on district and subcounty
endline <- endline %>%  mutate(clusterID = group_indices(., district, subcounty))
endline <- endline %>%  mutate(clusterID2 = group_indices(., district))
##dataset to be used for ancova
dta <- merge(endline, baseline[, -which(names(baseline)=="a21")], by="hhid")
endline$time <- 1
baseline$time <- 0
baseline$information <- 0
baseline$deliberation <- 0
baseline$district_baraza <- 0
baseline$information[baseline$treat=="info" | baseline$treat=="scbza"] <- 1 
baseline$deliberation[baseline$treat=="delib" | baseline$treat=="scbza"] <- 1 
baseline$district_baraza[baseline$treat=="dbza"] <- 1 


### merge in clusterID for standard error clustering in dif-in-dif
baseline <- merge(baseline, endline[c("hhid","clusterID","clusterID2")], by="hhid", all.y=T)
baseline_matching <- merge(baseline_matching, endline[c("hhid","clusterID","clusterID2")], by="hhid", all.y=T)

baseline <- baseline[c("information","deliberation","district_baraza","time","clusterID","clusterID2",baseline_outcomes )]
names(baseline) <- c("information","deliberation","district_baraza","time","clusterID","clusterID2",outcomes )


dta_long <- rbind(endline[c("information","deliberation","district_baraza","time", "clusterID","clusterID2",outcomes)], baseline[ c("information","deliberation","district_baraza","time", "clusterID","clusterID2",outcomes )])

### uncomment for heterogeneity re: time
if (hetero == 1) {
dta <- subset(dta, (time_dif>1.5) | time_dif == 0)
}
if (hetero == 2) {
dta <- subset(dta, j9 >= 5 )
}


###init arrays to store results
df_ancova <- array(NA,dim=c(6,5,length(outcomes)))
df_averages <- array(NA,dim=c(2,length(outcomes)))

for (i in 1:length(outcomes)) {
print(i)
# i <- 1

df_averages[1,i] <- mean(as.matrix(endline[outcomes[i]]), na.rm=T)
df_averages[2,i] <- sd(as.matrix(endline[outcomes[i]]), na.rm=T)


##ancova
## merge in baseline

ols <- lm(as.formula(paste(paste(outcomes[i],"information*deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta[dta$district_baraza == 0,]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID[dta$district_baraza == 0], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes, baseline_outcomes, subset(dta, district_baraza == 0) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <- RI_store$conf_2 
conf[3,4:5] <- RI_store$conf_3
res[2,5] <- RI_store$pval_2
res[3,5] <- RI_store$pval_3
}

df_ancova[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
df_ancova[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5], nobs(ols))


ols <- lm(as.formula(paste(paste(outcomes[i],"information:deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta[dta$district_baraza == 0 & (dta$information == dta$deliberation),]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID[dta$district_baraza == 0 & (dta$information == dta$deliberation)], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
conf[6,4:5] <- RI_store$conf_1
res[6,5] <- RI_store$pval_1
}

df_ancova[,1,i] <- c(res[6,1],res[6,2],res[6,5], conf[6,4],conf[6,5], nobs(ols))

ols <- lm(as.formula(paste(paste(outcomes[i],"district_baraza+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta[(dta$information == 0 & dta$deliberation==0) | dta$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID2[(dta$information == 0 & dta$deliberation==0) | dta$district_baraza == 1 ], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_dist(i,outcomes, baseline_outcomes, subset(dta, ((information == 0 & deliberation==0) | district_baraza == 1)) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <-  RI_store$conf
res[2,5] <- RI_store$pval
}
df_ancova[,4,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))


ols <- lm(as.formula(paste(paste(outcomes[i],"district_baraza+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta[(dta$information == 1 & dta$deliberation==1) | dta$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID2[(dta$information == 1 & dta$deliberation==1) | dta$district_baraza == 1 ], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_dist(i,outcomes, baseline_outcomes, subset(dta, ((information == 1 & deliberation==1) | district_baraza == 1)) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <-  RI_store$conf
res[2,5] <- RI_store$pval
}
df_ancova[,5,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
}


### save results
save_path <- ifelse(final_verion_swith, paste(path,"report/results/final", sep = "/"), paste(path,"report/results/", sep = "/"))
save_path <- ifelse(hetero ==1, paste(save_path,"hetero1", sep = "/"),  ifelse(hetero ==2, paste(save_path,"hetero2", sep = "/"), save_path))

save(df_ancova, file= paste(save_path,"df_ancova_additional.Rd", sep="/"))
save(df_averages, file= paste(save_path,"df_averages_additional.Rd", sep="/"))
save(baseline_desc, file= paste(save_path,"baseline_desc_additonal.Rd", sep="/"))



