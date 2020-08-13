#### run this file for the pre-registered (confirmatory) analysis
#### main output are the graphs

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

### set this switch to TRUE if you want to produce a final report - this will save results matrices in a static directory
final_verion_swith <- TRUE

RI_conf_switch <- TRUE
glob_repli <- 10000
glob_sig <- c(.025,.975) ### 5 percent conf intervals


path <- strsplit(getwd(), "/report")[[1]]
endline <- read.csv(paste(path,"data/public/endline.csv", sep="/"), stringsAsFactors = TRUE)
endline$a21 <- as.character(endline$region)
endline$region <- NULL
endline$interviewed <- TRUE
dta_plan <- read.csv(paste(path,"questionnaire/sampling_list_hh.csv", sep="/"), stringsAsFactors = TRUE)
endline <- merge(dta_plan, endline, by="hhid", all.x=T)
endline$interviewed[is.na(endline$interviewed)] <- FALSE
endline$a21 <- endline$a21.x
endline$a21.y <- NULL
endline$information <- NULL
endline$deliberation <- NULL
endline$district_baraza <- NULL

endline$attriter <- !endline$interviewed 



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

## function definitions 
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
#endline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/endline.csv", stringsAsFactors = TRUE)[10:403]

treats <- read.csv(paste(path,"questionnaire/final_list_5.csv", sep ="/"), stringsAsFactors = TRUE)
endline <- merge(treats, endline, by.x=c("district","subcounty"), by.y=c("a22","a23"))

#create unique ID for clustering based on district and subcounty
endline <- endline %>%  mutate(clusterID = group_indices(., district, subcounty))
endline <- endline %>%  mutate(clusterID2 = group_indices(., district))

outcomes <- "attriter"
baseline_outcomes <- NULL

###init arrays to store results
df_ancova <- array(NA,dim=c(6,5,length(outcomes)))
df_averages <- array(NA,dim=c(2,length(outcomes)))
### parallel computing for RI
cl <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
registerDoParallel(cl)
dta <- endline
for (i in 1:length(outcomes)) {
#print(i)
# i <- 1

df_averages[1,i] <- mean(as.matrix(endline[outcomes[i]]), na.rm=T)
df_averages[2,i] <- sd(as.matrix(endline[outcomes[i]]), na.rm=T)


##ancova
## merge in baseline

ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+a21",sep="~")), data=dta[dta$district_baraza == 0,]) 
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


ols <- lm(as.formula(paste(outcomes[i],"information:deliberation+a21",sep="~")), data=dta[dta$district_baraza == 0 & (dta$information == dta$deliberation),]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID[dta$district_baraza == 0 & (dta$information == dta$deliberation)], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
conf[5,4:5] <- RI_store$conf_1
res[5,5] <- RI_store$pval_1
}

df_ancova[,1,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5], nobs(ols))

ols <- lm(as.formula(paste(outcomes[i],"district_baraza+a21",sep="~"),), data=dta[(dta$information == 0 & dta$deliberation==0) | dta$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID2[(dta$information == 0 & dta$deliberation==0) | dta$district_baraza == 1 ], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_dist(i,outcomes, baseline_outcomes, subset(dta, ((information == 0 & deliberation==0) | district_baraza == 1)) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <-  RI_store$conf
res[2,5] <- RI_store$pval
}
df_ancova[,4,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))


ols <- lm(as.formula(paste(outcomes[i],"district_baraza+a21",sep="~")), data=dta[(dta$information == 1 & dta$deliberation==1) | dta$district_baraza == 1 ,]) 
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

save(df_ancova, file= paste(save_path,"df_ancova_attrition.Rd", sep="/"))
save(df_averages, file= paste(save_path,"df_averages_attrition.Rd", sep="/"))


