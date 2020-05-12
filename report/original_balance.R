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
set.seed(12345) #not needed for final version?

### this is executed in the /report subdirectory, need to ..
path <- strsplit(getwd(), "/report")[[1]]

### set this switch to TRUE if you want to produce a final report - this will save results matrices in a static directory
final_verion_swith <- FALSE
RI_conf_switch <- FALSE
glob_repli <- 100
glob_sig <- c(.025,.975) ### 5 percent conf intervals

########################################################### functions declarations #####################################################

trim <- function(var, dataset, trim_perc=.05) {
### function for triming a variable in a dataset - replaces with NA
 dataset[var][dataset[var] < quantile(dataset[var],c(trim_perc/2,1-(trim_perc/2)), na.rm=T)[1] | dataset[var] > quantile(dataset[var], c(trim_perc/2,1-(trim_perc/2)),na.rm=T)[2] ] <- NA
return(dataset)
}

FW_index <- function(indexer,revcols = NULL,data_orig) {
### function to make family wise index using covariance as weights (following http://cyrussamii.com/?p=2656)
### FW_index("messenger != 'ctrl' ", c("know_space", "know_combine", "know_weed"),dta)

#FW_index(c("baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1",  "baraza.D1.2"),revcols=c(5,6),data=endline)
#indexer <- c("baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1",  "baraza.D1.2")
#revcols <- c(5,6)
#data_orig <- endline

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


	dta_sim <- dta_sim %>%  mutate(clusterID = group_indices(., a22, a23))
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
	return(list(conf_1 = quantile(oper[,1],sig),conf_2 = quantile(oper[,2],sig),conf_3 = quantile(oper[,3],sig), pval_1= (sum(oper[,4])/nr_repl), pval_2= (sum(oper[,5])/nr_repl), pval_3= (sum(oper[,6])/nr_repl)))
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


	dta_sim <- dta_sim %>%  mutate(clusterID = group_indices(., a22))
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
	return(list(conf = quantile(oper[,1],sig),pval= sum(oper[,2])/nr_repl))
}

################################################################## end of funtions declarations

treats <- read.csv(paste(path,"questionnaire/final_list_5.csv", sep ="/"))

## baseline not needed in this first section, but used to generate fake data
baseline <- read.csv(paste(path,"data/public/baseline.csv", sep ="/"))
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
baseline$base_unprotected <- (as.numeric(baseline$c11a) %in%  c(10,13,14))
### is there are water committee
baseline$c10 <- as.numeric(baseline$c10=="Yes")
baseline$c12source <- log(baseline$c12source + sqrt(baseline$c12source ^ 2 + 1))
baseline <- trim("c12source", baseline)
baseline$qc15 <- log(baseline$qc15 + sqrt(baseline$qc15 ^ 2 + 1))
baseline <- trim("qc15", baseline)
baseline$a6 <- log(baseline$a6 + sqrt(baseline$a6 ^ 2 + 1))
baseline <- trim("a6", baseline)
baseline$pub_health_access <- as.numeric((baseline$feverd21_fever %in% c("HCII","HCIII","HCIV","Regional_referral_hospital")))
baseline$maternal_health_access <- as.numeric((baseline$delivery_birthd21_delivery_birth %in% c("HCII","HCIII","HCIV","Regional_referral_hospital")))
baseline$d31 <- as.numeric(baseline$d31=="Yes")
baseline$d43 <- NA
baseline$d43[!is.na(baseline$d43a)] <- baseline$d43a[!is.na(baseline$d43a)]
baseline$d43[!is.na(baseline$d43b)] <- baseline$d43b[!is.na(baseline$d43b)] 
baseline$d43 <- log(baseline$d43 + sqrt(baseline$d43 ^ 2 + 1))
baseline <- trim("d43", baseline)
## has anyone been sick
baseline$d11 <- as.numeric(baseline$d11=="Yes")

baseline$tot_sick[baseline$d11==0] <- 0 

#children in public schools
baseline$base_n_children <- rowSums(cbind(baseline$e2bupe,baseline$e2aupe,baseline$e2fuse,baseline$e2muse), na.rm=T)
baseline$e5 <- rowMeans(cbind(as.numeric(baseline$e5upe) , as.numeric(baseline$e5use)), na.rm=T) 
baseline$e5[is.na(baseline$e5upe) & is.na(baseline$e5use)] <- NA

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

baseline_matching <- merge(baseline,treats, by.x=c("a22","a23"), by.y=c("district","subcounty"))

baseline$information <- 0
baseline$deliberation <- 0
baseline$district_baraza <- 0
baseline$information[baseline$treat=="info" | baseline$treat=="scbza"] <- 1 
baseline$deliberation[baseline$treat=="delib" | baseline$treat=="scbza"] <- 1 
baseline$district_baraza[baseline$treat=="dbza"] <- 1 

baseline$thatched <- (baseline$a512 == "Grass leaf thatched")
baseline$trad_wall <- (baseline$a513 == "Wood_and_mud")


outcomes <- c("hhsize","agehead","femhead","head_sec","thatched","trad_wall","a6","b21","d31","base_n_children")




#create unique ID for clustering based on district and subcounty
baseline <- baseline %>%  mutate(clusterID = group_indices(., a22, a23))
baseline <- baseline %>%  mutate(clusterID2 = group_indices(., a22))




###init arrays to store results
df_ols <- array(NA,dim=c(6,5,length(outcomes)))
df_averages <- array(NA,dim=c(2,length(outcomes)))

cl <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
registerDoParallel(cl)

for (i in 1:length(outcomes)) {
print(i)

df_averages[1,i] <- mean(as.matrix(baseline[outcomes[i]]), na.rm=T)
df_averages[2,i] <- sd(as.matrix(baseline[outcomes[i]]), na.rm=T)

### simple difference and adjust se for clustered treatment assignment

ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+a21",sep="~")), data=baseline[baseline$district_baraza == 0,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID[baseline$district_baraza == 0], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes,NULL, dta_sim = subset(baseline, district_baraza == 0) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <- RI_store$conf_2 
conf[3,4:5] <- RI_store$conf_3
res[2,5] <- RI_store$pval_2
res[3,5] <- RI_store$pval_3
}

df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5], nobs(ols))

ols <- lm(as.formula(paste(outcomes[i],"information:deliberation+a21",sep="~")), data=baseline[baseline$district_baraza == 0 & (baseline$information == baseline$deliberation),]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID[baseline$district_baraza == 0 & (baseline$information == baseline$deliberation)], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
conf[5,4:5] <- RI_store$conf_1
res[5,5] <- RI_store$pval_1
}

df_ols[,1,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5], nobs(ols))

ols <- lm(as.formula(paste(outcomes[i],"district_baraza+a21",sep="~")), data=baseline[(baseline$information == 0 & baseline$deliberation==0) | baseline$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID2[(baseline$information == 0 & baseline$deliberation==0) | baseline$district_baraza == 1 ], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_dist(i,outcomes, NULL, subset(baseline, ((information == 0 & deliberation==0) | district_baraza == 1)) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <-  RI_store$conf
res[2,5] <- RI_store$pval
}
df_ols[,4,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))

ols <- lm(as.formula(paste(outcomes[i],"district_baraza+a21",sep="~")), data=baseline[(baseline$information == 1 & baseline$deliberation==1) | baseline$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID2[(baseline$information == 1 & baseline$deliberation==1) | baseline$district_baraza == 1 ], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_dist(i,outcomes, NULL, subset(baseline, ((information == 1 & deliberation==1) | district_baraza == 1)) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <-  RI_store$conf
res[2,5] <- RI_store$pval
}
df_ols[,5,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
}

### merge in treatments to drop the treated and then compare planned controls and untreated sub-counties
#for information treatment
baseline_treat <- merge(baseline, treats, all.x=T, by.y=c("district","subcounty"), by.x=c("a22","a23"))

baseline_treat <- subset(baseline_treat,district_baraza.x == 0)
### drop if information.y==1
baseline_treat$information.y[is.na(baseline_treat$information.y)] <- 0
baseline_treat_info <- subset(baseline_treat, information.y!=1)

##init arrays to store results
df_balance <- array(NA,dim=c(6,3,length(outcomes)))

baseline_treat_info$information <- baseline_treat_info$information.x
baseline_treat_info$deliberation <- baseline_treat_info$deliberation.x


for (i in 1:length(outcomes)) {
## simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+a21",sep="~")), data=baseline_treat_info) 
vcov_cluster <- vcovCR(ols, cluster = baseline_treat_info$clusterID, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes,NULL, dta_sim = baseline_treat_info , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <- RI_store$conf_2 
res[2,5] <- RI_store$pval_2
}


df_balance[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))

}

#for deliberation
baseline_treat <- merge(baseline, treats, all.x=T, by.y=c("district","subcounty"), by.x=c("a22","a23"))

baseline_treat <- subset(baseline_treat,district_baraza.x == 0)
baseline_treat$deliberation.y[is.na(baseline_treat$deliberation.y)] <- 0
baseline_treat_delib <- subset(baseline_treat, deliberation.y!=1)


baseline_treat_delib$information <- baseline_treat_delib$information.x
baseline_treat_delib$deliberation <- baseline_treat_delib$deliberation.x

for (i in 1:length(outcomes)) {
### simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+a21",sep="~")), data=baseline_treat_delib) 
vcov_cluster <- vcovCR(ols, cluster = baseline_treat_delib$clusterID, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes,NULL, dta_sim = baseline_treat_delib , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[3,4:5] <- RI_store$conf_3 
res[3,5] <- RI_store$pval_3
}

df_balance[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5], nobs(ols))
}

#for interaction
baseline_treat <- merge(baseline, treats, all.x=T, by.y=c("district","subcounty"), by.x=c("a22","a23"))

baseline_treat <- subset(baseline_treat,district_baraza.x == 0)
### drop if information.y==1 & information.y==1
baseline_treat$information.y[is.na(baseline_treat$information.y)] <- 0
baseline_treat_info <- subset(baseline_treat, information.y!=1 | deliberation.y!=1)

baseline_treat_info$information <- baseline_treat_info$information.x
baseline_treat_info$deliberation <- baseline_treat_info$deliberation.x

for (i in 1:length(outcomes)) {
## simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information.x:deliberation.x+a21",sep="~")), data=baseline_treat_info[baseline_treat_info$information.x ==  baseline_treat_info$deliberation.x, ]) 
vcov_cluster <- vcovCR(ols, cluster = baseline_treat_info$clusterID[baseline_treat_info$information.x ==  baseline_treat_info$deliberation.x], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes,NULL, dta_sim = baseline_treat_info , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[5,4:5] <- RI_store$conf_1 
res[5,5] <- RI_store$pval_1
}

df_balance[,1,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5], nobs(ols))

}


### redo balance table, but now between actual treated and matched controls
baseline$information <- NULL
baseline$deliberation <- NULL
baseline$district_baraza <- NULL


baseline <- merge(baseline, treats, all.y=T, by.y=c("district","subcounty"), by.x=c("a22","a23"))

###init arrays to store results
df_ols_end <- array(NA,dim=c(6,5,length(outcomes)))
df_averages_end <- array(NA,dim=c(2,length(outcomes)))


for (i in 1:length(outcomes)) {
print(i)

df_averages_end[1,i] <- mean(as.matrix(baseline[outcomes[i]]), na.rm=T)
df_averages_end[2,i] <- sd(as.matrix(baseline[outcomes[i]]), na.rm=T)

### simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+a21",sep="~")), data=baseline[baseline$district_baraza == 0,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID[baseline$district_baraza == 0], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes,NULL, dta_sim = subset(baseline, district_baraza == 0) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <- RI_store$conf_2 
conf[3,4:5] <- RI_store$conf_3
res[2,5] <- RI_store$pval_2
res[3,5] <- RI_store$pval_3
}

df_ols_end[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
df_ols_end[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5], nobs(ols))

ols <- lm(as.formula(paste(outcomes[i],"information:deliberation+a21",sep="~")), data=baseline[baseline$district_baraza == 0 & (baseline$information == baseline$deliberation),]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID[baseline$district_baraza == 0 & (baseline$information == baseline$deliberation)], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
conf[5,4:5] <- RI_store$conf_1
res[5,5] <- RI_store$pval_1
}


df_ols_end[,1,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5], nobs(ols))


ols <- lm(as.formula(paste(outcomes[i],"district_baraza+a21",sep="~")), data=baseline[(baseline$information == 0 & baseline$deliberation==0) | baseline$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID2[(baseline$information == 0 & baseline$deliberation==0) | baseline$district_baraza == 1 ], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_dist(i,outcomes, NULL, subset(baseline, ((information == 0 & deliberation==0) | district_baraza == 1)) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <-  RI_store$conf
res[2,5] <- RI_store$pval
}
df_ols_end[,4,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))

ols <- lm(as.formula(paste(outcomes[i],"district_baraza+a21",sep="~")), data=baseline[(baseline$information == 1 & baseline$deliberation==1) | baseline$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID2[(baseline$information == 1 & baseline$deliberation==1) | baseline$district_baraza == 1 ], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_dist(i,outcomes, NULL, subset(baseline, ((information == 1 & deliberation==1) | district_baraza == 1)) , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <-  RI_store$conf
res[2,5] <- RI_store$pval
}
df_ols_end[,5,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
}



### save results

 save(df_ols, file= paste(path,"report/results/df_ols_baseline.Rd", sep="/"))
 save(df_averages, file= paste(path,"report/results/df_averages_baseline.Rd", sep="/"))
 save(df_balance, file= paste(path,"report/results/df_balance_baseline.Rd", sep="/"))
 save(df_ols_end, file= paste(path,"report/results/df_ols_end.Rd", sep="/"))
 save(df_averages_end, file= paste(path,"report/results/df_averages_end.Rd", sep="/"))

if (final_verion_swith) { 
 save(df_ols, file= paste(path,"report/results/final/df_ols_baseline.Rd", sep="/"))
 save(df_averages, file= paste(path,"report/results/final/df_averages_baseline.Rd", sep="/"))
 save(df_balance, file= paste(path,"report/results/final/df_balance_baseline.Rd", sep="/"))
 save(df_ols_end, file= paste(path,"report/results/final/df_ols_end.Rd", sep="/"))
 save(df_averages_end, file= paste(path,"report/results/final/df_averages_end.Rd", sep="/"))
}

