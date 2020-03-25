rm(list=ls())
library(dplyr)
library(ggplot2)
library(MatchIt) 
library(multiwayvcov)
library(plm)
library(lmtest)
library(clubSandwich)
library(moments)
set.seed(12345) #not needed for final version?

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

################################################################## end of funtions declarations
if (Sys.info()['sysname'] =="Linux") {
path <- "/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline"
} else {
path <- "/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline"
}

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
baseline$base_unprotected <- as.numeric(( baseline$c11a %in%  c("Surface water","Bottled water","Cart with small tank","Unprotected dug well","Unprotected spring","Tanker truck"))    )
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
##dataset to be used for ancova


###init arrays to store results
df_ols <- array(NA,dim=c(6,4,length(outcomes)))
df_averages <- array(NA,dim=c(2,length(outcomes)))


for (i in 1:length(outcomes)) {
#print(i)

df_averages[1,i] <- mean(as.matrix(baseline[outcomes[i]]), na.rm=T)
df_averages[2,i] <- sd(as.matrix(baseline[outcomes[i]]), na.rm=T)

### simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+a21",sep="~")), data=baseline[baseline$district_baraza == 0,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID[baseline$district_baraza == 0], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5], nobs(ols))
df_ols[,1,i] <- c(res[7,1],res[7,2],res[7,5], conf[7,4],conf[7,5], nobs(ols))

ols <- lm(as.formula(paste(outcomes[i],"district_baraza+a21",sep="~")), data=baseline[(baseline$information == 1 & baseline$deliberation==1) | baseline$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID2[(baseline$information == 1 & baseline$deliberation==1) | baseline$district_baraza == 1 ], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
df_ols[,4,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
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

for (i in 1:length(outcomes)) {
## simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information.x*deliberation.x+a21",sep="~")), data=baseline_treat_info) 
vcov_cluster <- vcovCR(ols, cluster = baseline_treat_info$clusterID, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_balance[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))

}

#for deliberation
baseline_treat <- merge(baseline, treats, all.x=T, by.y=c("district","subcounty"), by.x=c("a22","a23"))

baseline_treat <- subset(baseline_treat,district_baraza.x == 0)
baseline_treat$deliberation.y[is.na(baseline_treat$deliberation.y)] <- 0
baseline_treat_delib <- subset(baseline_treat, deliberation.y!=1)

for (i in 1:length(outcomes)) {
### simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information.x*deliberation.x+a21",sep="~")), data=baseline_treat_delib) 
vcov_cluster <- vcovCR(ols, cluster = baseline_treat_delib$clusterID, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_balance[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5], nobs(ols))
}

#for interaction
baseline_treat <- merge(baseline, treats, all.x=T, by.y=c("district","subcounty"), by.x=c("a22","a23"))

baseline_treat <- subset(baseline_treat,district_baraza.x == 0)
### drop if information.y==1 & information.y==1
baseline_treat$information.y[is.na(baseline_treat$information.y)] <- 0
baseline_treat_info <- subset(baseline_treat, information.y!=1 | deliberation.y!=1)

for (i in 1:length(outcomes)) {
## simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information.x*deliberation.x+a21",sep="~")), data=baseline_treat_info) 
vcov_cluster <- vcovCR(ols, cluster = baseline_treat_info$clusterID, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_balance[,1,i] <- c(res[7,1],res[7,2],res[7,5], conf[7,4],conf[7,5], nobs(ols))

}


### redo balance table, but now between actual treated and matched controls
baseline$information <- NULL
baseline$deliberation <- NULL
baseline$district_baraza <- NULL


baseline <- merge(baseline, treats, all.y=T, by.y=c("district","subcounty"), by.x=c("a22","a23"))

###init arrays to store results
df_ols_end <- array(NA,dim=c(6,4,length(outcomes)))
df_averages_end <- array(NA,dim=c(2,length(outcomes)))


for (i in 1:length(outcomes)) {
#print(i)

df_averages_end[1,i] <- mean(as.matrix(baseline[outcomes[i]]), na.rm=T)
df_averages_end[2,i] <- sd(as.matrix(baseline[outcomes[i]]), na.rm=T)

### simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+a21",sep="~")), data=baseline[baseline$district_baraza == 0,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID[baseline$district_baraza == 0], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_ols_end[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
df_ols_end[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5], nobs(ols))
df_ols_end[,1,i] <- c(res[7,1],res[7,2],res[7,5], conf[7,4],conf[7,5], nobs(ols))

ols <- lm(as.formula(paste(outcomes[i],"district_baraza+a21",sep="~")), data=baseline[(baseline$information == 1 & baseline$deliberation==1) | baseline$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = baseline$clusterID2[(baseline$information == 1 & baseline$deliberation==1) | baseline$district_baraza == 1 ], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
df_ols_end[,4,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5], nobs(ols))
}




