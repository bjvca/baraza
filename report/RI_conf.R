rm(list=ls())
library(data.table)
dta <- read.csv("test.csv")
library(clubSandwich)
library(doParallel)
library(data.table)
library(dplyr)

dta <- read.csv("/home/bjvca/test.csv")

outcomes <- c("baraza.B2","baraza.B3","baraza.B4.1","inputs","baraza.B5.2","baraza.B5.3","ag_index","unprotected", "baraza.C1.2", "baraza.C1.3","baraza.C2.3","baraza.A6","infra_index","baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1.2",  "baraza.D4.6","baraza.D6","health_index","n_children","baraza.E5","baraza.E12","baraza.E14","baraza.E22","baraza.E32","baraza.E45","education_index", "pub_service_index")
baseline_outcomes <- c("b21","b31","b44","base_inputs","b5144","b5146","base_ag_index","base_unprotected","c12source", "qc15","c10","a6","base_infra_index","pub_health_access","maternal_health_access","d31","d43","tot_sick","wait_time","d61","base_health_index","base_n_children","e5","e12", "e14","e22","e32","e45","base_education_index","base_pub_service_index")


i <- 1


ols <- lm(as.formula(paste(paste(outcomes[i],"information*deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta[dta$district_baraza == 0 ,]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID[dta$district_baraza == 0], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

ols <- lm(as.formula(paste(paste(outcomes[i],"information:deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta[dta$district_baraza == 0 & (dta$information == dta$deliberation),]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID[dta$district_baraza == 0 & (dta$information == dta$deliberation)], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

		

cl <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
registerDoParallel(cl)

## function definitions 
RI_conf_sc <- function(formula,var_interest,  dta_sim , nr_repl = 10000, sig = c(.025,.975)) {

formula1 <- as.formula(paste(paste(outcomes[i],"information:deliberation+a21",sep="~"),baseline_outcomes[i],sep="+"))
formula2 <- as.formula(paste(paste(outcomes[i],"information*deliberation+a21",sep="~"),baseline_outcomes[i],sep="+"))
dta_sim <- subset(dta, district_baraza == 0)
var_interest <- "information:deliberation"
nr_repl <- 1000
sig <- c(.025,.975)

dta_sim <- dta_sim %>%  mutate(clusterID = group_indices(., district, subcounty))

ols_1 <- lm(formula1, data=dta_sim[dta_sim$deliberation == dta_sim$information,]) 
ols_2 <- lm(formula2, data=dta_sim) 
dta_sim$dep <- as.numeric(unlist(dta_sim[as.character(formula1[[2]])]))

treat_nrs <- table(data.frame(aggregate(dta_sim[c("information","deliberation")], list(dta_sim$clusterID),mean))[,2:3])

#### this is hard...


dta_sim$pot_out_0_1 <- NA
dta_sim$pot_out_0_1[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0]
dta_sim$pot_out_0_1[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] - coef(ols_1)["information:deliberation"]

dta_sim$pot_out_1_1 <- NA
dta_sim$pot_out_1_1[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1]
dta_sim$pot_out_1_1[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] + coef(ols_1)["information:deliberation"]


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

oper <- foreach (repl = 1:nr_repl,.combine=cbind) %dopar% {
	#do the permuations	
print(repl)	
	perm_treat <- data.frame(cbind(sample(c(rep("C",treat_nrs[1,1]), rep("D",treat_nrs[1,2]), rep("I",treat_nrs[2,1]),rep("B",treat_nrs[2,2]))),names(table(dta_sim$clusterID))))  
	names(perm_treat) <- c("perm_treat","clusterID")
	dta_perm <- merge(dta_sim, perm_treat, by.x="clusterID", by.y="clusterID")

	dta_perm$information <- ifelse(dta_perm$perm_treat %in% c("I","B"), 1, 0)
	dta_perm$deliberation <- ifelse(dta_perm$perm_treat %in% c("D","B"), 1, 0)



	dta_perm$dep_2 <- NA
dta_perm$dep_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1] <- dta_perm$pot_out_11_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1]
dta_perm$dep_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0] <- dta_perm$pot_out_00_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0] 
dta_perm$dep_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==0] <- dta_perm$pot_out_10_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==0] 
dta_perm$dep_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==1] <- dta_perm$pot_out_01_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==1] 

dta_perm[outcomes[i]] <- dta_perm$dep_2
return(coef(lm(formula2, data=dta_perm))["information:deliberation"])

}
return(quantile(oper,sig))
}

RI_conf_sc(as.formula(paste(paste(outcomes[i],"information*deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")),"information", dta_sim= subset(dta, district_baraza == 0))


