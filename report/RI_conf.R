rm(list=ls())
library(data.table)
dta <- read.csv("test.csv")
library(clubSandwich)
library(doParallel)
library(data.table)

outcomes <- c("baraza.B2","baraza.B3","baraza.B4.1","inputs","baraza.B5.2","baraza.B5.3","ag_index","unprotected", "baraza.C1.2", "baraza.C1.3","baraza.C2.3","baraza.A6","infra_index","baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1.2",  "baraza.D4.6","baraza.D6","health_index","n_children","baraza.E5","baraza.E12","baraza.E14","baraza.E22","baraza.E32","baraza.E45","education_index", "pub_service_index")
baseline_outcomes <- c("b21","b31","b44","base_inputs","b5144","b5146","base_ag_index","base_unprotected","c12source", "qc15","c10","a6","base_infra_index","pub_health_access","maternal_health_access","d31","d43","tot_sick","wait_time","d61","base_health_index","base_n_children","e5","e12", "e14","e22","e32","e45","base_education_index","base_pub_service_index")


i <- 1

ols <- lm(as.formula(paste(paste(outcomes[i],"information:deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta[dta$district_baraza == 0 & (dta$information == dta$deliberation),]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID[dta$district_baraza == 0 & (dta$information == dta$deliberation)], type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)


dta <- subset(dta, district_baraza == 0 & (information == deliberation))

dta_sim <- dta
dta_sim$dep <- as.numeric(unlist(dta_sim[outcomes[i]]))

	
 		

cl <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
registerDoParallel(cl)

treat_nr <- sum(data.frame(aggregate(dta$information, list(dta$subcounty),mean))[,2])
scs <- levels(factor(dta$subcounty))
dta_sim$pot_out_0 <- ifelse(dta_sim$information ==0, dta_sim$dep, dta_sim$dep - coef(ols)[6])
dta_sim$pot_out_1 <- ifelse(dta_sim$information ==1, dta_sim$dep, dta_sim$dep + coef(ols)[6])


nr_repl <- 10000

oper <- foreach (repl = 1:nr_repl,.combine=cbind,.packages = c("data.table")) %dopar% {
dta_sim$information <- ifelse(dta_sim$subcounty %in% sample(scs,treat_nr), 1, 0)
dta_sim$dep <- ifelse(dta_sim$information ==1, dta_sim$pot_out_1, dta_sim$pot_out_0)

return(coef(lm(as.formula(paste("dep~information+a21",baseline_outcomes[i],sep="+")), data=dta_sim))[2])

}

quantile(oper,c(.025,.975))


