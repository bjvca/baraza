rm(list=ls())
library(data.table)
library(doParallel)
library(ggplot2)
set.seed(12345)
#################################### functions ###############################
trim <- function(var, dataset, trim_perc=.05) {
### function for triming a variable in a dataset - replaces with NA
 dataset[var][dataset[var] < quantile(dataset[var],c(trim_perc/2,1-(trim_perc/2)), na.rm=T)[1] | dataset[var] > quantile(dataset[var], c(trim_perc/2,1-(trim_perc/2)),na.rm=T)[2] ] <- NA
return(dataset)
}
#############################################

cl <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
registerDoParallel(cl)


##hyper parameters
alpha <- 0.05 # Standard significance level 
sims <- 1000 # Number of simulations to conduct for each N 
#mde <- seq(from=.01, to=.1, by=.001) # The effect sizes we'll be considering 
mde <- seq(from=.05, to=.2, by=.0025)


baseline_outcomes <- c("base_unprotected","qc15","c10","a6","pub_health_access","maternal_health_access","d31","d43","d11","wait_time","d61","base_n_children","e5","e12", "e14","e22","e32","e45")
#baseline_outcomes <- c("b31","b44","base_inputs","b5144","b5146","c12source")
#baseline_outcomes <- c("c12source")
#baseline_outcomes <- c("b21")
for (outcome_index in 1:length(baseline_outcomes)) {
#outcome_index <- 1
print(baseline_outcomes[outcome_index])


#baseline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/baseline.csv")
#treats <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/final_list_5.csv")
#wget https://www.dropbox.com/s/sakp13112o1to6u/baseline.csv?dl=0
baseline <- read.csv("baseline.csv")
#wget https://www.dropbox.com/s/bsvv2ggud2g7yrj/final_list_5.csv?dl=0
treats <- read.csv("final_list_5.csv")

baseline$a23 <- as.character(baseline$a23)
baseline$a23[baseline$a23 == "LUWERO TC"] <- "LUWERO_TC"
baseline$a23[baseline$a23 == "SEMBABULE TC"] <- "SEMBABULE_TC"
baseline$a23[baseline$a23 == "RAKAI TC"] <- "RAKAI_TC"
baseline$a23[baseline$a23 == "NTUSI"] <- "NTUUSI"

baseline$b21 <-  as.numeric(baseline$b21=="Yes")
baseline$b31 <- as.numeric(baseline$b31=="Yes")
baseline$b44 <-  as.numeric(baseline$b44=="Yes")
baseline$b44[is.na(baseline$b44)] <- 0
baseline$base_inputs <- as.numeric(baseline$used_seed=="Yes" | baseline$used_fert=="Yes")
baseline$b5144 <- as.numeric(baseline$b5144=="Yes")
baseline$b5146 <- as.numeric(baseline$b5146=="Yes")
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

baseline$wait_time <- baseline$d410hh*60+ baseline$d410mm
baseline$wait_time <- log(baseline$wait_time + sqrt(baseline$wait_time ^ 2 + 1))
baseline <- trim("wait_time", baseline)

baseline$d61 <- as.numeric(baseline$d61=="Yes")

#children in public schools
baseline$base_n_children <- rowSums(cbind(baseline$e2bupe,baseline$e2aupe,baseline$e2fuse,baseline$e2muse), na.rm=T)
baseline$e5 <- rowMeans(cbind(as.numeric(baseline$e5upe) , as.numeric(baseline$e5use)), na.rm=T) 
baseline$e5[is.na(baseline$e5upe) & is.na(baseline$e5use)] <- NA
baseline$e5 <- log(baseline$e5 + sqrt(baseline$e5 ^ 2 + 1))
baseline <- trim("e5", baseline)



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


baseline_saved <- baseline
treats_saved <- treats

###power calcs for information treatment
### drop district barazas
treats <- subset(treats, district_baraza == 0 )
#baseline <- subset(baseline, district_baraza == 0 )
## should loop somewhere here


baseline$outcome <- unlist(baseline[baseline_outcomes[outcome_index]])
#baseline$outcome <- baseline$c12source

baseline_orig <- baseline[c("outcome","a21","a22","a23")]


power <- rep(NA, length(mde)) # Empty object to collect simulation estimates 
powers2 <- rep(NA, length(mde)) # Empty object to collect simulation estimates
powers3 <- rep(NA, length(mde)) # Empty object to collect simulation estimates
powers4 <- rep(NA, length(mde)) # Empty object to collect simulation estimates

#for sc baraza
#### Outer loop to vary MDE #### 
for (j in 1:length(mde)){ 
#print(j/length(mde)*100)
## for binary outcome
if (max(baseline_orig$outcome, na.rm=T) <= 1) {
 baseline_orig$Y1 <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=mde[j]) >= 1)
 baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=.45) >= 1)      
#### for cont outcome
} else {
baseline_orig$Y1 <- baseline_orig$outcome +  mde[j] #treatment potential outcome   
baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rnorm(length(baseline_orig$outcome ), mean = 0, sd = sd(baseline_orig$outcome, na.rm=T)))                                 
}
  significant.experiments <- rep(NA, sims) # Empty object to count significant experiments 

                              
  #### Inner loop to conduct experiments "sims" times over for each N #### 
  oper <- foreach (i = 1:sims,.combine=cbind) %dopar% {
 	treats$information <- sample(treats$information)
	treats$deliberation <- sample(treats$deliberation)
 	baseline_sim <- merge(baseline_orig,treats, by.x=c("a22","a23"), by.y=c("district","subcounty"))
       	baseline_sim$Y.sim <- baseline_sim$Y1*baseline_sim$information*baseline_sim$deliberation + baseline_sim$outcome*(1-baseline_sim$information*baseline_sim$deliberation) # Reveal 
        fit.sim <- lm(Y.sim ~ information:deliberation +base_out, data=baseline_sim[baseline_sim$information == baseline_sim$deliberation,]) # Do analysis (Simple regression) 
        p.value <- summary(fit.sim)$coefficients[3,4] # Extract p-values 
        return(p.value <= alpha) # Determine significance according to p <= 0.05
        }
  powers4[j] <- mean(oper) # store average success rate (power) for each N 
  } 

# for information
#### Outer loop to vary MDE #### 
for (j in 1:length(mde)){ 
#print(j/length(mde)*100)
## for binary outcome
if (max(baseline_orig$outcome, na.rm=T) <= 1) {
 baseline_orig$Y1 <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=mde[j]) >= 1)
 baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=.45) >= 1)      
#### for cont outcome
} else {
baseline_orig$Y1 <- baseline_orig$outcome +  mde[j] #treatment potential outcome   
baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rnorm(length(baseline_orig$outcome ), mean = 0, sd = sd(baseline_orig$outcome, na.rm=T)))                                 
}
  significant.experiments <- rep(NA, sims) # Empty object to count significant experiments 

                              
  #### Inner loop to conduct experiments "sims" times over for each N #### 
  oper <- foreach (i = 1:sims,.combine=cbind) %dopar% {
 	treats$information <- sample(treats$information)
 	baseline_sim <- merge(baseline_orig,treats, by.x=c("a22","a23"), by.y=c("district","subcounty"))
       	baseline_sim$Y.sim <- baseline_sim$Y1*baseline_sim$information + baseline_sim$outcome*(1-baseline_sim$information) # Reveal 
        fit.sim <- lm(Y.sim ~ information +base_out, data=baseline_sim) # Do analysis (Simple regression) 
        p.value <- summary(fit.sim)$coefficients[2,4] # Extract p-values 
        return(p.value <= alpha) # Determine significance according to p <= 0.05
        }
  power[j] <- mean(oper) # store average success rate (power) for each N 
  } 

###now for deliberation

#### Outer loop to vary MDE #### 
for (j in 1:length(mde)){ 
#print(j/length(mde)*100)
## for binary outcome
if (max(baseline_orig$outcome, na.rm=T) <= 1) {
 baseline_orig$Y1 <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=mde[j]) >= 1)
 baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=.45) >= 1)      
#### for cont outcome
} else {
baseline_orig$Y1 <- baseline_orig$outcome +  mde[j] #treatment potential outcome   
baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rnorm(length(baseline_orig$outcome ), mean = 0, sd = sd(baseline_orig$outcome, na.rm=T)))                                 
}          
  significant.experiments <- rep(NA, sims) # Empty object to count significant experiments 

                              
  #### Inner loop to conduct experiments "sims" times over for each N #### 
  oper <- foreach (i = 1:sims,.combine=cbind) %dopar% {

 		treats$deliberation <- sample(treats$deliberation)
 	baseline_sim <- merge(baseline_orig,treats, by.x=c("a22","a23"), by.y=c("district","subcounty"))
       	baseline_sim$Y.sim <- baseline_sim$Y1*baseline_sim$deliberation + baseline_sim$outcome*(1-baseline_sim$deliberation) # Reveal 
        fit.sim <- lm(Y.sim ~ deliberation +base_out, data=baseline_sim) # Do analysis (Simple regression) 
        p.value <- summary(fit.sim)$coefficients[2,4] # Extract p-values 
        return(p.value <= alpha)
        }
  powers2[j] <- mean(oper) # store average success rate (power) for each N 
  } 

##now compare district to sc baraza

treats <- subset(treats_saved, (information == 1 & deliberation ==1 ) | district_baraza == 1 )
baseline <- merge(baseline_saved, treats,  by.x=c("a22","a23"), by.y=c("district","subcounty"))
treats <- aggregate(treats$district_baraza,list(treats$district), max)
names(treats) <- c("district","district_baraza")
## should loop somewhere here
baseline$outcome <-  unlist(baseline[baseline_outcomes[outcome_index]])
#baseline$outcome <- baseline$c12source

baseline_orig <- baseline[c("outcome","a21","a22","a23")]



#### Outer loop to vary MDE #### 
for (j in 1:length(mde)){ 
## for binary outcome
if (max(baseline_orig$outcome, na.rm=T) <= 1) {
 baseline_orig$Y1 <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=mde[j]) >= 1)
 baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=.45) >= 1)      
#### for cont outcome
} else {
baseline_orig$Y1 <- baseline_orig$outcome +  mde[j] #treatment potential outcome   
baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rnorm(length(baseline_orig$outcome ), mean = 0, sd = sd(baseline_orig$outcome, na.rm=T)))                                 
}
significant.experiments <- rep(NA, sims) # Empty object to count significant experiments 
                              
  #### Inner loop to conduct experiments "sims" times over for each N #### 
  oper <- foreach (i = 1:sims,.combine=cbind) %dopar% {
		
		### we need to permute the treatment assignment at the cluster level - how to do that???
 		treats$district_baraza <- sample(treats$district_baraza)
 	baseline_sim <- merge(baseline_orig,treats, by.x=c("a22"), by.y=c("district"))
       	baseline_sim$Y.sim <- baseline_sim$Y1*baseline_sim$district_baraza + baseline_sim$outcome*(1-baseline_sim$district_baraza) # Reveal outcomes according to assignment 
        #fit.sim <- lm(Y.sim ~ information, data=baseline_sim) # Do analysis (Simple regression) 
        fit.sim <- lm(Y.sim ~ district_baraza +base_out, data=baseline_sim) # Do analysis (Simple regression) 


#test <- data.frame(cbind(baseline_sim$base_out,baseline_sim$information, baseline_sim$deliberation,  baseline_sim$a21))
#names(test) <- c("Y.sim","information","deliberation","a21")
#test$time <- 0
#test2 <- data.frame(cbind(baseline_sim$Y.sim,baseline_sim$information, baseline_sim$deliberation,  baseline_sim$a21))
#names(test2) <- c("Y.sim","information","deliberation","a21")
#test2$time <- 1
#baseline_sim <- rbind(test,test2)

#  fit.sim3 <- lm(Y.sim ~ information*time, data=baseline_sim) # Do analysis (Simple regression) 
        p.value <- summary(fit.sim)$coefficients[2,4] # Extract p-values 
        return(p.value <= alpha) # Determine significance according to p <= 0.05
#        p.value2 <- summary(fit.sim2)$coefficients[2,4] # Extract p-values 
#        significant.experiments2[i] <- (p.value2 <= alpha) # Determine significance according to p <= 0.05
#        p.value3 <- summary(fit.sim3)$coefficients[4,4] # Extract p-values 
#        significant.experiments3[i] <- (p.value3 <= alpha) # Determine significance according to p <= 0.05

        }
#  powers[j] <- mean(significant.experiments) # store average success rate (power) for each N 
#  powers2[j] <- mean(significant.experiments2) # store average success rate (power) for each N 
  powers3[j] <- mean(oper) # store average success rate (power) for each N 
  } 



#png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_extension.png", units="px", height=3200, width= 3200, res=600)
#png("power_extension.png", units="px", height=3200, width= 3200, res=600)
# plot(mde, power, type="l", ylim=c(0,1), lwd=2.5)
#lines(mde,powers2,col="green")
#lines(mde,powers3,col="blue")
#abline(h=.8, col="red")
#dev.off()

df1 <- data.frame(mde,power)
df1$hypo <- "information"
df2 <-  data.frame(mde,powers2)
names(df2) <- c("mde","power")
df2$hypo <- "deliberation"
df3 <-  data.frame(mde,powers3)
names(df3)<- c("mde","power")
df3$hypo <- "level"
df4 <-  data.frame(mde,powers4)
names(df4)<- c("mde","power")
df4$hypo <- "sc_baraza"

df <- rbind(df1,df2,df3,df4)
write.csv(df, file = paste(baseline_outcomes[outcome_index],"csv",sep="."))
}

### could not figure out how to easily automate this so do this manually:

##to copy entire folder from AWS to local machine, use:
scp -i "bjornkey.pem" ubuntu@ec2-3-249-40-192.eu-west-1.compute.amazonaws.com:/home/ubuntu/* "/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/"
## in a local terminal

###f***ing R: want to save this as png but for some reason I can not use dynamic names with he png command
library(ggplot2)
df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/b21.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_extension.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/b31.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_field_visits.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/b44.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_NAADS_in_village.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/base_inputs.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_inputs.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/b5144.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_advice_market_comm.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/b5146.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_advice_market_coop.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

### infra

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/base_unprotected.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_unprotected.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()


df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/c12source.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_dist_water.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/qc15.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_wait_water.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()


df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/c10.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_water_user_com.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()


df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/a6.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_dist_road.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/d31.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_vht.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()


df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/d43.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_dist_gov.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()


df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/pub_health_access.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_pub_health_access.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()


df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/maternal_health_access.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_maternal_health_access.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/d11.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_school_work_missed.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()


df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/wait_time_health.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_wait_time.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

base_n_children.csv 

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/base_n_children.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_base_n_children.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

base_n_children
e12.csv  

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/e12.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_complete_fence.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()


df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/e5.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_dist_school.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/e14.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_water_school.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/e22.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_SMC_school.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/e32.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_info_SMC_school.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

df <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/e45.csv")
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_inspection_school.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B","#000000"))  + scale_linetype_manual(values=c("solid","dashed", "twodash","dotted"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

   


