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

#baseline$c12source <- log(baseline$c12source + sqrt(baseline$c12source ^ 2 + 1))
#baseline <- trim("c12source", baseline)

###power calcs for information treatment
### drop district barazas
treats <- subset(treats, district_baraza == 0 )
baseline <- subset(baseline, district_baraza == 0 )
## should loop somewhere here
baseline$outcome <- baseline$b21
#baseline$outcome <- baseline$c12source

baseline_orig <- baseline[c("outcome","a21","a22","a23")]



mde <- seq(from=.01, to=.06, by=.001) # The effect sizes we'll be considering 
power <- rep(NA, length(mde)) # Empty object to collect simulation estimates 
powers2 <- rep(NA, length(mde)) # Empty object to collect simulation estimates
powers3 <- rep(NA, length(mde)) # Empty object to collect simulation estimates
alpha <- 0.05 # Standard significance level 
sims <- 2000 # Number of simulations to conduct for each N 

#### Outer loop to vary the number of subjects #### 
for (j in 1:length(mde)){ 
print(j/length(mde)*100)
## for binary outcome
 baseline_orig$Y1 <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=mde[j]) >= 1)
 baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=.45) >= 1)      
#### for cont outcome
#baseline_orig$Y1 <- baseline_orig$outcome +  mde[j] #treatment potential outcome   
#baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rnorm(length(baseline_orig$outcome ), mean = 0, sd = sd(baseline_orig$outcome, na.rm=T)))                                 
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

#### Outer loop to vary the number of subjects #### 
for (j in 1:length(mde)){ 
print(j/length(mde)*100)
 ## for binary outcme
  baseline_orig$Y1 <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=mde[j]) >= 1)
  baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=.45) >= 1)      
### for cont outcome
#baseline_orig$Y1 <- baseline_orig$outcome +  mde[j] #treatment potential outcome   
#baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rnorm(length(baseline_orig$outcome ), mean = 0, sd = sd(baseline_orig$outcome, na.rm=T)))           
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




#baseline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/baseline.csv")
#treats <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/final_list_5.csv")
baseline <- read.csv("baseline.csv")
treats <- read.csv("final_list_5.csv")

baseline$a23 <- as.character(baseline$a23)
baseline$a23[baseline$a23 == "LUWERO TC"] <- "LUWERO_TC"
baseline$a23[baseline$a23 == "SEMBABULE TC"] <- "SEMBABULE_TC"
baseline$a23[baseline$a23 == "RAKAI TC"] <- "RAKAI_TC"
baseline$a23[baseline$a23 == "NTUSI"] <- "NTUUSI"

baseline$b21 <-  as.numeric(baseline$b21=="Yes")

#baseline$c12source <- log(baseline$c12source + sqrt(baseline$c12source ^ 2 + 1))
#baseline <- trim("c12source", baseline)


###power calcs for information treatment
### drop district barazas
treats <- subset(treats, (information == 1 & deliberation ==1 ) | district_baraza == 1 )
baseline <- merge(baseline, treats,  by.x=c("a22","a23"), by.y=c("district","subcounty"))
treats <- aggregate(treats$district_baraza,list(treats$district), max)
names(treats) <- c("district","district_baraza")
## should loop somewhere here
baseline$outcome <- baseline$b21
#baseline$outcome <- baseline$c12source

baseline_orig <- baseline[c("outcome","a21","a22","a23")]



#### Outer loop to vary the number of subjects #### 
for (j in 1:length(mde)){ 
print(j/length(mde)*100)  
  ## for binary outcme
  baseline_orig$Y1 <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=mde[j]) >= 1)
  baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rbinom(n=length(baseline_orig$outcome), size=1, prob=.45) >= 1)      
#### for cont outcome
#baseline_orig$Y1 <- baseline_orig$outcome +  mde[j] #treatment potential outcome   
#baseline_orig$base_out <- as.numeric(baseline_orig$outcome +  rnorm(length(baseline_orig$outcome ), mean = 0, sd = sd(baseline_orig$outcome, na.rm=T)))  
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

df <- rbind(df1,df2,df3)
png("power_extension.png", units="px", height=3200, width= 3200, res=600)
ggplot(df, aes(x = mde, y = power, group = hypo)) +  geom_line(aes(color=hypo, linetype=hypo), size=1)  + scale_color_manual(values=c("#CCCCCC", "#6E8DAB", "#104E8B"))  + scale_linetype_manual(values=c("solid","dashed", "twodash"))+ geom_hline(yintercept = .8, colour =  "red", size=1)
dev.off()

##to copy from AWS to local machine, use:
scp -i "bjornkey.pem" ubuntu@ec2-52-31-191-248.eu-west-1.compute.amazonaws.com:/home/ubuntu/power_extension.png "/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/power_extension.png"


