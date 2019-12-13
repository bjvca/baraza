### to solve - this also gives me the same results regardless of random seed...


rm(list=ls())
library(dplyr)
library(ggplot2)
library(MatchIt) 
library(multiwayvcov)
library(plm)
library(lmtest)
library(clubSandwich)
library(moments)
set.seed(66666) #not needed for final version?

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

### wrapper function to make summary forest plots
credplot.gg <- function(d,units, hypo, axlabs, lim){
 # d is a data frame with 4 columns
 # d$x gives variable names
 # d$y gives center point
 # d$ylo gives lower limits
 # d$yhi gives upper limits
 require(ggplot2)
 p <- ggplot(d, aes(x=x, y=y, ymin=ylo, ymax=yhi, colour=as.factor(grp)))+
 geom_pointrange(position=position_dodge(-.4), size=1)+
 geom_hline(yintercept = 0, linetype=2)+
 coord_flip(ylim = c(-lim,lim))+
 xlab('') + ylab(units)+ labs(title=hypo)  + theme_minimal()+ theme(axis.text=element_text(size=18),
        axis.title=element_text(size=14,face="bold"),legend.text=element_text(size=18), plot.title = element_text(size=22,hjust = 0.5), legend.title=element_blank())+
    geom_errorbar(aes(ymin=ylo, ymax=yhi),position=position_dodge(-.4),width=0,cex=2.5) + scale_colour_manual(values = c("#CCCCCC", "#6E8DAB", "#104E8B")) + scale_x_discrete(labels=axlabs)
 return(p)
}

################################################################## end of funtions declarations

# takes raw data (baseline and endline), makes it anonymous and pust in into the data/public folder, ready to be analysed by the code chucks below
#source("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/anonyize.R")


### this is for the dummy endline for the mock report - I read in a dummy endline of 3 households just to get the correct variable names
endline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/endline.csv")[10:403]
### I then merge with the sampling list to basically create an empty endline dataset
### and merge in the treatments
list <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/sampling_list_hh.csv")[c("hhid","a21","a22","a23")]
endline <- merge(list, endline, by="hhid", all.x=T)
treats <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/final_list.csv")
endline <- merge(treats, endline, by.x=c("district","subcounty"), by.y=c("a22","a23"))

## baseline not needed in this first section, but used to generate fake data
baseline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/baseline.csv")
baseline$a23 <- as.character(baseline$a23)
baseline$a23[baseline$a23 == "LUWERO TC"] <- "LUWERO_TC"
baseline$a23[baseline$a23 == "SEMBABULE TC"] <- "SEMBABULE_TC"

baseline$b21 <-  as.numeric(baseline$b21=="Yes")
baseline$b31 <-  as.numeric(baseline$b31=="Yes")
baseline$b44 <-  as.numeric(baseline$b44=="Yes")
baseline$b44[is.na(baseline$b44)] <- 0
baseline$base_inputs <- as.numeric(baseline$used_seed=="Yes" | baseline$used_fert=="Yes")
baseline$b5144 <- as.numeric(baseline$b5144=="Yes")
baseline$b5146 <- as.numeric(baseline$b5146=="Yes")
##use of unprotected water sources in dry season
baseline$base_unprotected <- as.numeric(( baseline$c11a %in%  c("Rain water","Surface water","Tube well or borehole","Unprotected dug well","Unprotected spring"))    )
### is there are water committee
baseline$c10 <- as.numeric(baseline$c10=="Yes")
baseline$c12source <- log(baseline$c12source + sqrt(baseline$c12source ^ 2 + 1))
baseline <- trim("c12source", baseline)
baseline$qc15 <- log(baseline$qc15 + sqrt(baseline$qc15 ^ 2 + 1))
baseline <- trim("qc15", baseline)
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

#

##need to take logs?
baseline$log_farmsize <- log(baseline$farmsize + sqrt(baseline$farmsize ^ 2 + 1))
baseline$log_farmsize[is.infinite(baseline$log_farmsize)] <- NA
baseline <- trim("farmsize", baseline)

baseline$ironroof <- as.numeric(baseline$a512 =="Corrugated iron sheets")
baseline$improved_wall <- as.numeric(baseline$a513 %in% c("Mud_bricks_burnt_bricks","Concrete_blocks") )
baseline$head_sec <- as.numeric(baseline$a36) > 15

baseline_matching <- merge(baseline,treats, by.x=c("a22","a23"), by.y=c("district","subcounty"))

#define endline variables
endline$baraza.B3 <- 0
endline$baraza.B3 <- endline$baraza.B3 ==1 |  endline$baraza.B3.3 ==1
endline$inputs <- 0
endline$inputs <- as.numeric(endline$baraza.B1==1 | endline$baraza.B1.5==1) 
endline$unprotected <- (as.numeric(endline$baraza.C1 %in% c(4,6,8,9,12)) )

### here we simulate endline variables - remove if endline data is in
endline$baraza.B2  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b21 == 1, na.rm=T))
###visits

endline$baraza.B3 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b31 == 1, na.rm=T))
###naads in village
endline$baraza.B4.1  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b44 == 1, na.rm=T))
###simulate an effect on this one
endline$inputs <- rbinom(n=length(endline$inputs),size=1,prob=mean(baseline$base_inputs, na.rm=T)) 
endline$baraza.B5.2  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b5144 == 1, na.rm=T))
endline$baraza.B5.3  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b5146 == 1, na.rm=T))

endline$unprotected <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$base_unprotected ==1, na.rm=T))
endline$baraza.C1.2 <- sample(baseline$c12source,dim(endline)[1])  ### this needs to be inverse hypersine transformed and trimmed in final version
#endline$baraza.C1.2 <-  log(endline$baraza.C1.2 + sqrt(endline$baraza.C1.2 ^ 2 + 1))
#endline <- trim("baraza.C1.2",endline)
endline$baraza.C1.3 <- sample(baseline$qc15,dim(endline)[1]) ### this needs to be inverse hypersine transformed and trimmed in final version
#endline$baraza.C1.3 <- log(endline$baraza.C1.3 + sqrt(endline$baraza.C1.3 ^ 2 + 1))
#endline <- trim("baraza.C1.3",endline)
endline$baraza.C2.3 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$c10 ==1, na.rm=T))
### access to public health if sick
endline$baraza.D2 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$pub_health_access, na.rm=T))
endline$baraza.D2.4 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$maternal_health_access , na.rm=T))
endline$baraza.D3 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$pub_health_access, na.rm=T))
endline$baraza.D4.2 <- sample(baseline$d43 ,dim(endline)[1])  ### this needs to be inverse hypersine transformed and trimmed in final version
#endline$baraza.D4.2 <- log(endline$baraza.D4.2 + sqrt(endline$baraza.D4.2 ^ 2 + 1))
#endline <- trim("baraza.D4.2",endline)
#health outcome - less people sick 
endline$baraza.D1  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$d11, na.rm=T))
endline$baraza.D1.2 <- sample(baseline$tot_sick ,dim(endline)[1]) 
endline$baraza.D1.2[is.na(endline$baraza.D1.2)] <- 0


##ag

#1 access to extension, 
#2 visit extension officer or demo site, 
#2 farmer in naads supported group, 
#4 uses inputs (seed or fert), 
#5 support in marketing from Village procurement committe/Village farmers forum/Village farmers forum executive; 
#6 support in marketing from Cooperative

endline <- endline[!duplicated(endline$hhid),]
# 7 #make an ag index
endline <- FW_index(c("baraza.B2","baraza.B3","baraza.B4.1","inputs","baraza.B5.2","baraza.B5.3"),data=endline)
names(endline)[names(endline) == 'index'] <- 'ag_index'
baseline <- FW_index(c("b21","b31","b44","base_inputs","b5144","b5146"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_ag_index'
baseline_matching <- FW_index(c("b21","b31","b44","base_inputs","b5144","b5146"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_ag_index'

##make an infrastructure index - note the revcols argument as the first 3 outcomes are "more is worse"
#8 unprotected water source in dry season, 
#9 distance to water source in dry season, 
#10 waiting time, 
#11 is there a water commitee?
#12  #index
endline <- FW_index(c("unprotected", "baraza.C1.2", "baraza.C1.3","baraza.C2.3"),revcols=c(1,2,3),data=endline)
names(endline)[names(endline) == 'index'] <- 'infra_index'
baseline <- FW_index(c("base_unprotected","c12source", "qc15","c10"),revcols=c(1,2,3),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_infra_index'
baseline_matching <- FW_index(c("base_unprotected","c12source", "qc15","c10"),revcols=c(1,2,3),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_infra_index'

##make a health index
#13 pub health access
#14 maternal health acess
#15 is there a VHT?
#16 distance to gvt health facility
##17 has been sick
###18 number of days sick
## 19 index
endline <- FW_index(c("baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1",  "baraza.D1.2"),revcols=c(4,5,6),data=endline)
names(endline)[names(endline) == 'index'] <- 'health_index'
baseline <- FW_index(c("pub_health_access","maternal_health_access","d31","d43","d11","tot_sick"),revcols=c(4,5,6),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_health_index'
baseline_matching <- FW_index(c("pub_health_access","maternal_health_access","d31","d43","d11","tot_sick"),revcols=c(4,5,6),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_health_index'

#20 make an index of indices
endline <- FW_index(c("ag_index","infra_index","health_index"),data=endline)
names(endline)[names(endline) == 'index'] <- 'pub_service_index'
baseline <- FW_index(c("base_ag_index","base_infra_index","base_health_index"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_pub_service_index'
baseline_matching <- FW_index(c("base_ag_index","base_infra_index","base_health_index"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_pub_service_index'



outcomes <- c("baraza.B2","baraza.B3","baraza.B4.1","inputs","baraza.B5.2","baraza.B5.3","ag_index","unprotected", "baraza.C1.2", "baraza.C1.3","baraza.C2.3","infra_index","baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1",  "baraza.D1.2","health_index","pub_service_index")
baseline_outcomes <- c("b21","b31","b44","base_inputs","b5144","b5146","base_ag_index","base_unprotected","c12source", "qc15","c10","base_infra_index","pub_health_access","maternal_health_access","d31","d43","d11","tot_sick","base_health_index","base_pub_service_index")


#create unique ID for clustering based on district and subcounty
endline <- endline %>%  mutate(clusterID = group_indices(., district, subcounty))
##dataset to be used for ancova
dta <- merge(endline, baseline[, -which(names(baseline)=="a21")], by="hhid")
endline$time <- 1
baseline$time <- 0
baseline$information <- 0
baseline$deliberation <- 0
baseline$information[baseline$treat=="info" | baseline$treat=="scbza"] <- 1 
baseline$deliberation[baseline$treat=="delib" | baseline$treat=="scbza"] <- 1 
#I am droping the district baraza's here... I completely forgot about those...
baseline <- subset(baseline, treat != "dbza")
### merge in clusterID for standard error clustering in dif-in-dif
baseline <- merge(baseline, endline[c("hhid","clusterID")], by="hhid", all.y=T)

baseline <- baseline[c("information","deliberation","time",baseline_outcomes )]
names(baseline) <- c("information","deliberation","time","clusterID",outcomes )


dta_long <- rbind(endline[c("information","deliberation","time", "clusterID",outcomes)], baseline[ c("information","deliberation","time",outcomes )])


###init arrays to store results
df_ols <- array(NA,dim=c(6,3,length(outcomes)))
df_ancova <- array(NA,dim=c(6,3,length(outcomes)))
df_dif_in_dif <- array(NA,dim=c(6,3,length(outcomes)))
df_matcher <- array(NA,dim=c(6,3,length(outcomes)))
df_averages <- array(NA,dim=c(2,length(outcomes)))

for (i in 1:length(outcomes)) {
print(i)

df_averages[1,i] <- mean(as.matrix(endline[outcomes[i]]), na.rm=T)
df_averages[2,i] <- sd(as.matrix(endline[outcomes[i]]), na.rm=T)

### simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+a21",sep="~")), data=endline) 
vcov_cluster <- vcovCR(ols, cluster = endline$clusterID, type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_ols[,1,i] <- c(res[2,1],res[2,2],res[2,4], conf[2,4],conf[2,5], nobs(ols))
df_ols[,2,i] <- c(res[3,1],res[3,2],res[3,4], conf[3,4],conf[3,5], nobs(ols))
df_ols[,3,i] <- c(res[7,1],res[7,2],res[7,4], conf[7,4],conf[7,5], nobs(ols))

##ancova
## merge in baseline

ols <- lm(as.formula(paste(paste(outcomes[i],"information*deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID, type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_ancova[,1,i] <- c(res[2,1],res[2,2],res[2,4], conf[2,4],conf[2,5], nobs(ols))
df_ancova[,2,i] <- c(res[3,1],res[3,2],res[3,4], conf[3,4],conf[3,5], nobs(ols))
df_ancova[,3,i] <- c(res[8,1],res[8,2],res[8,4], conf[8,4],conf[8,5], nobs(ols))


## dif-in-dif
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=dta_long)
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID, type = "CR0")
res <- coef(summary(lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=dta_long) ))
conf <- confint(lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=dta_long))

df_dif_in_dif[,1,i] <- c(res[6,1],res[6,2],res[6,4], conf[6,1], conf[6,2], nobs(ols))
df_dif_in_dif[,2,i] <- c(res[7,1],res[7,2],res[7,4], conf[7,1], conf[7,2], nobs(ols))
df_dif_in_dif[,3,i] <- c(res[8,1],res[8,2],res[8,4], conf[8,1], conf[8,2], nobs(ols))


###matched dif-in-dif

baseline_complete <- baseline_matching[complete.cases(baseline_matching[c("information", "deliberation","hhsize","femhead","agehead","log_farmsize","ironroof","improved_wall","has_phone","head_sec","a26a","a26b","hhid",baseline_outcomes[i])]),] 
baseline_complete <- baseline_complete[c("information","deliberation","hhsize","femhead","agehead","log_farmsize","ironroof","improved_wall","has_phone","head_sec","a26a","a26b","hhid",baseline_outcomes[i])]

####matching for information
nearest.match <- matchit(formula = information ~ hhsize + femhead + agehead + log_farmsize + ironroof + improved_wall + has_phone +head_sec+a26a+a26b,  data =baseline_complete ,method = "nearest",distance = "logit")
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)

### this is the matched data -  we now need to extract these ids from the endline and stack base and endline to do a dif-in-dif

matched.endline <- endline[endline$hhid %in% matched.baseline$hhid,]
matched.baseline$time <- 0
matched.endline$time <- 1

matched.baseline <- matched.baseline[c("information","deliberation","time",baseline_outcomes[i] )]
names(matched.baseline) <- c("information","deliberation","time",outcomes[i] )


matched.dta_long <- rbind(matched.endline[c("information","deliberation","time", outcomes[i])], matched.baseline[ c("information","deliberation","time",outcomes[i] )])

ols <- lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=matched.dta_long)
res <- coef(summary(lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=matched.dta_long) ))
conf <- confint(lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=matched.dta_long) )

df_matcher[,1,i] <- c(res[6,1],res[6,2],res[6,4], conf[6,1], conf[6,2], nobs(ols))

####matching for deliberation
nearest.match <- matchit(formula = deliberation ~ hhsize + femhead + agehead + log_farmsize + ironroof + improved_wall + has_phone +head_sec+a26a+a26b,  data =baseline_complete ,method = "nearest",distance = "logit")
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)

### this is the matched data -  we now need to extract these ids from the endline and stack base and endline to do a dif-in-dif

matched.endline <- endline[endline$hhid %in% matched.baseline$hhid,]
matched.baseline$time <- 0
matched.endline$time <- 1

matched.baseline <- matched.baseline[c("information","deliberation","time",baseline_outcomes[i] )]
names(matched.baseline) <- c("information","deliberation","time",outcomes[i] )


matched.dta_long <- rbind(matched.endline[c("information","deliberation","time", outcomes[i])], matched.baseline[ c("information","deliberation","time",outcomes[i] )])
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=matched.dta_long)
res <- coef(summary(lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=matched.dta_long) ))
conf <- confint(lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=matched.dta_long))

df_matcher[,2,i] <- c(res[7,1],res[7,2],res[7,4], conf[7,1], conf[7,2], nobs(ols))

####matching for interaction
nearest.match <- matchit(formula = (information*deliberation) ~ hhsize + femhead + agehead + log_farmsize + ironroof + improved_wall + has_phone +head_sec+a26a+a26b,  data =baseline_complete ,method = "nearest",distance = "logit")
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)

### this is the matched data -  we now need to extract these ids from the endline and stack base and endline to do a dif-in-dif

matched.endline <- endline[endline$hhid %in% matched.baseline$hhid,]
matched.baseline$time <- 0
matched.endline$time <- 1

matched.baseline <- matched.baseline[c("information","deliberation","time",baseline_outcomes[i] )]
names(matched.baseline) <- c("information","deliberation","time",outcomes[i] )


matched.dta_long <- rbind(matched.endline[c("information","deliberation","time", outcomes[i])], matched.baseline[ c("information","deliberation","time",outcomes[i] )])
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=matched.dta_long)
res <- coef(summary(lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=matched.dta_long) ))
conf <- confint(lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=matched.dta_long) )
df_matcher[,3,i] <- c(res[8,1],res[8,2],res[8,4], conf[8,1], conf[8,2], nobs(ols))

}

### create data.frame to plot - make sure you get correct i's for the indices; last one is overall index
d_plot <- data.frame(rbind(df_matcher[c(1,4,5),1,7],df_matcher[c(1,4,5),2,7],df_matcher[c(1,4,5),3,7]))
d_plot <- rbind(d_plot,data.frame(rbind(df_matcher[c(1,4,5),1,12],df_matcher[c(1,4,5),2,12],df_matcher[c(1,4,5),3,12])))
d_plot <- rbind(d_plot,data.frame(rbind(df_matcher[c(1,4,5),1,19],df_matcher[c(1,4,5),2,19],df_matcher[c(1,4,5),3,19])))
d_plot <- rbind(d_plot, data.frame(rbind(c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA))))
d_plot <- rbind(d_plot,data.frame(rbind(df_matcher[c(1,4,5),1,20],df_matcher[c(1,4,5),2,20],df_matcher[c(1,4,5),3,20])))


names(d_plot) <- c("y","ylo","yhi")
rep(1:4, times=3, each=3)
d_plot$x <- rep(c("agricuture","infrastructure","health","","index"), each=3)
d_plot$grp <- rep(c("info","delib","both"), times=5)
d_plot$grp <-  factor(d_plot$grp , levels=c("info","delib","both"))
d_plot$x <-  factor(d_plot$x, levels=rev((c("agricuture","infrastructure","health","","index"))))
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/impact_summary.png", units="px", height=3200, width= 6400, res=600)
credplot.gg(d_plot,'SDs','',levels(d_plot$x),.2)
dev.off()
 # d is a data frame with 4 columns
 # d$x gives variable names
 # d$y gives center point
 # d$ylo gives lower limits
 # d$yhi gives upper limits


