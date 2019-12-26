rm(list=ls())
library(dplyr)
library(ggplot2)
library(MatchIt) 
library(multiwayvcov)
library(plm)
library(lmtest)
library(clubSandwich)
library(moments)
set.seed(54321) #not needed for final version?

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

# takes raw data (baseline and endline), makes it anonymous and puts in into the data/public folder, ready to be analysed by the code chucks below
#source("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/anonyize.R")
### for the mock report, I use a dummy endline - I read in a dummy endline of 3 households just to get the correct variable names
endline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/endline.csv")[10:403]

### I then merge with the sampling list to basically create an empty endline dataset
### and merge in the treatments
list <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/sampling_list_hh.csv")[c("hhid","a21","a22","a23")]
endline <- merge(list, endline, by="hhid", all.x=T)
treats <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/final_list_5.csv")
endline <- merge(treats, endline, by.x=c("district","subcounty"), by.y=c("a22","a23"))



## baseline not needed in this first section, but used to generate fake data
baseline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/baseline.csv")
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
baseline$base_unprotected <- as.numeric(( baseline$c11a %in%  c("Rain water","Surface water","Tube well or borehole","Unprotected dug well","Unprotected spring"))    )
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
endline$baraza.A6 <- sample(baseline$a6,dim(endline)[1]) ### this needs to be inverse hypersine transformed and trimmed in final version

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
endline$baraza.D4.6 <- sample(baseline$wait_time ,dim(endline)[1]) 
endline$baraza.D6  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$d61, na.rm=T))

##edu
endline$n_children <- rowSums(cbind(endline$baraza.E1.2,endline$baraza.E2.1), na.rm=T) 
endline$n_children <- sample(baseline$base_n_children ,dim(endline)[1]) 
endline$baraza.E1  <- sample(baseline$e5 ,dim(endline)[1]) 
endline$baraza.E1.4 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e12, na.rm=T))
endline$baraza.E1.6 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e14, na.rm=T))
endline$baraza.E1.10 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e22, na.rm=T))
endline$baraza.E1.13 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e32, na.rm=T))
endline$baraza.E1.18 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e45, na.rm=T))
##ag

#1 access to extension, 
#2 visit extension officer or demo site, 
#2 farmer in naads supported group, 
#4 uses inputs (seed or fert), 
#5 support in marketing from Village procurement committe/Village farmers forum/Village farmers forum executive; 
#6 support in marketing from Cooperative

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
endline <- FW_index(c("unprotected", "baraza.C1.2", "baraza.C1.3","baraza.C2.3","baraza.A6"),revcols=c(1,2,3),data=endline)
names(endline)[names(endline) == 'index'] <- 'infra_index'
baseline <- FW_index(c("base_unprotected","c12source", "qc15","c10","a6"),revcols=c(1,2,3),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_infra_index'
baseline_matching <- FW_index(c("base_unprotected","c12source", "qc15","c10","a6"),revcols=c(1,2,3),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_infra_index'

##make a health index
#13 pub health access
#14 maternal health acess
#15 is there a VHT?
#16 distance to gvt health facility
##17 has been sick
###18 number of days sick
## 19 index
endline <- FW_index(c("baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1",  "baraza.D4.6", "baraza.D6"),revcols=c(4,5,6,7),data=endline)
names(endline)[names(endline) == 'index'] <- 'health_index'
baseline <- FW_index(c("pub_health_access","maternal_health_access","d31","d43","d11","wait_time","d61"),revcols=c(4,5,6,7),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_health_index'
baseline_matching <- FW_index(c("pub_health_access","maternal_health_access","d31","d43","d11","wait_time","d61"),revcols=c(4,5,6,7),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_health_index'

###make and education index
endline <- FW_index(c("n_children","baraza.E1","baraza.E1.4","baraza.E1.6","baraza.E1.10","baraza.E1.13","baraza.E1.18"),revcols=c(2),data=endline)
names(endline)[names(endline) == 'index'] <- 'education_index'
baseline <- FW_index(c("base_n_children","e5","e12", "e14","e22","e32","e45"),revcols=c(2),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_education_index'
baseline_matching <- FW_index(c("base_n_children","e5","e12", "e14","e22","e32","e45"),revcols=c(2),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_education_index'

#20 make an index of indices
endline <- FW_index(c("ag_index","infra_index","health_index","education_index"),data=endline)
names(endline)[names(endline) == 'index'] <- 'pub_service_index'
baseline <- FW_index(c("base_ag_index","base_infra_index","base_health_index","base_education_index"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_pub_service_index'
baseline_matching <- FW_index(c("base_ag_index","base_infra_index","base_health_index","base_education_index"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_pub_service_index'



outcomes <- c("baraza.B2","baraza.B3","baraza.B4.1","inputs","baraza.B5.2","baraza.B5.3","ag_index","unprotected", "baraza.C1.2", "baraza.C1.3","baraza.C2.3","baraza.A6","infra_index","baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1",  "baraza.D4.6","baraza.D6","health_index","n_children","baraza.E1","baraza.E1.4","baraza.E1.6","baraza.E1.10","baraza.E1.13","baraza.E1.18","education_index", "pub_service_index")
baseline_outcomes <- c("b21","b31","b44","base_inputs","b5144","b5146","base_ag_index","base_unprotected","c12source", "qc15","c10","a6","base_infra_index","pub_health_access","maternal_health_access","d31","d43","d11","wait_time","d61","base_health_index","base_n_children","e5","e12", "e14","e22","e32","e45","base_education_index","base_pub_service_index")

##      outcomes            baseline_outcomes    

## agriculture   
## [1,] "baraza.B2"         "b21"                   	"Was visited by an extension officer at home"
## [2,] "baraza.B3"         "b31"                   	"Visited extension office or demo site"
## [3,] "baraza.B4.1"       "b44"                   	"NAADS/OWC in village"
## [4,] "inputs"            "base_inputs"           	"uses modern inputs"
## [5,] "baraza.B5.2"       "b5144"                 	"marketing help village committee"
## [6,] "baraza.B5.3"       "b5146"                 	"marketing help cooperative"
## [7,] "ag_index"          "base_ag_index"         

##infrastructure
## [8,] "unprotected"       "base_unprotected"      	"Household uses unprotected water source during dry season"
## [9,] "baraza.C1.2"       "c12source"             	"How far in Km is this water source from your household?"
##[10,] "baraza.C1.3"       "qc15"                  	"On average, how long do you have to wait to collect water during the dry season? (mins)"
##[11,] "baraza.C2.3"       "c10"             		"Is there a Water User Committee in this village?"
##[12,] "baraza.A6"         "a6"			"Distance to nearest all weather road (km)"             
##[13,] "infra_index"       "base_infra_index"   

##health   
##[14,] "baraza.D2"         "pub_health_access"     	"Seek treatment for fever in public health facility"
##[15,] "baraza.D2.4"       "maternal_health_access"	"Go to public health facility to give birth"
##[16,] "baraza.D3"         "d31"                   	"Is there a VHT in village"
##[17,] "baraza.D4.2"       "d43"                   	"Distance to nearest govt health facility"
##[18,] "baraza.D1"         "d11"                   	"Were any household members unable to work or go to school due to an illness in the past one year?"
##[19,] "baraza.D4.6"       "wait_time"           	"How long did you have to wait before being attended (in min)"
#20	"baraza.D6" 		"d61"			"Have you visited a traditional health practitioner in the last year?"
   
##[21,] "health_index"      "base_health_index"     

##educations
##[22,] "n_children"	    "base_n_children"		"Number of children in UPS or USE"
##23"baraza.E1"		e5				"Distance to UPE school
##24"baraza.E1.4" "e12" 			"Complete boundary fence
##25"baraza.E1.6" "e14" 			"Has water facility
##26"baraza.E1.10 "e22" 			"Has SMC
##27"baraza.E1.13 "e32" 			"Informed about SMC
##28"baraza.E1.18 "e45" 				"Inspectors visited schools

##29 "education_index"      "base_education_index"  

##30 "pub_service_index" "base_pub_service_index"


#create unique ID for clustering based on district and subcounty
endline <- endline %>%  mutate(clusterID = group_indices(., district, subcounty))
endline <- endline %>%  mutate(clusterID2 = group_indices(., district))

###drop all district barazas:
endline <- subset(endline, district_baraza == 0)

endline <- merge(endline, baseline[c("hhid","treat")], by.x=c("hhid"), by.y=c("hhid"))
endline$information_planned <- 0
endline$information_planned[endline$treat == "info" | endline$treat == "scbza"] <- 1
endline$deliberation_planned <- 0
endline$deliberation_planned[endline$treat == "delib" | endline$treat == "scbza"] <- 1

endline_info <- subset(endline, information == 0)

##init arrays to store results
df_balance <- array(NA,dim=c(6,2,length(outcomes)))

for (i in 1:length(outcomes)) {
## simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information_planned+a21",sep="~")), data=endline_info) 
vcov_cluster <- vcovCR(ols, cluster = endline_info$clusterID, type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_balance[,1,i] <- c(res[2,1],res[2,2],res[2,4], conf[2,4],conf[2,5], nobs(ols))

}

#for deliberation

endline_delib <- subset(endline, deliberation == 0)


for (i in 1:length(outcomes)) {
### simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"deliberation_planned+a21",sep="~")), data=endline_delib) 
vcov_cluster <- vcovCR(ols, cluster = endline_delib$clusterID, type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_balance[,2,i] <- c(res[3,1],res[3,2],res[3,4], conf[3,4],conf[3,5], nobs(ols))
}

#
