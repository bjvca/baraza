rm(list=ls())
library(dplyr)
library(ggplot2)
library(MatchIt) 
library(multiwayvcov)
library(plm)
library(lmtest)
library(clubSandwich)
library(moments)
set.seed(123456789) #not needed for final version?

# takes raw data (baseline and endline), makes it anonymous and puts in into the data/public folder, ready to be analysed by the code chucks below
source("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/anonyize.R")
endline$a21 <- as.character(endline$region)
endline$region <- NULL
### EDITS SHOULD BE MADE HERE FOR FINAL VERSION###########################################################################################
### to do in final version: #
### - check skewness before taking IHS.
### - variables in endline that were transformed should also be transformed in baseline
### - number of days sick can only be determined if all data is in, as we need to establish the max number of sick household members
### there should be no duplicates in this dataset
endline <- endline[!duplicated(endline$hhid),]

#endline$district_baraza[endline$subcounty=="CEGERE"] <- 1
#endline$deliberation[endline$subcounty=="INOMO"] <- 1
#endline$information[endline$subcounty=="INOMO"] <- 0
endline$a21[sample(1:dim(endline)[1],50)] <- "Central"
endline$a21 <- as.factor(endline$a21)

### also search and delete below:

###################################################################
#baseline$district_baraza[baseline$a23 == "KAKOBA"] <- 1
####################################################################
##########################################################################################################################################






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

################################################################## end of funtions declarations

### for the mock report, I use a dummy endline - I read in a dummy endline of 3 households just to get the correct variable names
#endline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/endline.csv")[10:403]

#### I then merge with the sampling list to basically create an empty endline dataset
#### and merge in the treatments
#list <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/sampling_list_hh.csv")[c("hhid","a21","a22","a23")]
#endline <- merge(list, endline, by="hhid", all.x=T)
#treats <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/final_list_5.csv")
#endline <- merge(treats, endline, by.x=c("district","subcounty"), by.y=c("a22","a23"))



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
baseline_desc <- baseline
baseline_matching <- merge(baseline,treats, by.x=c("a22","a23"), by.y=c("district","subcounty"))

#define endline variables -
#endline$baraza.B3 <- 0
endline$baraza.B3 <- endline$baraza.B3 ==1 |  endline$baraza.B3.3 ==1
#endline$inputs <- 0
endline$inputs <- as.numeric(endline$baraza.B1==1 | endline$baraza.B1.5==1) 
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
endline$baraza.D1.2 <- sample(baseline$tot_sick[!is.na(baseline$tot_sick)] ,dim(endline)[1]) 
endline$baraza.D1.2[is.na(endline$baraza.D1.2)] <- 0

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

##endline$baraza.E1.6 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e14, na.rm=T))
##endline$baraza.E1.10 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e22, na.rm=T))
##endline$baraza.E1.13 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e32, na.rm=T))
##endline$baraza.E1.18 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$e45, na.rm=T))
##ag


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
baseline <- FW_index(c("base_unprotected","c12source", "qc15","c10","a6"),revcols=c(1,2,3,5),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_infra_index'
baseline_matching <- FW_index(c("base_unprotected","c12source", "qc15","c10","a6"),revcols=c(1,2,3,5),data=baseline_matching)
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
#endline <- FW_index(c("n_children","baraza.E1","baraza.E1.4","baraza.E1.6","baraza.E1.10","baraza.E1.13","baraza.E1.18"),revcols=c(2),data=endline)
#names(endline)[names(endline) == 'index'] <- 'education_index'
#baseline <- FW_index(c("base_n_children","e5","e12", "e14","e22","e32","e45"),revcols=c(2),data=baseline)
#names(baseline)[names(baseline) == 'index'] <- 'base_education_index'
#baseline_matching <- FW_index(c("base_n_children","e5","e12", "e14","e22","e32","e45"),revcols=c(2),data=baseline_matching)
#names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_education_index'

endline <- FW_index(c("n_children","baraza.E5","baraza.E12","baraza.E14","baraza.E22"),revcols=c(2),data=endline)
names(endline)[names(endline) == 'index'] <- 'education_index'
baseline <- FW_index(c("base_n_children","e5","e12", "e14","e22"),revcols=c(2),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_education_index'
baseline_matching <- FW_index(c("base_n_children","e5","e12", "e14","e22"),revcols=c(2),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_education_index'


#20 make an index of indices
endline <- FW_index(c("ag_index","infra_index","health_index","education_index"),data=endline)
names(endline)[names(endline) == 'index'] <- 'pub_service_index'
baseline <- FW_index(c("base_ag_index","base_infra_index","base_health_index","base_education_index"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_pub_service_index'
baseline_matching <- FW_index(c("base_ag_index","base_infra_index","base_health_index","base_education_index"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_pub_service_index'



outcomes <- c("baraza.B2","baraza.B3","baraza.B4.1","inputs","baraza.B5.2","baraza.B5.3","ag_index","unprotected", "baraza.C1.2", "baraza.C1.3","baraza.C2.3","baraza.A6","infra_index","baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1",  "baraza.D4.6","baraza.D6","health_index","n_children","baraza.E5","baraza.E12","baraza.E14","baraza.E22","baraza.E32","baraza.E45","education_index", "pub_service_index")
baseline_outcomes <- c("b21","b31","b44","base_inputs","b5144","b5146","base_ag_index","base_unprotected","c12source", "qc15","c10","a6","base_infra_index","pub_health_access","maternal_health_access","d31","d43","d11","wait_time","d61","base_health_index","base_n_children","e5","e12", "e14","e22","e32","e45","base_education_index","base_pub_service_index")

##      outcomes            baseline_outcomes    

## agriculture   
## [1,] "baraza.B2"         "b21"                   	Was visited by extension officer at home (yes/no)
## [2,] "baraza.B3"         "b31"                   	Visited training or demonstration site (yes/no)
## [3,] "baraza.B4.1"       "b44"                   	NAADS or OWC in village (yes/no)
## [4,] "inputs"            "base_inputs"           	Uses modern inputs (improved seed or fertilizer) (yes/no)
## [5,] "baraza.B5.2"       "b5144"                 	Support in marketing from village procurement committe (yes/no)
## [6,] "baraza.B5.3"       "b5146"                 	Support in marketing from cooperative (yes/no)
## [7,] "ag_index"          "base_ag_index"         

##infrastructure
## [8,] "unprotected"       "base_unprotected"      	Household uses unprotected water source during dry season (yes/no)
## [9,] "baraza.C1.2"       "c12source"             	Distance to water source
##[10,] "baraza.C1.3"       "qc15"                  	Average waiting time at source (min)
##[11,] "baraza.C2.3"       "c10"             		Is there a Water User Committee in this village? (yes/no)
##[12,] "baraza.A6"         "a6"			Distance to nearest all weather road (km)             
##[13,] "infra_index"       "base_infra_index"   

##health   
##[14,] "baraza.D2"         "pub_health_access"     	Seek treatment for fever in public health facility (1=yes) 
##[15,] "baraza.D2.4"       "maternal_health_access"	Go to public health facility to give birth (1=yes)  
##[16,] "baraza.D3"         "d31"                   	Is there a VHT in village? (1=yes) 
##[17,] "baraza.D4.2"       "d43"                   	Distance to nearest govt health facility (km) 
##[18,] "baraza.D1"         "d11"                   	Were days work/school missed due to illness? (1=yes)  
##[19,] "baraza.D4.6"       "wait_time"           	Waiting time before being attended (min)  
#20	"baraza.D6" 		"d61"			Has visited traditional health practitioner? (1=yes)  
   
##[21,] "health_index"      "base_health_index"     

##educations
##[22,] "n_children"	    "base_n_children"		Number of children in UPS or USE 
##23"baraza.E1"		e5				Distance to public school (km)  
##24"baraza.E1.4" "e12" 				Has complete boundary fence (1=yes) 
##25"baraza.E1.6" "e14" 				Has water facility (1=yes) 
##26"baraza.E1.10 "e22" 				Has School Management Committee (1=yes)  
##27"baraza.E1.13 "e32" 				Is informed about School Management Committee (1=yes)  
##28"baraza.E1.18 "e45" 				Inspectors visited schools (1=yes)  

##29 "education_index"      "base_education_index"  

##30 "pub_service_index" "base_pub_service_index"


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


###init arrays to store results
df_ols <- array(NA,dim=c(6,4,length(outcomes)))
df_ancova <- array(NA,dim=c(6,4,length(outcomes)))
df_dif_in_dif <- array(NA,dim=c(6,4,length(outcomes)))
df_matcher <- array(NA,dim=c(6,4,length(outcomes)))
df_averages <- array(NA,dim=c(2,length(outcomes)))

for (i in 1:length(outcomes)) {
#print(i)
# i <- 1

df_averages[1,i] <- mean(as.matrix(endline[outcomes[i]]), na.rm=T)
df_averages[2,i] <- sd(as.matrix(endline[outcomes[i]]), na.rm=T)

### simple difference and adjust se for clustered treatment assignment
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+a21",sep="~")), data=endline[endline$district_baraza == 0,]) 
vcov_cluster <- vcovCR(ols, cluster = endline$clusterID[endline$district_baraza == 0], type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,4], conf[2,4],conf[2,5], nobs(ols))
df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,4], conf[3,4],conf[3,5], nobs(ols))
df_ols[,1,i] <- c(res[7,1],res[7,2],res[7,4], conf[7,4],conf[7,5], nobs(ols))

ols <- lm(as.formula(paste(outcomes[i],"district_baraza+a21",sep="~")), data=endline[(endline$information == 1 & endline$deliberation==1) | endline$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = endline$clusterID2[(endline$information == 1 & endline$deliberation==1) | endline$district_baraza == 1 ], type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
df_ols[,4,i] <- c(res[2,1],res[2,2],res[2,4], conf[2,4],conf[2,5], nobs(ols))

##ancova
## merge in baseline

ols <- lm(as.formula(paste(paste(outcomes[i],"information*deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta[dta$district_baraza == 0,]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID[dta$district_baraza == 0], type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_ancova[,2,i] <- c(res[2,1],res[2,2],res[2,4], conf[2,4],conf[2,5], nobs(ols))
df_ancova[,3,i] <- c(res[3,1],res[3,2],res[3,4], conf[3,4],conf[3,5], nobs(ols))
df_ancova[,1,i] <- c(res[7,1],res[7,2],res[7,4], conf[7,4],conf[7,5], nobs(ols))


ols <- lm(as.formula(paste(paste(outcomes[i],"district_baraza+a21",sep="~"),baseline_outcomes[i],sep="+")), data=dta[(dta$information == 1 & dta$deliberation==1) | dta$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = dta$clusterID2[(dta$information == 1 & dta$deliberation==1) | dta$district_baraza == 1 ], type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
df_ancova[,4,i] <- c(res[2,1],res[2,2],res[2,4], conf[2,4],conf[2,5], nobs(ols))


## dif-in-dif
ols <- lm(as.formula(paste(outcomes[i],"information*deliberation*time",sep="~")), data=dta_long[dta_long$district_baraza == 0,])
vcov_cluster <- vcovCR(ols, cluster = dta_long$clusterID[dta_long$district_baraza == 0], type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)

df_dif_in_dif[,2,i] <- c(res[6,1],res[6,2],res[6,4], conf[6,4], conf[6,5], nobs(ols))
df_dif_in_dif[,3,i] <- c(res[7,1],res[7,2],res[7,4], conf[7,4], conf[7,5], nobs(ols))
df_dif_in_dif[,1,i] <- c(res[8,1],res[8,2],res[8,4], conf[8,4], conf[8,5], nobs(ols))


ols <- lm(as.formula(paste(outcomes[i],"district_baraza*time",sep="~")), data=dta_long[(dta_long$information == 1 & dta_long$deliberation==1) | dta_long$district_baraza == 1 ,]) 
vcov_cluster <- vcovCR(ols, cluster = dta_long$clusterID2[(dta_long$information == 1 & dta_long$deliberation==1) | dta_long$district_baraza == 1 ], type = "CR0")
res <- coeftest(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
df_dif_in_dif[,4,i] <- c(res[4,1],res[4,2],res[4,4], conf[4,4], conf[4,5], nobs(ols))


}

### create data.frame to plot - make sure you get correct i's for the indices; last one is overall index
d_plot <- data.frame(rbind(df_ancova[c(1,4,5),1,7],df_ancova[c(1,4,5),2,7],df_ancova[c(1,4,5),3,7],df_ancova[c(1,4,5),4,7]))
d_plot <- rbind(d_plot,data.frame(rbind(df_ancova[c(1,4,5),1,13],df_ancova[c(1,4,5),2,13],df_ancova[c(1,4,5),3,13],df_ancova[c(1,4,5),4,13])))
d_plot <- rbind(d_plot,data.frame(rbind(df_ancova[c(1,4,5),1,21],df_ancova[c(1,4,5),2,21],df_ancova[c(1,4,5),3,21],df_ancova[c(1,4,5),4,21])))
d_plot <- rbind(d_plot,data.frame(rbind(df_ancova[c(1,4,5),1,29],df_ancova[c(1,4,5),2,29],df_ancova[c(1,4,5),3,29],df_ancova[c(1,4,5),4,29])))
d_plot <- rbind(d_plot, data.frame(rbind(c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA))))
d_plot <- rbind(d_plot,data.frame(rbind(df_ancova[c(1,4,5),1,30],df_ancova[c(1,4,5),2,30],df_ancova[c(1,4,5),3,30],df_ancova[c(1,4,5),4,30])))


names(d_plot) <- c("y","ylo","yhi")

d_plot$x <- rep(c("agricuture","infrastructure","health","education","","index"), each=4)
d_plot$grp <- rep(c("sc baraza","info","delib","level"), times=6)
d_plot$grp <-  factor(d_plot$grp , levels=c("sc baraza","info","delib","level"))
d_plot$x <-  factor(d_plot$x, levels=rev((c("agricuture","infrastructure","health","education","","index"))))
png("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/figure/impact_summary_ancova.png", units="px", height=3200, width= 6400, res=600)
credplot.gg(d_plot,'SDs','',levels(d_plot$x),.3)
dev.off()

