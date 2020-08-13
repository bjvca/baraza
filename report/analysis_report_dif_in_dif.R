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

### this is executed in the /report subdirectory, need to ..
path <- strsplit(getwd(), "/report")[[1]]

### set this switch to TRUE if you want to produce a final report - this will save results matrices in a static directory
final_verion_swith <- TRUE
RI_conf_switch <- FALSE

glob_repli <- 1000
glob_sig <- c(.025,.975) ### 5 percent conf intervals

# takes raw data (baseline and endline), makes it anonymous and puts in into the data/public folder, ready to be analysed by the code chucks below
#source("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/cleaning.R")
#source("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/anonyize.R")
endline <- read.csv(paste(path,"data/public/endline.csv", sep="/"), stringsAsFactors = TRUE)
endline$a21 <- as.character(endline$region)
endline$region <- NULL
### EDITS SHOULD BE MADE HERE FOR FINAL VERSION###########################################################################################
### to do in final version: #
### - check skewness before taking IHS.
### - variables in endline that were transformed should also be transformed in baseline
### - number of days sick can only be determined if all data is in, as we need to establish the max number of sick household members
### there should be no duplicates in this dataset
endline <- endline[!duplicated(endline$hhid),]

endline$a21 <- as.factor(endline$a21)

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

#### for the mock report, I use a dummy endline - I read in a dummy endline of 3 households just to get the correct variable names
#endline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/endline.csv", stringsAsFactors = TRUE)[10:403]

#### I then merge with the sampling list to basically create an empty endline dataset
#### and merge in the treatments
#list <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/sampling_list_hh.csv", stringsAsFactors = TRUE)[c("hhid","a21","a22","a23")]
#endline <- merge(list, endline, by="hhid", all.x=T)
treats <- read.csv(paste(path,"questionnaire/final_list_5.csv", sep ="/"), stringsAsFactors = TRUE)
#endline <- merge(treats, endline, by.x=c("district","subcounty"), by.y=c("a22","a23"))



## baseline not needed in this first section, but used to generate fake data
baseline <- read.csv(paste(path,"data/public/baseline.csv",sep="/"), stringsAsFactors = TRUE)
baseline$a23 <- as.character(baseline$a23)
baseline$a23[baseline$a23 == "LUWERO TC"] <- "LUWERO_TC"
baseline$a23[baseline$a23 == "SEMBABULE TC"] <- "SEMBABULE_TC"
baseline$a23[baseline$a23 == "RAKAI TC"] <- "RAKAI_TC"
baseline$a23[baseline$a23 == "NTUSI"] <- "NTUUSI"


baseline$b21 <-  as.numeric(baseline$b21=="Yes")
baseline$b31 <-  as.numeric(baseline$b31=="Yes")
baseline$b41 <-  as.numeric(baseline$b41=="Yes")
baseline$b41[is.na(baseline$b41)] <- 0
baseline$b44 <-  as.numeric(baseline$b44=="Yes")
baseline$b44[is.na(baseline$b44)] <- 0
baseline$base_inputs <- as.numeric(baseline$used_seed=="Yes" | baseline$used_fert=="Yes")
baseline$used_seed <- baseline$used_seed=="Yes"
baseline$used_fert <- baseline$used_fert=="Yes"
baseline$used_chem <- baseline$used_chem=="Yes"
baseline$b5144 <- as.numeric(baseline$b5144=="Yes")
baseline$b5146 <- as.numeric(baseline$b5146=="Yes")
##use of unprotected water sources in dry season
###this was changed post registration to follow https://www.who.int/water_sanitation_health/monitoring/jmp2012/key_terms/en/ guidelines on what is considered improved, that also considers rainwater a protected source
baseline$base_unprotected <- (as.numeric(baseline$c11a) %in%  c(10,13,14))
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
baseline$not_work[baseline$d11==0] <- 0 
baseline$not_school[baseline$d11==0] <- 0 
baseline$not_work_school <- baseline$not_work + baseline$not_school

baseline$tot_sick <- log(baseline$tot_sick + sqrt(baseline$tot_sick ^ 2 + 1))
baseline <- trim("tot_sick", baseline)

baseline$not_work_school <- log(baseline$not_work_school + sqrt(baseline$not_work_school ^ 2 + 1))
baseline <- trim("not_work_school", baseline)

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

baseline$e13 <- rowSums(cbind(as.numeric(baseline$e13upe == "Yes") , as.numeric(baseline$e13use == "Yes")), na.rm=T) > 0
baseline$e13[is.na(baseline$e13upe) & is.na(baseline$e13use)] <- NA

baseline$e14 <- rowSums(cbind(as.numeric(baseline$e14upe == "Yes") , as.numeric(baseline$e14use == "Yes")), na.rm=T) > 0
baseline$e14[is.na(baseline$e14upe) & is.na(baseline$e14use)] <- NA

baseline$e18 <- rowSums(cbind(as.numeric(baseline$e18upe == 1) , as.numeric(baseline$e18use == 1)), na.rm=T) > 0
baseline$e18[is.na(baseline$e18upe) & is.na(baseline$e18use)] <- NA

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

baseline$d411 <- log(baseline$d411 + sqrt(baseline$d411 ^ 2 + 1))
baseline <- trim("d411", baseline)

baseline$ironroof <- as.numeric(baseline$a512 =="Corrugated iron sheets")
baseline$improved_wall <- as.numeric(baseline$a513 %in% c("Mud_bricks_burnt_bricks","Concrete_blocks") )
baseline$head_sec <- as.numeric(baseline$a36) > 15

baseline$f241.LC1.election <- baseline$f241.LC1.election == "Yes"   
baseline$f241.LC3.election <- baseline$f241.LC3.election == "Yes"   
baseline$f241.LC5.election <- baseline$f241.LC5.election == "Yes"   
baseline$f241.Pesidential <- baseline$f241.Pesidential == "Yes"   
baseline$f241.Parliamentary <- baseline$f241.Parliamentary == "Yes"  
baseline$f241.Party.leader <- baseline$f241.Party.leader == "Yes"


baseline$f2301 <- as.numeric(baseline$f2301) %in% c(6,7)
baseline$f2303 <- as.numeric(baseline$f2303) %in% c(1,2,3,6,7)
baseline$f2307  <- as.numeric(baseline$f2307) %in% c(1,3,6,7) 
baseline$f2309 <- as.numeric(baseline$f2309) %in% c(1,2,3,6,7)
baseline$f2310 <- as.numeric(baseline$f2310) %in% c(1,2,3,6,7)

baseline$f21 <- baseline$f21 == "Yes"

baseline$roof <- baseline$a512 %in% c("Corrugated iron sheets", "Tiles")
baseline$wall <- baseline$a513 == "Mud_bricks_burnt_bricks"
baseline$b314 <- baseline$b314 == "Yes"
baseline$b316 <- baseline$b316 == "Yes"
baseline$b316[is.na(baseline$b316)] <- FALSE
baseline$used_livestock_tech <- baseline$used_livestock_tech == "Yes"
baseline$d416 <- baseline$d416 == "Yes"
baseline$d419 <- baseline$d419 == "Yes"

baseline$b320 <- baseline$b320 == "Decided by extension agents/forum members without any consultation"
baseline$qc16 <- (baseline$qc16 == "Satisfied"| baseline$qc16 == "Very satisfied")
baseline$d420 <- (baseline$d420 == "Satisfied" | baseline$d420 == "Very satisfied")
baseline$c4 <- baseline$c4 %in% c("Boil","Use chlorine/bleach")
baseline$c11 <- baseline$c11 == "Yes"
baseline$c11[is.na(baseline$c11 )] <- FALSE 
baseline$c13 <- baseline$c13 == "Yes"
baseline$c13[is.na(baseline$c13 )] <- FALSE
baseline$d32 <- baseline$d32 == "Yes"
baseline$d32[is.na(baseline$d32 )] <- FALSE
baseline$d315 <- baseline$d315 == "Yes"
baseline$d315[is.na(baseline$d315 )] <- FALSE
baseline$base_doctor <- baseline$d49 %in% c("Doctor","In-charge")
baseline$base_doctor[is.na(baseline$d49)] <- NA
baseline$base_paid_health <- baseline$d413 == "Yes"
baseline$base_paid_health[is.na(baseline$d413)] <- NA
baseline$d426 <- (baseline$d426 == "Yes")

endline$baraza.D4.7 <- as.numeric(as.character(endline$baraza.D4.7))
endline$baraza.D4.7[endline$baraza.D4.7 == 999] <- NA
endline$baraza.D4.7 <-  log(endline$baraza.D4.7 + sqrt(endline$baraza.D4.7 ^ 2 + 1))
endline <- trim("baraza.D4.7",endline)


baseline_desc <- baseline
baseline_matching <- merge(baseline,treats, by.x=c("a22","a23"), by.y=c("district","subcounty"))

#define endline variables -
#endline$baraza.B3 <- 0
endline$baraza.B3 <- endline$baraza.B3 ==1 |  endline$baraza.B3.3 ==1
#endline$inputs <- 0
endline$inputs <- as.numeric(endline$baraza.B1==1 | endline$baraza.B1.5==1) 
endline$baraza.B1 <- endline$baraza.B1==1
endline$baraza.B1.5 <- endline$baraza.B1.5==1
endline$baraza.B1.9 <- endline$baraza.B1.9==1
endline$baraza.B1.13 <- endline$baraza.B1.13==1
endline$baraza.B3.4 <- endline$baraza.B3.4==1
endline$baraza.B3.5 <- endline$baraza.B3.5==1

###this was changed post registration to follow https://www.who.int/water_sanitation_health/monitoring/jmp2012/key_terms/en/ guidelines on what is considered improved, that also considers rainwater a protected source
#baseline$base_unprotected <- as.numeric(( baseline$c11a %in%  c("Surface water","Bottled water","Cart with small tank","Unprotected dug well","Unprotected spring","Tanker truck"))    )
### is there are water committee
endline$unprotected <- (as.numeric(endline$baraza.C1 %in% c(5,7,11)) )

### here we simulate endline variables - remove if endline data is in
#endline$baraza.B2  <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b21 == 1, na.rm=T))
endline$baraza.B2  <- endline$baraza.B2  == 1
###visits
### here we simulate endline variables - remove if endline data is in
#endline$baraza.B3 <- rbinom(n=dim(endline)[1],size=1,prob=mean(baseline$b31 == 1, na.rm=T))
###naads in village
### here we simulate endline variables - remove if endline data is in
endline$baraza.B4 <- endline$baraza.B4 == 1
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

members <- paste(paste("baraza.labour", 1:15, sep="."),".D1.2", sep=".")

endline[members] <- lapply(endline[members], function(x) as.numeric(as.character(x)) )
endline[members] <- lapply(endline[members], function(x) replace(x, x == 999,NA) )
endline$baraza.D1.2 <- rowSums(endline[members], na.rm=T)

members <- paste(paste("baraza.labour", 1:15, sep="."),".D1.3", sep=".")

endline[members] <- lapply(endline[members], function(x) as.numeric(as.character(x)) )
endline[members] <- lapply(endline[members], function(x) replace(x, x == 999,NA) )
endline$baraza.D1.3 <- rowSums(endline[members], na.rm=T)


endline$baraza.D1.2 <- log(endline$baraza.D1.2 + sqrt(endline$baraza.D1.2 ^ 2 + 1))
endline <- trim("baraza.D1.2", endline)

endline$baraza.D1.3 <- log(endline$baraza.D1.3 + sqrt(endline$baraza.D1.3 ^ 2 + 1))
endline <- trim("baraza.D1.3", endline)

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

endline$baraza.E13 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.5)) == 1 , as.numeric(as.character(endline$baraza.E2.5))==1), na.rm=T) > 0
endline$baraza.E13[is.na(as.numeric(as.character(endline$baraza.E1.5)) ) & is.na(as.numeric(as.character(endline$baraza.E2.5)) )] <- NA


endline$baraza.E14 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.6)) == 1 , as.numeric(as.character(endline$baraza.E2.6))==1), na.rm=T) > 0
endline$baraza.E14[is.na(as.numeric(as.character(endline$baraza.E1.6)) ) & is.na(as.numeric(as.character(endline$baraza.E2.6)) )] <- NA

endline$baraza.E18 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.9)) == 1 , as.numeric(as.character(endline$baraza.E2.9))==1), na.rm=T) > 0
endline$baraza.E18[is.na(as.numeric(as.character(endline$baraza.E1.9)) ) & is.na(as.numeric(as.character(endline$baraza.E2.9)) )] <- NA

endline$baraza.E22 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.10)) == 1 , as.numeric(as.character(endline$baraza.E2.10))==1), na.rm=T) > 0
endline$baraza.E22[is.na(as.numeric(as.character(endline$baraza.E1.10)) ) & is.na(as.numeric(as.character(endline$baraza.E2.10)) )] <- NA

endline$baraza.E32 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.13)) == 1 , as.numeric(as.character(endline$baraza.E2.13))==1), na.rm=T) > 0
endline$baraza.E32[is.na(as.numeric(as.character(endline$baraza.E1.13)) ) & is.na(as.numeric(as.character(endline$baraza.E2.13)) )] <- NA

endline$baraza.E45 <- rowSums(cbind(as.numeric(as.character(endline$baraza.E1.18)) == 1 , as.numeric(as.character(endline$baraza.E2.18))==1), na.rm=T) > 0
endline$baraza.E45[is.na(as.numeric(as.character(endline$baraza.E1.18)) ) & is.na(as.numeric(as.character(endline$baraza.E2.18)) )] <- NA

### assorted outcomes
### type of roof
endline$baraza.roof <- endline$baraza.roof %in% 1:2

endline$baraza.wall <- endline$baraza.wall == 2

endline$seed_OWC <- endline$baraza.B1.6.1 == "True"
endline$baraza.B3.20.3 <- endline$baraza.B3.20.3 =="True"
endline$baraza.D4.11 <- as.numeric(as.character(endline$baraza.D4.11))
endline$baraza.D4.11 <- endline$baraza.D4.11=="1"
endline$baraza.D4.12 <- as.numeric(as.character(endline$baraza.D4.12))
endline$baraza.D4.12 <- endline$baraza.D4.12=="1"
endline$baraza.D4.13 <- as.numeric(as.character(endline$baraza.D4.13))
endline$baraza.D4.13 <- endline$baraza.D4.13 <= 2
endline$baraza.C1.4 <- endline$baraza.C1.4<=2
endline$baraza.C2.1 <- endline$baraza.C2.1 %in%  c(1:2)
endline$baraza.C2.4 <- endline$baraza.C2.4 == 1
endline$baraza.C2.5  <- endline$baraza.C2.5 == 1

endline$baraza.D3.1 <- endline$baraza.D3.1==1
endline$baraza.D3.3 <- endline$baraza.D3.3==1

endline$doctor <- FALSE
endline$doctor <- endline$baraza.D4.5.1 == "True" | endline$baraza.D4.5.7 == "True" 
endline$doctor[endline$baraza.D4.5.1 == "n/a" & endline$baraza.D4.5.7 == "n/a"] <- NA

endline$paid_health <- endline$baraza.D4.10 == 1
endline$paid_health[endline$baraza.D4.10 == "n/a"] <- NA

endline$baraza.D4.14 <- as.numeric(as.character(endline$baraza.D4.14))
endline$baraza.D4.14 <- endline$baraza.D4.14 == 1

#PLEASE DESCRIBE YOUR PARTICIPATION IN DIFFERENT TYPES OF ELECTIONS:
#F1	F1. In this household, are there any members who currently hold any political/traditional positions?
endline$baraza.F1 <- (endline$baraza.F1 == 1)


#F2	F2. In the last elections did you participate in the  LCI elections? 
#F2.1	F2.1 In the last elections did you participate in the  LC3 elections? 
#F2.2	F2.2  In the last elections did you participate in the LC5  elections?   
#F2.3	F2.3 In the last elections did you participate in  the Presidential elections?  
#F2.4	F2.4  In the last elections did you participate in the Parliamentary election?
#F2.5	F2.5 In the last elections did you participate in the  Party leaders elections? 

endline$baraza.part.F2 <- endline$baraza.part.F2 == 1 
endline$baraza.part.F2.1 <- endline$baraza.part.F2.1 == 1 
endline$baraza.part.F2.2 <- endline$baraza.part.F2.2 == 1
endline$baraza.part.F2.3 <- endline$baraza.part.F2.3 == 1
endline$baraza.part.F2.4 <- endline$baraza.part.F2.4 == 1 
endline$baraza.part.F2.5 <- endline$baraza.part.F2.5 == 1


##ag

#1 


# 7 #make an pol index
endline <- FW_index(c("baraza.F1","baraza.part.F2","baraza.part.F2.1","baraza.part.F2.2","baraza.part.F2.3","baraza.part.F2.4","baraza.part.F2.5"),data=endline)
names(endline)[names(endline) == 'index'] <- 'pol_index'
baseline <- FW_index(c("f21","f241.LC1.election", "f241.LC3.election","f241.LC5.election", "f241.Pesidential", "f241.Parliamentary", "f241.Party.leader"  ),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_pol_index'
baseline_matching <- FW_index(c("f21","f241.LC1.election", "f241.LC3.election","f241.LC5.election", "f241.Pesidential", "f241.Parliamentary", "f241.Party.leader"  ),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_pol_index'

##F1.1 When was the last time that you spoke personally with the LC 1 Chairperson, for a reason relating to service provision in agriculture, health, education, water or roads?
##F1.2 When was the last time that you spoke personally with the Subcounty Chief, for a reason relating to service provision in agriculture, health, education, water or roads?
##F1.3 When was the last time that you spoke personally with the Head teacher/ SMC member, for a reason relating to service provision in agriculture, health, education, water or roads?
##F1.4  When was the last time that you spoke personally with the Health Unit Management Committee (HUMC) Member, for a reason relating to service provision in agriculture, health, education, water or roads?
##F1.5  When was the last time that you spoke personally with the Water committee member, for a reason relating to service provision in agriculture, health, education, water or roads?
endline$baraza.F1.1 <- (endline$baraza.F1.1<=2) # last month
endline$baraza.F1.2 <- (endline$baraza.F1.2<=5) # last year 
endline$baraza.F1.3 <- (endline$baraza.F1.3<=4) # last 1/2 year 
endline$baraza.F1.4 <- (endline$baraza.F1.4<=5) # last year 
endline$baraza.F1.5 <- (endline$baraza.F1.5<=5) # last year 

###make a contact index
endline <- FW_index(c("baraza.F1.1", "baraza.F1.2", "baraza.F1.3","baraza.F1.4","baraza.F1.5"),data=endline)
names(endline)[names(endline) == 'index'] <- 'contact_index'
baseline <- FW_index(c("f2301","f2303", "f2307","f2309","f2310"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_contact_index'
baseline_matching <- FW_index(c("f2301","f2303", "f2307","f2309","f2310"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_contact_index'

###priority rankings
#H1	H1.  Access to a drinking water source is a serious problem. (where 1 means you completely disagree and 10 means you completely agree) 
#H2	H2. Drinking water is usually dirty. (where 1 means you completely disagree and 10 means you completely agree)    
#H3	H3.  Access to a government health centre or hospital is a serious problem. (where 1 means you completely disagree and 10 means you completely agree)
#H4	H4  Government health centres or hospitals do not have relevant medicines. (where 1 means you completely disagree and 10 means you completely agree)
#H5	H5 Staff at government health centres or hospitals are rude to patients. (where 1 means you completely disagree and 10 means you completely agree)
#H6	H6.  Medical staff at government health centres or hospitals are often absent. (where 1 means you completely disagree and 10 means you completely agree)
#H7	H7 Access to a government primary school is a serious problem. (where 1 means you completely disagree and 10 means you completely agree)
#H8	H8 Teachers in government schools are often absent. (where 1 means you completely disagree and 10 means you completely agree)
#H9	H9 Children’s learning outcomes in government schools are poor. (where 1 means you completely disagree and 10 means you completely agree)  
#H10	H10.  Availability/ Access to all-weather roads is a serious problem (where 1 means you completely disagree and 10 means you completely agree)
#H11	H11  Agricultural inputs supplied by the government are of poor quality. (where 1 means you completely disagree and 10 means you completely agree)
#H11b	H11  Agricultural inputs supplied by the government are delivered in time. (where 1 means you completely disagree and 10 means you completely agree)
#H12	H12 There is lack of transparency in how farmers are selected to receive agricultural inputs from the government. (where 1 means you completely disagree and 10 means you completely agree)
#H13	H13  Agricultural extension agents rarely visit. (where 1 means you completely disagree and 10 means you completely agree)
#H14	H14.  Agricultural extension agents are not aware of enterprises or agricultural inputs relevant to farmers. (where 1 means you completely disagree and 10 means you completely agree)

endline$baraza.H1[endline$baraza.H1 == 11] <- NA
endline$baraza.H2[endline$baraza.H2 == 11] <- NA
endline$baraza.H3[endline$baraza.H3 == 11] <- NA
endline$baraza.H4[endline$baraza.H4 == 11] <- NA
endline$baraza.H5[endline$baraza.H5 == 11] <- NA
endline$baraza.H6[endline$baraza.H6 == 11] <- NA
endline$baraza.H7[endline$baraza.H7 == 11] <- NA
endline$baraza.H8[endline$baraza.H8 == 11] <- NA
endline$baraza.H9[endline$baraza.H9 == 11] <- NA
endline$baraza.H10[endline$baraza.H10 == 11] <- NA
endline$baraza.H11[endline$baraza.H11 == 11] <- NA
endline$baraza.H12[endline$baraza.H12 == 11] <- NA
endline$baraza.H13[endline$baraza.H13 == 11] <- NA
endline$baraza.H14[endline$baraza.H14 == 11] <- NA 

endline <- FW_index(c("baraza.H1","baraza.H2","baraza.H3","baraza.H4", "baraza.H5",  "baraza.H6", "baraza.H7", "baraza.H8", "baraza.H9", "baraza.H10", "baraza.H11", "baraza.H12", "baraza.H13","baraza.H14"),data=endline)
names(endline)[names(endline) == 'index'] <- 'priority_index'
baseline <- FW_index(c("i1","i2","i3","i4","i5","i6","i7","i8","i9","i10","i11","i12","i13","i14"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_priority_index'
baseline_matching <- FW_index(c("i1","i2","i3","i4","i5","i6","i7","i8","i9","i10","i11","i12","i13","i14"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_priority_index'

### contributions
#G1	G1: Did you ever make any contributions to the school in the last two years? 
#G1.1	G1.1 Did you ever make any contributions to the health centre in the last two years? 
#G1.2	G1.2 Did you ever make any contributions to the road/ bridge in the last two years? 
#G1.3	G1.3 Did you ever make any contributions to the drinking water facility in the last two years? 
#G1.4	G1.4 Did you ever make any contributions to the dam/ irrigation facility in the last two years? 
#G1.5	G1.5 Did you  make any contributions to any other government or community building/ structure in the last two years? 
endline$baraza.G1k <- endline$baraza.G1 %in% c(1,3)
endline$baraza.G1c <- endline$baraza.G1 %in% c(2,3)

endline$baraza.G1.1k <- endline$baraza.G1.1 %in% c(1,3)
endline$baraza.G1.1c <- endline$baraza.G1.1 %in% c(2,3)

endline$baraza.G1.2k <- endline$baraza.G1.2 %in% c(1,3)
endline$baraza.G1.2c <- endline$baraza.G1.2 %in% c(2,3)

endline$baraza.G1.3k <- endline$baraza.G1.3 %in% c(1,3)
endline$baraza.G1.3c <- endline$baraza.G1.3 %in% c(2,3)

endline$baraza.G1.4k <- endline$baraza.G1.4 %in% c(1,3)
endline$baraza.G1.4c <- endline$baraza.G1.4 %in% c(2,3)

endline$baraza.G1.5k <- endline$baraza.G1.5 %in% c(1,3)
endline$baraza.G1.5c <- endline$baraza.G1.5 %in% c(2,3)

####make a contribution  index
endline <- FW_index(c("baraza.G1k","baraza.G1.1k","baraza.G1.2k","baraza.G1.3k","baraza.G1.4k","baraza.G1.5k"),data=endline)
names(endline)[names(endline) == 'index'] <- 'in_kind_index'
baseline <- FW_index(c("cschoolk", "chealthk","croadk","cwaterk",  "cdamk",    "cbuildk"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_in_kind_index'
baseline_matching <- FW_index(c("cschoolk", "chealthk","croadk","cwaterk",  "cdamk",    "cbuildk"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_in_kind_index'

endline <- FW_index(c("baraza.G1c","baraza.G1.1c","baraza.G1.2c","baraza.G1.3c","baraza.G1.4c","baraza.G1.5c"),data=endline)
names(endline)[names(endline) == 'index'] <- 'in_cash_index'
baseline <- FW_index(c("cschoolc", "chealthc","croadc","cwaterc",  "cdamc",    "cbuildc"),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_in_cash_index'
baseline_matching <- FW_index(c("cschoolc", "chealthc","croadc","cwaterc",  "cdamc",    "cbuildc"),data=baseline_matching)
names(baseline_matching)[names(baseline_matching) == 'index'] <- 'base_in_cash_index'


#1# roof	A7. Type of roof of the household  "roof"
#2# wall	A8. Type of wall of household "wall"
#note5	SECTION B: AGRICULTURE
#3# B1		B1. Did the household use inorganic fertilizers (DAP, Urea, NPK, Foliar,TSP, SSP, MOP) during the last 12 months?  "used_fert"
#4# B1.5	B1.5 Did the household use improved seeds during the last 12 months? "used_seed"
#5# seed_OWC	From whom did you buy or obtain these improved seeds? "seed NARO"
#6# B1.9	B1.9 Did the household use any agro-chemicals (pesticides/ herbicides/ fungicides/acaricides) during the last 12 months? "used_chem"
#7 ###B1.13	B1.13 Did the household use improved livestock methods (breeds/feeds/drugs/artificial insemination) during the last 12 months?  "used_livestock_tech" 
#8#B2	B2. Did an expert (e.g. crop or livestock extension agent, or community based facilitator or another experienced farmer) visit your home in the last 12 months? "b21"
#9# B3	B3. During the last 12 months, did you or someone in the household visit an extension office or a meeting/training organized by an extension officer? "b31"
#10# B3.4	B3.4  Are there any agricultural enterprises, improved technologies or inputs you would like to adopt?
#11# B3.5	B3.5  Are the extension agents/farmer forum members in the village/parish aware of this need?
#12# B3.20	B3.20 How are new agricultural enterprises that are being promoted by the government decided? (neg: True if decided without consultation)
#13#B4 B4. Are there any farmer associations/groups in this village?
#14# B4.1 B4.1 Are any of these farmer groups supported by Naads or Operation Wealth Creation?
#15# B5.2B5.2  Did you receive any help in marketing your produce from the Village procurement committe/Village farmers forum/Village farmers forum executive in the last 12 months? 
#16# B5.3 B5.3 Did you receive any help in marketing your produce from the Cooperative/Association in the last 12 months?

# Section C: water
#17# [8,] "unprotected"       "base_unprotected"      	Household uses unprotected water source during dry season (yes/no)
##18 ## [9,] "baraza.C1.2"       "c12source"             	Distance to water source
#19#[10,] "baraza.C1.3"       "qc15"                  	Average waiting time at source (min)
#20#[11,] "baraza.C2.3"       "c10"             		Is there a Water User Committee in this village? (yes/no)
#21# C1.4 How satisfied are you with the qualiy of water available at this water source during the dry season? == very satisfied or satisfied "qc16"
#22#C2.1 How do you treat water before drinking?  (boil or treat) "c4"
#23#C2.4  Are you or any member of the household a member of the Committee? 1 = yes, nas are 0 "c11"
#24#C2.5 Does the water committee hold public meetings? 1 = yes, nas are 0 "c13"

### Section D: health

#25#[14,] "baraza.D2"         "pub_health_access"     	Seek treatment for fever in public health facility (1=yes) 
#26#[15,] "baraza.D2.4"       "maternal_health_access"	Go to public health facility to give birth (1=yes)  
#27#[16,] "baraza.D3"         "d31"                   	Is there a VHT in village? (1=yes) 
#28# 31 D3.1 Are you or any member of the household part of the VHT? d32
#29# D3.3 Did the VHT organise any public meetings in your village in the last 12 months? d315
#30#[17,] "baraza.D4.2"       "d43"                   	Distance to nearest govt health facility (km) 
#31  	baraza.D1 Any members sick?
#32#[18,] "baraza.D1.2"	"tot_sick"		Number of days ill
#33"baraza.D1.3" number of days school/work missed due to illness
#34#[19,] "baraza.D4.6"       "wait_time"           	Waiting time before being attended (min)  
#35	"baraza.D6" 		"d61"			Has visited traditional health practitioner? (1=yes)  
#36 doctor was examined by in-charge/doctor 
#37 baraza.D4.7 time for examination d411"
#38# paid anything
#39	baraza.D4.11		received meds d416
#40 baraza.D4.12 had to buy meds d419
#41 baraza.D4.13	satisfied with services at hospital d420
#42"baraza.D4.14"MHU at govt facility d426

##43 "n_children"	    "base_n_children"		Number of children in UPS or USE 
##44"baraza.E1"		e5				Distance to public school (km)  
##45"baraza.E1.4" "e12" 				Has complete boundary fence (1=yes) 
##46"baraza.E1.5" "e13" 				Has electricity (1=yes) 
##47"baraza.E1.6" "e14" 				Has water facility (1=yes) 
##48 baraza.E18 Were any Parent Teacher Association (PTA) meetings held in that primary UPE school during the last 12 months? 
##49"baraza.E1.10 "e22" 				Has School Management Committee (1=yes)  
##50"baraza.E1.13 "e32" 				Is informed about School Management Committee (1=yes)  
##51"baraza.E1.18 "e45" 				Inspectors visited schools (1=yes)
##52 baraza.A6"   Distance to nearest all weather road (km)       "a6

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
endline <- FW_index(c("unprotected", "baraza.C1.2", "baraza.C1.3","baraza.C2.3","baraza.A6"),revcols=c(1,2,3,5),data=endline)
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
##17 number of days sick
###18 wait time
## 19 index
endline <- FW_index(c("baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1.2",  "baraza.D4.6", "baraza.D6"),revcols=c(4,5,6,7),data=endline)
names(endline)[names(endline) == 'index'] <- 'health_index'
baseline <- FW_index(c("pub_health_access","maternal_health_access","d31","d43","tot_sick","wait_time","d61"),revcols=c(4,5,6,7),data=baseline)
names(baseline)[names(baseline) == 'index'] <- 'base_health_index'
baseline_matching <- FW_index(c("pub_health_access","maternal_health_access","d31","d43","tot_sick","wait_time","d61"),revcols=c(4,5,6,7),data=baseline_matching)
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


outcomes <- c("baraza.B2","baraza.B3","baraza.B4.1","inputs","baraza.B5.2","baraza.B5.3","ag_index","unprotected", "baraza.C1.2", "baraza.C1.3","baraza.C2.3","baraza.A6","infra_index","baraza.D2","baraza.D2.4","baraza.D3","baraza.D4.2", "baraza.D1.2",  "baraza.D4.6","baraza.D6","health_index","n_children","baraza.E5","baraza.E12","baraza.E14","baraza.E22","baraza.E32","baraza.E45","education_index", "pub_service_index")
baseline_outcomes <- c("b21","b31","b44","base_inputs","b5144","b5146","base_ag_index","base_unprotected","c12source", "qc15","c10","a6","base_infra_index","pub_health_access","maternal_health_access","d31","d43","tot_sick","wait_time","d61","base_health_index","base_n_children","e5","e12", "e14","e22","e32","e45","base_education_index","base_pub_service_index")

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
df_matcher <- array(NA,dim=c(6,5,length(outcomes)))
df_averages <- array(NA,dim=c(2,length(outcomes)))

for (i in 1:length(outcomes)) {
print(i)
# i <- 1



###matched dif-in-dif

baseline_complete <- baseline_matching[complete.cases(baseline_matching[c("information", "deliberation","district_baraza","hhsize","femhead","agehead","log_farmsize","ironroof","improved_wall","has_phone","head_sec","a26a","a26b","hhid","clusterID","clusterID2",baseline_outcomes[i])]),] 
baseline_complete <- baseline_complete[c("information","deliberation","district_baraza","hhsize","femhead","agehead","log_farmsize","ironroof","improved_wall","has_phone","head_sec","a26a","a26b","hhid","clusterID","clusterID2",baseline_outcomes[i])]

####matching for information
a26a_cut <- seq(min(baseline_complete$a26a),max(baseline_complete$a26a),by=1)
a26b_cut <- seq(min(baseline_complete$a26b),max(baseline_complete$a26b),by=1)
agehead_cut <- seq(min(baseline_complete$agehead),max(baseline_complete$agehead),by=20)
log_farmsize_cut <- seq(min(baseline_complete$log_farmsize),max(baseline_complete$log_farmsize),by=1)
log_farmsize_cut <- seq(min(baseline_complete$log_farmsize),max(baseline_complete$log_farmsize),by=1)
#log_farmsize_cut <- c(0,1.7,6)
my.cutpoints <- list(a26a=a26a_cut, a26b=a26b_cut, agehead = agehead_cut, log_farmsize=log_farmsize_cut)

nearest.match <- matchit(formula = information ~ hhsize + femhead + agehead  + head_sec+ log_farmsize + ironroof + improved_wall+has_phone+a26a+a26b,  data =baseline_complete[baseline_complete$district_baraza == 0,] ,method = "cem", cutpoints = my.cutpoints )

summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)

### this is the matched data -  we now need to extract these ids from the endline and stack base and endline to do a dif-in-dif

matched.merged <-   merge(matched.baseline, endline[c("hhid","a21","district","subcounty",outcomes[i])], by="hhid")

df_averages[1,i] <- mean(as.matrix(matched.merged[outcomes[i]]), na.rm=T)
df_averages[2,i] <- sd(as.matrix(matched.merged[outcomes[i]]), na.rm=T)

ols <- lm(as.formula(paste(paste(outcomes[i],"information*deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=matched.merged, weights= weights) 
vcov_cluster <- vcovCR(ols, cluster = matched.merged$clusterID, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes, baseline_outcomes, matched.merged , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <- RI_store$conf_2 
res[2,5] <- RI_store$pval_2
}

df_matcher[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4], conf[2,5], nobs(ols))

####matching for deliberation
nearest.match <- matchit(formula = deliberation ~  hhsize + femhead + agehead  + head_sec+ log_farmsize + ironroof + improved_wall+has_phone+a26a+a26b,   data =baseline_complete[baseline_complete$district_baraza == 0,]  ,method = "cem", cutpoints = my.cutpoints )
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)

matched.merged <-   merge(matched.baseline, endline[c("hhid","a21","district","subcounty",outcomes[i])], by="hhid")


ols <- lm(as.formula(paste(paste(outcomes[i],"information*deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=matched.merged, weights= weights) 
vcov_cluster <- vcovCR(ols, cluster = matched.merged$clusterID, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes, baseline_outcomes, matched.merged , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[3,4:5] <- RI_store$conf_3 
res[3,5] <- RI_store$pval_3
}

### this is the matched data -  we now need to extract these ids from the endline and stack base and endline to do a dif-in-dif
df_matcher[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4], conf[3,5], nobs(ols))


####matching for interaction
nearest.match <- matchit(formula = deliberation*information ~ hhsize + femhead + agehead  + head_sec+ log_farmsize + ironroof + improved_wall+has_phone+a26a+a26b,  data =baseline_complete[baseline_complete$district_baraza == 0 & (baseline_complete$information == baseline_complete$deliberation),]  ,method = "cem", cutpoints = my.cutpoints )
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)
matched.merged <-   merge(matched.baseline, endline[c("hhid","a21","district","subcounty",outcomes[i])], by="hhid")

####
ols <- lm(as.formula(paste(paste(outcomes[i],"information:deliberation+a21",sep="~"),baseline_outcomes[i],sep="+")), data=matched.merged, weights= weights) 
vcov_cluster <- vcovCR(ols, cluster = matched.merged$clusterID, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc_custom(i,outcomes, baseline_outcomes, matched.merged , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[6,4:5] <- RI_store$conf_1
res[6,5] <- RI_store$pval_1
}

df_matcher[,1,i] <- c(res[6,1],res[6,2],res[6,5], conf[6,4],conf[6,5], nobs(ols))


####matching for district baraza

nearest.match <- matchit(formula = district_baraza ~hhsize + femhead + agehead  + head_sec+ log_farmsize + ironroof + improved_wall+has_phone,  data =baseline_complete[(baseline_complete$information == 0 & baseline_complete$deliberation==0) | baseline_complete$district_baraza == 1 ,]  ,method = "cem", cutpoints = my.cutpoints )
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)
matched.merged <-   merge(matched.baseline, endline[c("hhid","a21","district","subcounty",outcomes[i])], by="hhid")


ols <- lm(as.formula(paste(paste(outcomes[i],"district_baraza+a21",sep="~"),baseline_outcomes[i],sep="+")), data=matched.merged, weights= weights ) 
vcov_cluster <- vcovCR(ols, cluster = matched.merged$clusterID2, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_dist(i,outcomes, baseline_outcomes, matched.merged , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <-  RI_store$conf
res[2,5] <- RI_store$pval
}

df_matcher[,4,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4], conf[2,5], nobs(ols))

####matching for district baraza (dif with sc baraza)

nearest.match <- matchit(formula = district_baraza ~hhsize + femhead + agehead  + head_sec+ log_farmsize + ironroof + improved_wall+has_phone,  data =baseline_complete[(baseline_complete$information == 1 & baseline_complete$deliberation==1) | baseline_complete$district_baraza == 1 ,]  ,method = "cem", cutpoints = my.cutpoints )
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)
matched.merged <-   merge(matched.baseline, endline[c("hhid","a21","district","subcounty",outcomes[i])], by="hhid")


ols <- lm(as.formula(paste(paste(outcomes[i],"district_baraza+a21",sep="~"),baseline_outcomes[i],sep="+")), data=matched.merged, weights= weights ) 
vcov_cluster <- vcovCR(ols, cluster = matched.merged$clusterID2, type = "CR0")
res <- coef_test(ols, vcov_cluster)
conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_dist(i,outcomes, baseline_outcomes, matched.merged , ctrls = "a21", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <-  RI_store$conf
res[2,5] <- RI_store$pval
}

df_matcher[,5,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4], conf[2,5], nobs(ols))

}

### create data.frame to plot - make sure you get correct i's for the indices; last one is overall index
d_plot <- data.frame(rbind(df_matcher[c(1,4,5),1,7],df_matcher[c(1,4,5),2,7],df_matcher[c(1,4,5),3,7],df_matcher[c(1,4,5),5,7]))
d_plot <- rbind(d_plot,data.frame(rbind(df_matcher[c(1,4,5),1,13],df_matcher[c(1,4,5),2,13],df_matcher[c(1,4,5),3,13],df_matcher[c(1,4,5),5,13])))
d_plot <- rbind(d_plot,data.frame(rbind(df_matcher[c(1,4,5),1,21],df_matcher[c(1,4,5),2,21],df_matcher[c(1,4,5),3,21],df_matcher[c(1,4,5),5,21])))
d_plot <- rbind(d_plot,data.frame(rbind(df_matcher[c(1,4,5),1,29],df_matcher[c(1,4,5),2,29],df_matcher[c(1,4,5),3,29],df_matcher[c(1,4,5),5,29])))
d_plot <- rbind(d_plot, data.frame(rbind(c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA),c(NA,NA,NA))))
d_plot <- rbind(d_plot,data.frame(rbind(df_matcher[c(1,4,5),1,30],df_matcher[c(1,4,5),2,30],df_matcher[c(1,4,5),3,30],df_matcher[c(1,4,5),5,30])))


names(d_plot) <- c("y","ylo","yhi")

d_plot$x <- rep(c("agricuture","infrastructure","health","education","","index"), each=4)
d_plot$grp <- rep(c("sc baraza","info","delib","level"), times=6)
d_plot$grp <-  factor(d_plot$grp , levels=c("sc baraza","info","delib","level"))
d_plot$x <-  factor(d_plot$x, levels=rev((c("agricuture","infrastructure","health","education","","index"))))


### save results
save_path <- ifelse(final_verion_swith, paste(path,"report/results/final", sep = "/"), paste(path,"report/results/", sep = "/"))


save(df_matcher, file= paste(save_path,"df_matcher.Rd", sep="/"))
save(df_averages, file= paste(save_path,"df_averages_matcher.Rd", sep="/"))


png(paste(save_path,"impact_summary_matcher.png",sep = "/"), units="px", height=3200, width= 6400, res=600)
print(credplot.gg(d_plot,'SDs','',levels(d_plot$x),.3))
dev.off()





