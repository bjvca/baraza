rm(list=ls())

library(dplyr)
library(ggplot2)
library(MatchIt) 
library(multiwayvcov)
library(plm)
library(lmtest)
library(clubSandwich)
library(moments)

if (Sys.info()['sysname'] =="Windows") {
path <- "C:/users/u0127963/Desktop/PhD/baraza"
} else {
path <- "/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline"
}

#load sc endline and baseline data
sc_endline <- read.csv(paste(path,"data/public/sc_level_endline.csv", sep ="/"))
sc_baseline <- read.csv(paste(path,"data/public/sc_level_baseline.csv", sep ="/"))
sc_baseline$actor <- NA
sc_baseline$actor[sc_baseline$designation %in% c("LC3_chair","Vice Chairperson LCIII")] <- "politician"
sc_baseline$actor[sc_baseline$designation %in% c("Subcounty Chief/Town Clerk","Parish Chief","Acting Subcounty Chief/Town Clerk","Deputy subcounty chief/Town Clerk","Community Development Officer")] <- "civil servant"
sc_endline$actor <- NA
sc_endline$actor[sc_endline$baraza.B1 %in% c(1,3,5,6) ]  <- "civil servant"
sc_endline$actor[sc_endline$baraza.B1 %in% c(2,4) ]  <- "politician"

#look at what is typed in designation if B1 == "other"
#baraza.other_b1
sc_endline[sc_endline$baraza.B1==96,c("baraza.other_b1","scID")]

sc_endline$actor[sc_endline$scID == 9] <- "civil servant" #                         Secretary for Production    9
sc_endline$actor[sc_endline$scID == 28] <- "civil servant" #        Secretary for General Purpose and Works   28
sc_endline$actor[sc_endline$scID == 44] <- "politician" #                       Vice Chairperson LC 3   44
sc_endline$actor[sc_endline$scID == 58] <- "politician" #                        VICE LC III CHAIRPERSON   58
sc_endline$actor[sc_endline$scID == 61] <- "politician" #                          Vice LC3 chairperson   61
sc_endline$actor[sc_endline$scID == 66] <- "politician" #                         Vice Chairperson LC 3   66
sc_endline$actor[sc_endline$scID == 68] <- "politician" #                                  Vice chairperson   68
sc_endline$actor[sc_endline$scID == 74] <- "civil servant" #                     Parish chief   74
sc_endline$actor[sc_endline$scID == 76] <- "politician" #                             Lc 3 VICE CHAIRPERSON   76
sc_endline$actor[sc_endline$scID == 79] <- "politician" #                               Vice chairman L c 3   79
sc_endline$actor[sc_endline$scID == 91] <- "politician" #                Secretary to the chairman LC III   91
sc_endline$actor[sc_endline$scID == 101] <- "politician" #                            Vice Chairperson LC3  101
sc_endline$actor[sc_endline$scID == 107] <- "civil servant" #                         Town board treasurer  107
sc_endline$actor[sc_endline$scID == 109] <- "politician" #                               Vice Chairperson   109
sc_endline$actor[sc_endline$scID == 110] <- "civil servant" #                                   Town agent  110
sc_endline$actor[sc_endline$scID == 112] <- "politician" #                        Concillor Bunasaka parish  112
sc_endline$actor[sc_endline$scID == 116] <- "politician" #               Head of finance (Area councillor).  116
sc_endline$actor[sc_endline$scID == 120]  <- "civil servant" # OFFICER FOR AGRICULTURE AND FISHERIES DEPARTMENT  120
sc_endline$actor[sc_endline$scID == 122] <- "politician" #                            Sub County Councilor  122
sc_endline$actor[sc_endline$scID == 134] <- "politician" #                                      LC3 VICE    134
sc_endline$actor[sc_endline$scID == 150] <- "politician" #                       Secretary Social services.  150
sc_endline$actor[sc_endline$scID == 187] <- "politician" #                          VICE CHAIRPERSON LCIII  187
sc_endline$actor[sc_endline$scID == 197] <- "politician" #    Interim L.C 3  Chairperson  Bumbo T. Council  197
sc_endline$actor[sc_endline$scID == 199]  <- "civil servant"  #                              Health Assistant  199
sc_endline$actor[sc_endline$scID == 200] <- "politician" #                          Vice  Chairperson  LC3  200
sc_endline$actor[sc_endline$scID == 215] <- "politician" #                                       Councillor  215
sc_endline$actor[sc_endline$scID == 221] <- "politician" #                                      Councillor   221
sc_endline$actor[sc_endline$scID == 225] <- "politician" #                                   LC3 Councillor  225
sc_endline$actor[sc_endline$scID == 232] <- "politician" #                                Youth   Councilor  232
sc_endline$actor[sc_endline$scID == 234] <- "politician" #                     FINANCIAL COMMITTEE CHAIRMAN  234
sc_endline$actor[sc_endline$scID == 249] <- "politician" #                    Councilor  Bufumbo  Subcounty  249
sc_endline$actor[sc_endline$scID == 252] <- "politician" #                            S/c vice Chairperson  252
sc_endline$actor[sc_endline$scID == 254] <- "politician" #                                     Vice Chair  254
sc_endline$actor[sc_endline$scID == 256] <- "politician" #                                    Councilor PWD  256
sc_endline$actor[sc_endline$scID == 258] <- "civil servant" #                                    Town clerk  258

#load treatment assignment
treats <- read.csv(paste(path,"questionnaire/final_list_5.csv", sep ="/"))
#merge in treatments
sc_endline <- merge(treats, sc_endline, by.x=c("district","subcounty"), by.y=c("district","sub"))

#setdiff(sc_endline$subcounty,sc_baseline$subcounty)
#[1] "HARUGALI"     "KIGULYA": these two were apparenty not inteviewed during baseline? 
#    "NTUUSI"       "SEMBABULE_TC":: these two have spelling errors

sc_baseline$subcounty <- as.character(sc_baseline$subcounty)
sc_baseline$subcounty[sc_baseline$subcounty == "SEMBABULE TC"] <- "SEMBABULE_TC"
sc_baseline$subcounty[sc_baseline$subcounty ==  "NTUSI"] <- "NTUUSI"

sc_merged <- merge(sc_baseline,sc_endline,by.x = c("subcounty","actor"), by.y=c("subcounty","actor"))

#drop sc with district treatment
sc_merged_no_dist <- subset(sc_merged, district_baraza == 0)

#make indicator that pools I, D and ID
sc_merged_no_dist$ind_treat <- (sc_merged_no_dist$information == 1 | sc_merged_no_dist$deliberation == 1)

########LOOPS########

#SECTION D: SUBCOUNTY'S BASIC INFORMATION#
outcomes_D <- c("baraza.km.D1","baraza.km.D2","baraza.km.D3","baraza.km.D4a","baraza.km.D4b")
baseline_outcomes_D <- c("d12","d13","d14","d15a","d15b")

df_ols <- array(NA,dim=c(6,4,length(outcomes_D)))

sc_merged[outcomes_D] <- lapply(sc_merged[outcomes_D], function(x) replace(x, x == 999, NA) )

for (i in 1:length(outcomes_D)) {
  print(i)
  
  ols <- lm(as.formula(paste(paste(outcomes_D[i],"information*deliberation+region.x",sep="~"),baseline_outcomes_D[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$subcounty[sc_merged$district_baraza == 0], type = "CR0")
  
  res <- coef_test(ols, vcov_cluster)
  
  #conf <- conf_int(ols, vcov_cluster)
  
  print(res)
  #print(conf)
  
  #df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5],
                    #nobs(ols))
  #df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5],
                    #nobs(ols))
  
  ols <- lm(as.formula(paste(paste(outcomes_D[i],"information:deliberation+region.x",sep="~"),baseline_outcomes_D[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$subcounty[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  #conf <- conf_int(ols, vcov_cluster)
  
  print(res)
  #print(conf)
  
  #df_ols[,1,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5],
  #nobs(ols))
  
}

#SECTION E: GOVERNMENT AND COMMUNITY BODIES IN THE SUBCOUNTY#
outcomes_E <- c("baraza.production.E1a","baraza.health.E2a","baraza.gender1.E3a","baraza.works.E4a","baraza.finance1.E5a","baraza.E7")
baseline_outcomes_E <- c("e11a","e11b","e11c","e11d","e11e","e13")

df_ols <- array(NA,dim=c(6,4,length(outcomes_E)))

sc_merged[outcomes_E] <- lapply(sc_merged[outcomes_E], function(x) replace(x, x == 999, NA) )

sc_merged$baraza.production.E1a_binary <- (sc_merged$baraza.production.E1a == 1)
sc_merged$baraza.health.E2a_binary <- (sc_merged$baraza.health.E2a == 1)
sc_merged$baraza.gender1.E3a_binary <- (sc_merged$baraza.gender1.E3a == 1)
sc_merged$baraza.works.E4a_binary <- (sc_merged$baraza.works.E4a == 1)
sc_merged$baraza.finance1.E5a_binary <- (sc_merged$baraza.finance1.E5a == 1)
sc_merged$baraza.E7[sc_merged$baraza.E7==96] <- NA
#median is 2
sc_merged$baraza.E7_binary <- (sc_merged$baraza.E7 > 2)

sc_merged$e11a_binary <- (sc_merged$e11a == "Yes")
sc_merged$e11b_binary <- (sc_merged$e11b == "Yes")
sc_merged$e11c_binary <- (sc_merged$e11c == "Yes")
sc_merged$e11d_binary <- (sc_merged$e11d == "Yes")
sc_merged$e11e_binary <- (sc_merged$e11e == "Yes")
sc_merged$e13_binary <- (sc_merged$e13 %in% c("After every two months", "Monthly","Quarterly"))

outcomes_E_binary <- c("baraza.production.E1a_binary","baraza.health.E2a_binary","baraza.gender1.E3a_binary","baraza.works.E4a_binary","baraza.finance1.E5a_binary","baraza.E7_binary")
baseline_outcomes_E_binary <- c("e11a_binary","e11b_binary","e11c_binary","e11d_binary","e11e_binary","e13_binary")


for (i in 1:length(outcomes_E_binary)) {
  print(i)
  
  ols <- lm(as.formula(paste(paste(outcomes_E_binary[i],"information*deliberation+region.x",sep="~"),baseline_outcomes_E_binary[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$subcounty[sc_merged$district_baraza == 0], type = "CR0")
  
  res <- coef_test(ols, vcov_cluster)
  
  #conf <- conf_int(ols, vcov_cluster)
  
  print(res)
  #print(conf)
  
  #df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5],
  #nobs(ols))
  #df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5],
  #nobs(ols))
  
  ols <- lm(as.formula(paste(paste(outcomes_E_binary[i],"information:deliberation+region.x",sep="~"),baseline_outcomes_E_binary[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$subcounty[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  #conf <- conf_int(ols, vcov_cluster)
  
  print(res)
  #print(conf)
  
  #df_ols[,1,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5],
  #nobs(ols))
  
}

#SECTION F: COMMUNITY MEETINGS HELD IN THE SUBCOUNTY#
sc_merged$f14a[sc_merged$f14a==2015] <- NA

outcomes_F <- c("baraza.meeting.F1","baraza.meeting.F2","baraza.meeting.F3","baraza.meeting.F4")
baseline_outcomes_F <- c("f14a","f14b","f14c","f14d")

df_ols <- array(NA,dim=c(6,4,length(outcomes_F)))

sc_merged[outcomes_F] <- lapply(sc_merged[outcomes_F], function(x) replace(x, x == 999, NA) )

for (i in 1:length(outcomes_F)) {
  print(i)
  
  ols <- lm(as.formula(paste(paste(outcomes_F[i],"information*deliberation+region.x",sep="~"),baseline_outcomes_F[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$subcounty[sc_merged$district_baraza == 0], type = "CR0")
  
  res <- coef_test(ols, vcov_cluster)
  
  #conf <- conf_int(ols, vcov_cluster)
  
  print(res)
  #print(conf)
  
  #df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5],
  #nobs(ols))
  #df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5],
  #nobs(ols))
  
  ols <- lm(as.formula(paste(paste(outcomes_F[i],"information:deliberation+region.x",sep="~"),baseline_outcomes_F[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$subcounty[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  #conf <- conf_int(ols, vcov_cluster)
  
  print(res)
  #print(conf)
  
  #df_ols[,1,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5],
  #nobs(ols))
  
}