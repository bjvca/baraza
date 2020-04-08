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

#look at what is typed if B1 == "other"
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

sc_merged <- merge(sc_baseline,sc_endline,by.x = c("district","subcounty","actor"), by.y=c("district","subcounty","actor"))

########RECODING########
#SECTION E: GOVERNMENT AND COMMUNITY BODIES IN THE SUBCOUNTY#
sc_merged$baraza.production.E1a_binary <- (sc_merged$baraza.production.E1a == 1)
sc_merged$baraza.health.E2a_binary <- (sc_merged$baraza.health.E2a == 1)
sc_merged$baraza.gender1.E3a_binary <- (sc_merged$baraza.gender1.E3a == 1)
sc_merged$baraza.works.E4a_binary <- (sc_merged$baraza.works.E4a == 1)
sc_merged$baraza.finance1.E5a_binary <- (sc_merged$baraza.finance1.E5a == 1)
sc_merged$baraza.E7[sc_merged$baraza.E7==96] <- NA
sc_merged$baraza.E7_binary <- (sc_merged$baraza.E7 > 2) #median is 2

sc_merged$e11a_binary <- (sc_merged$e11a == "Yes")
sc_merged$e11b_binary <- (sc_merged$e11b == "Yes")
sc_merged$e11c_binary <- (sc_merged$e11c == "Yes")
sc_merged$e11d_binary <- (sc_merged$e11d == "Yes")
sc_merged$e11e_binary <- (sc_merged$e11e == "Yes")
sc_merged$e13_binary <- (sc_merged$e13 %in% c("After every two months", "Monthly","Quarterly"))

#SECTION F: COMMUNITY MEETINGS HELD IN THE SUBCOUNTY#
sc_merged$f14a[sc_merged$f14a==2015] <- NA

#SUBSECTION H1 - HEALTH#
sc_merged$baraza.H3 <- as.numeric(as.character(sc_merged$baraza.H3))
sc_merged$baraza.H4 <- as.numeric(as.character(sc_merged$baraza.H4))
sc_merged$baraza.H5 <- as.numeric(as.character(sc_merged$baraza.H5))
sc_merged$baraza.H6 <- as.numeric(as.character(sc_merged$baraza.H6))
sc_merged$baraza.H7 <- as.numeric(as.character(sc_merged$baraza.H7))
sc_merged$baraza.H8 <- as.numeric(as.character(sc_merged$baraza.H8))
sc_merged$baraza.H9 <- as.numeric(as.character(sc_merged$baraza.H9))
sc_merged$baraza.H10 <- as.numeric(as.character(sc_merged$baraza.H10))
sc_merged$baraza.H11 <- as.numeric(as.character(sc_merged$baraza.H11))
sc_merged$baraza.H12 <- as.numeric(as.character(sc_merged$baraza.H12))
sc_merged$baraza.H13 <- as.numeric(as.character(sc_merged$baraza.H13))
sc_merged$baraza.H15 <- as.numeric(as.character(sc_merged$baraza.H15))
sc_merged$baraza.H16 <- as.numeric(as.character(sc_merged$baraza.H16))
sc_merged$baraza.H17 <- as.numeric(as.character(sc_merged$baraza.H17))
sc_merged$baraza.H18 <- as.numeric(as.character(sc_merged$baraza.H18))
sc_merged$baraza.H19 <- as.numeric(as.character(sc_merged$baraza.H19))
sc_merged$baraza.H20 <- as.numeric(as.character(sc_merged$baraza.H20))
sc_merged$baraza.H21 <- as.numeric(as.character(sc_merged$baraza.H21))
sc_merged$baraza.H22 <- as.numeric(as.character(sc_merged$baraza.H22))
sc_merged$baraza.H23 <- as.numeric(as.character(sc_merged$baraza.H23))
sc_merged$baraza.H24 <- as.numeric(as.character(sc_merged$baraza.H24))
sc_merged$baraza.H25 <- as.numeric(as.character(sc_merged$baraza.H25))
sc_merged$baraza.H26 <- as.numeric(as.character(sc_merged$baraza.H26))
sc_merged$baraza.H27 <- as.numeric(as.character(sc_merged$baraza.H27))
sc_merged$baraza.H28 <- as.numeric(as.character(sc_merged$baraza.H28))
sc_merged$baraza.H29 <- as.numeric(as.character(sc_merged$baraza.H29))
sc_merged$baraza.H30 <- as.numeric(as.character(sc_merged$baraza.H30))
sc_merged$baraza.H31 <- as.numeric(as.character(sc_merged$baraza.H31))
sc_merged$baraza.H32 <- as.numeric(as.character(sc_merged$baraza.H32))

sc_merged$d15a[sc_merged$d15a==0.60000002] <- 60

sc_merged$sum_h182a_h183a <- (sc_merged$h182a + sc_merged$h183a)
sc_merged$sum_h182b_h183b <- (sc_merged$h182b + sc_merged$h183b)
sc_merged$sum_h182d_h183d <- (sc_merged$h182d + sc_merged$h183d)

sc_merged$sum_h1123a_h1124a <- (sc_merged$h1123a + sc_merged$h1124a)
sc_merged$sum_h1123b_h1124b <- (sc_merged$h1123b + sc_merged$h1124b)
sc_merged$sum_h1123c_h1124c <- (sc_merged$h1123c + sc_merged$h1124c)

sc_merged <- sc_merged %>%  mutate(clusterID = group_indices(., district, subcounty))

########LOOPS########
outcomes <- c("baraza.km.D1","baraza.km.D2","baraza.km.D3","baraza.km.D4a","baraza.km.D4b","baraza.production.E1a_binary","baraza.health.E2a_binary","baraza.gender1.E3a_binary","baraza.works.E4a_binary","baraza.finance1.E5a_binary","baraza.E7_binary","baraza.meeting.F1","baraza.meeting.F2","baraza.meeting.F3","baraza.meeting.F4")
baseline_outcomes <- c("d12","d13","d14","d15a","d15b","e11a_binary","e11b_binary","e11c_binary","e11d_binary","e11e_binary","e13_binary","f14a","f14b","f14c","f14d")
#df_ols <- array(NA,dim=c(6,3,length(outcomes)))
#df_ols <- array(NA,dim=c(7,3,length(outcomes)))
df_ols <- array(NA,dim=c(3,3,length(outcomes)))

sc_merged[outcomes] <- lapply(sc_merged[outcomes], function(x) replace(x, x == 999, NA) )
sc_merged[outcomes] <- lapply(sc_merged[outcomes], function(x) replace(x, x == "n/a", NA) )

for (i in 1:length(outcomes)) {
  print(i)
  ols <- lm(as.formula(paste(paste(outcomes[i],"information*deliberation+region.x",sep="~"),baseline_outcomes[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,5],res[2,6],conf[2,4],conf[2,5],nobs(ols))
  #df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,5],conf[2,4],conf[2,5],nobs(ols))
  df_ols[,2,i] <- c(res[2,1],res[2,5],nobs(ols))
  #df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,5],res[3,6],conf[3,4],conf[3,5],nobs(ols))
  #df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,5],conf[3,4],conf[3,5],nobs(ols))
  df_ols[,3,i] <- c(res[3,1],res[3,5],nobs(ols))
  ols <- lm(as.formula(paste(paste(outcomes[i],"information:deliberation+region.x",sep="~"),baseline_outcomes[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols[,1,i] <- c(res[6,1],res[6,2],res[6,5],res[6,6],conf[6,4],conf[6,5],nobs(ols))
  #df_ols[,1,i] <- c(res[6,1],res[6,2],res[6,5],conf[6,4],conf[6,5],nobs(ols))
  df_ols[,1,i] <- c(res[6,1],res[6,5],nobs(ols))
}

#loop for health variables
sc_merged$h181a <- replace(sc_merged$h181a, is.na(sc_merged$h181a), 0)
sc_merged$h181b <- replace(sc_merged$h181b, is.na(sc_merged$h181b), 0)
sc_merged$h181d <- replace(sc_merged$h181d, is.na(sc_merged$h181d), 0)

outcomes_health <- c("baraza.H1","baraza.H2","baraza.H2b","baraza.H3","baraza.H4","baraza.H5","baraza.H6","baraza.H7","baraza.H8","baraza.H9","baraza.H10","baraza.H11","baraza.H12","baraza.H13","baraza.H18","baraza.H19","baraza.H20","baraza.H21","baraza.H22","baraza.H23","baraza.H24","baraza.H25","baraza.H26","baraza.H27","baraza.H28","baraza.H29","baraza.H30","baraza.H31","baraza.H32")
baseline_outcomes_health <- c("h11","h12","h15","h181a","h181b","h181d","sum_h182a_h183a","sum_h182b_h183b","sum_h182d_h183d","h184a","h184b","h184d","h110","h111","h1121a","h1121b","h1121c","h1122a","h1122b","h1122c","sum_h1123a_h1124a","sum_h1123b_h1124b","sum_h1123c_h1124c","h1125a","h1125b","h1125c","h1126a","h1126b","h1126c")
#df_ols_health <- array(NA,dim=c(6,10,length(outcomes_health)))
df_ols_health <- array(NA,dim=c(3,9,length(outcomes_health)))
sc_merged[outcomes_health] <- lapply(sc_merged[outcomes_health], function(x) replace(x, x == 999, NA) )

for (i in 1:length(outcomes_health)) {
  #A: NAs remain NAs
  print(i)
  ols <- lm(as.formula(paste(paste(outcomes_health[i],"information*deliberation+region.x",sep="~"),baseline_outcomes_health[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols_health[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5],nobs(ols))
  df_ols_health[,2,i] <- c(res[2,1],res[2,5],nobs(ols))
  #df_ols_health[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5],nobs(ols))
  df_ols_health[,3,i] <- c(res[3,1],res[3,5],nobs(ols))
  
  ols <- lm(as.formula(paste(paste(outcomes_health[i],"information:deliberation+region.x",sep="~"),baseline_outcomes_health[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols_health[,1,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5],nobs(ols))
  df_ols_health[,1,i] <- c(res[6,1],res[6,5],nobs(ols))
  
  #B: recode NAs in baseline data as 0
  sc_merged[baseline_outcomes_health[i]] <- replace(sc_merged[baseline_outcomes_health[i]], is.na(sc_merged[baseline_outcomes_health[i]]), 0)
  
  print(i)
  ols <- lm(as.formula(paste(paste(outcomes_health[i],"information*deliberation+region.x",sep="~"),baseline_outcomes_health[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols_health[,5,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5],nobs(ols))
  df_ols_health[,5,i] <- c(res[2,1],res[2,5],nobs(ols))
  #df_ols_health[,6,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5],nobs(ols))
  df_ols_health[,6,i] <- c(res[3,1],res[3,5],nobs(ols))
  ols <- lm(as.formula(paste(paste(outcomes_health[i],"information:deliberation+region.x",sep="~"),baseline_outcomes_health[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols_health[,4,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5],nobs(ols))
  df_ols_health[,4,i] <- c(res[6,1],res[6,5],nobs(ols))
    
  #C: recode NAs in baseline and endline data as 0
  sc_merged[outcomes_health[i]] <- replace(sc_merged[outcomes_health[i]], is.na(sc_merged[outcomes_health[i]]), 0)
  
  print(i)
  ols <- lm(as.formula(paste(paste(outcomes_health[i],"information*deliberation+region.x",sep="~"),baseline_outcomes_health[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols_health[,8,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5],nobs(ols))
  df_ols_health[,8,i] <- c(res[2,1],res[2,5],nobs(ols))
  #df_ols_health[,9,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5],nobs(ols))
  df_ols_health[,9,i] <- c(res[3,1],res[3,5],nobs(ols))
  
  ols <- lm(as.formula(paste(paste(outcomes_health[i],"information:deliberation+region.x",sep="~"),baseline_outcomes_health[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols_health[,7,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5],nobs(ols))
  df_ols_health[,7,i] <- c(res[6,1],res[6,5],nobs(ols))
  
}

#waiting for Bjorns answer
#outcomes_H1h <- c("baraza.H15","baraza.H16","baraza.H17")
#baseline_outcomes_H1h <- c("?","?","?")
