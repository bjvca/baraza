rm(list=ls())
library(dplyr)
library(ggplot2)
library(MatchIt) 
library(multiwayvcov)
library(plm)
library(lmtest)
library(clubSandwich)
library(moments)
library(doParallel)

if (Sys.info()['sysname'] =="Windows") {
path <- "C:/users/u0127963/Desktop/PhD/baraza"
} else {
path <- "/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline"
}

RI_conf_switch <- FALSE
glob_repli <- 1000
glob_sig <- c(.025,.975) ### 5 percent conf intervals

########################################################## function definitions #########################################################################################
RI_conf_sc <- function(i,outcomes, baseline_outcomes, dta_sim , ctrls = NULL, nr_repl = 1000, sig = c(.025,.975)) {
### a function to esimate confidence intervals using randomization inference following Gerber and Green pages 66-71 and 83.
#RI_conf_dist(2,outcomes, baseline_outcomes, subset(dta, ((information == 1 & deliberation==1) | district_baraza == 1)) , ctrls = "a21", nr_repl = 1000, sig = c(.025,.975))
#dta_sim <- dta[dta$district_baraza == 0 ,]
#ctrls <- "a21"
#nr_repl <- 1000
#sig <-  c(.025,.975)

	if (is.null(baseline_outcomes)) {
		formula1 <- as.formula(paste(outcomes[i],paste("information:deliberation",ctrls,sep="+"),sep="~"))
		formula2 <- as.formula(paste(outcomes[i],paste("information*deliberation",ctrls,sep="+"),sep="~"))
	} else {
		formula1 <- as.formula(paste(paste(outcomes[i],paste("information:deliberation",ctrls,sep="+"),sep="~"),baseline_outcomes[i],sep="+"))
		formula2 <- as.formula(paste(paste(outcomes[i],paste("information*deliberation",ctrls,sep="+"),sep="~"),baseline_outcomes[i],sep="+"))
	}


	dta_sim <- dta_sim %>%  mutate(clusterID = group_indices(., district, subcounty))
	### get ATEs for two different models
	ols_1 <- lm(formula1, data=dta_sim[dta_sim$deliberation == dta_sim$information,]) 
	ols_2 <- lm(formula2, data=dta_sim) 
	dta_sim$dep <- as.numeric(unlist(dta_sim[as.character(formula1[[2]])]))

	treat_nrs <- table(data.frame(aggregate(dta_sim[c("information","deliberation")], list(dta_sim$clusterID),mean))[,2:3])

	### calculate potential outcomes
	### for model 1 (sc effect)
	dta_sim$pot_out_0_1 <- NA
	dta_sim$pot_out_0_1[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0]
	dta_sim$pot_out_0_1[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] - coef(ols_1)["information:deliberation"]

	dta_sim$pot_out_1_1 <- NA
	dta_sim$pot_out_1_1[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1] <- dta_sim$dep[dta_sim["information"] == 1 & dta_sim["deliberation"] == 1]
	dta_sim$pot_out_1_1[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] <- dta_sim$dep[dta_sim["information"] == 0 & dta_sim["deliberation"] == 0] + coef(ols_1)["information:deliberation"]

	### for model 2 (factorial design with )
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

	oper <- foreach (repl = 1:nr_repl,.combine=rbind) %dopar% {
		#do the permuations	
		perm_treat <- data.frame(cbind(sample(c(rep("C",treat_nrs[1,1]), rep("D",treat_nrs[1,2]), rep("I",treat_nrs[2,1]),rep("B",treat_nrs[2,2]))),names(table(dta_sim$clusterID))))  
		names(perm_treat) <- c("perm_treat","clusterID")
		dta_perm <- merge(dta_sim, perm_treat, by.x="clusterID", by.y="clusterID")
		dta_perm$information <- ifelse(dta_perm$perm_treat %in% c("I","B"), 1, 0)
		dta_perm$deliberation <- ifelse(dta_perm$perm_treat %in% c("D","B"), 1, 0)

		dta_perm$dep_1 <- NA
		dta_perm$dep_1[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1] <- dta_perm$pot_out_1_1[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1]
		dta_perm$dep_1[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0] <- dta_perm$pot_out_0_1[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0]

		dta_perm$dep_2 <- NA
		dta_perm$dep_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1] <- dta_perm$pot_out_11_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==1]
		dta_perm$dep_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0] <- dta_perm$pot_out_00_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==0] 
		dta_perm$dep_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==0] <- dta_perm$pot_out_10_2[dta_perm["information"] ==1 & dta_perm["deliberation"] ==0] 
		dta_perm$dep_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==1] <- dta_perm$pot_out_01_2[dta_perm["information"] ==0 & dta_perm["deliberation"] ==1] 

### p-value
		exceed1 <- coef(lm(formula1, data=dta_perm))["information:deliberation"] > abs(coef(ols_1)["information:deliberation"])
		exceed2 <- coef(lm(formula2, data=dta_perm))["information"] > abs(coef(ols_2)["information"])
		exceed3 <- coef(lm(formula2, data=dta_perm))["deliberation"] > abs(coef(ols_2)["deliberation"])


		dta_perm[outcomes[i]] <- dta_perm$dep_1
		r1 <-coef(lm(formula1, data=dta_perm[dta_perm$deliberation == dta_perm$information,]))["information:deliberation"]

		dta_perm[outcomes[i]] <- dta_perm$dep_2
		r2 <-coef(lm(formula2, data=dta_perm))["information"]
		r3 <- coef(lm(formula2, data=dta_perm))["deliberation"]
		oper <- return(c(r1,r2,r3, exceed1, exceed2, exceed3))
	}
	return(list(conf_1 = quantile(oper[,1],sig),conf_2 = quantile(oper[,2],sig),conf_3 = quantile(oper[,3],sig), pval_1= (sum(oper[,4])/nr_repl)*2, pval_2= (sum(oper[,5])/nr_repl)*2, pval_3= (sum(oper[,6])/nr_repl)*2))
	}

########################################################################################### end functions definitions ###############################################################

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

sc_merged <- sc_merged %>%  mutate(clusterID = group_indices(., district, subcounty))

########RECODING########
#SECTION D: SUBCOUNTY'S BASIC INFORMATION#
sc_merged$d15a[sc_merged$d15a==0.60000002] <- 60

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
sc_merged$baraza.H33 <- as.numeric(as.character(sc_merged$baraza.H33))
sc_merged$baraza.H34 <- as.numeric(as.character(sc_merged$baraza.H34))
sc_merged$baraza.H35 <- as.numeric(as.character(sc_merged$baraza.H35))
sc_merged$baraza.H36 <- as.numeric(as.character(sc_merged$baraza.H36))
sc_merged$baraza.H37 <- as.numeric(as.character(sc_merged$baraza.H37))
sc_merged$baraza.H38 <- as.numeric(as.character(sc_merged$baraza.H38))
sc_merged$baraza.H39 <- as.numeric(as.character(sc_merged$baraza.H39))
sc_merged$baraza.H40 <- as.numeric(as.character(sc_merged$baraza.H40))
sc_merged$baraza.H41 <- as.numeric(as.character(sc_merged$baraza.H41))
sc_merged$baraza.H42 <- as.numeric(as.character(sc_merged$baraza.H42))
sc_merged$baraza.H43 <- as.numeric(as.character(sc_merged$baraza.H43))
sc_merged$baraza.H44 <- as.numeric(as.character(sc_merged$baraza.H44))
sc_merged$baraza.H45 <- as.numeric(as.character(sc_merged$baraza.H45))
sc_merged$baraza.H46 <- as.numeric(as.character(sc_merged$baraza.H46))
sc_merged$baraza.H47 <- as.numeric(as.character(sc_merged$baraza.H47))
sc_merged$baraza.H49 <- as.numeric(as.character(sc_merged$baraza.H49))
sc_merged$baraza.H50 <- as.numeric(as.character(sc_merged$baraza.H50))
sc_merged$baraza.H51 <- as.numeric(as.character(sc_merged$baraza.H51))
sc_merged$baraza.H52 <- as.numeric(as.character(sc_merged$baraza.H52))
sc_merged$baraza.H53 <- as.numeric(as.character(sc_merged$baraza.H53))
sc_merged$baraza.H54 <- as.numeric(as.character(sc_merged$baraza.H54))
sc_merged$baraza.H55 <- as.numeric(as.character(sc_merged$baraza.H55))
sc_merged$baraza.H56 <- as.numeric(as.character(sc_merged$baraza.H56))
sc_merged$baraza.H57 <- as.numeric(as.character(sc_merged$baraza.H57))
sc_merged$baraza.H58 <- as.numeric(as.character(sc_merged$baraza.H58))
sc_merged$baraza.H59 <- as.numeric(as.character(sc_merged$baraza.H59))
sc_merged$baraza.H60 <- as.numeric(as.character(sc_merged$baraza.H60))
sc_merged$baraza.H61 <- as.numeric(as.character(sc_merged$baraza.H61))
sc_merged$baraza.H62 <- as.numeric(as.character(sc_merged$baraza.H62))
sc_merged$baraza.H63 <- as.numeric(as.character(sc_merged$baraza.H63))
sc_merged$baraza.H64 <- as.numeric(as.character(sc_merged$baraza.H64))
sc_merged$baraza.H68 <- as.numeric(as.character(sc_merged$baraza.H68))
sc_merged$baraza.H70 <- as.numeric(as.character(sc_merged$baraza.H70))
sc_merged$baraza.H71 <- as.numeric(as.character(sc_merged$baraza.H71))

#baraza.H5 has 115 NA's and is excluded from analysis
sc_merged$baraza.H24[sc_merged$baraza.H24==49] <- NA
sc_merged$baraza.H25[sc_merged$baraza.H25==49] <- NA
sc_merged$h181a <- replace(sc_merged$h181a, is.na(sc_merged$h181a), 0)
sc_merged$h181b <- replace(sc_merged$h181b, is.na(sc_merged$h181b), 0)
sc_merged$h181d <- replace(sc_merged$h181d, is.na(sc_merged$h181d), 0)
sc_merged$sum_h182a_h183a <- (sc_merged$h182a + sc_merged$h183a)
sc_merged$sum_h182b_h183b <- (sc_merged$h182b + sc_merged$h183b)
sc_merged$sum_h182d_h183d <- (sc_merged$h182d + sc_merged$h183d)
sc_merged$sum_h1123a_h1124a <- (sc_merged$h1123a + sc_merged$h1124a)
sc_merged$sum_h1123b_h1124b <- (sc_merged$h1123b + sc_merged$h1124b)
sc_merged$sum_h1123c_h1124c <- (sc_merged$h1123c + sc_merged$h1124c)
sc_merged$h1131f[sc_merged$h1131f=="Na"] <- NA
sc_merged$h1131h[sc_merged$h1131h=="I1"] <- 1
#baraza.H48 excluded from analysis because no HC2 has isolation room
sc_merged$"h1132f" <- replace(sc_merged$"h1132f", is.na(sc_merged$"h1132f"), 0) #otherwise loop does not run because 59 1's and 177 NA's
sc_merged$"h1132j" <- replace(sc_merged$"h1132j", is.na(sc_merged$"h1132j"), 0) #otherwise loop does not run
sc_merged$"h1132k" <- replace(sc_merged$"h1132k", is.na(sc_merged$"h1132k"), 0) #otherwise loop does not run
sc_merged$"h1132o" <- replace(sc_merged$"h1132o", is.na(sc_merged$"h1132o"), 0) #otherwise loop does not run
sc_merged$"h1132p" <- replace(sc_merged$"h1132p", is.na(sc_merged$"h1132p"), 0) #otherwise loop does not run
sc_merged$baraza.H65[sc_merged$baraza.H65==98] <- NA
sc_merged$baraza.H65_binary <- (sc_merged$baraza.H65 == 1)
sc_merged$h1171_binary <- (sc_merged$h1171 == "yes")
sc_merged$h119a[sc_merged$h119a==0.25] <- 25
sc_merged$h119b[sc_merged$h119b==0.2] <- 20
sc_merged$baraza.H69[sc_merged$baraza.H69==98] <- NA
sc_merged$baraza.H69_binary <- (sc_merged$baraza.H69 == 1)
sc_merged$h121_binary <- (sc_merged$h121 == "yes")
sc_merged$sum_h1221_to_h1224 <- (sc_merged$h1221 + sc_merged$h1222 + sc_merged$h1223 + sc_merged$h1224)
#baraza.H8, baraza.H11, baraza.H41, baraza.H42, baraza.H43, baraza.H44, baraza.H47 excluded from analysis because too many NA's

#SUBSECTION: WATER INFRASTRUCTURE#
sc_merged$baraza.H77[sc_merged$baraza.H77==98] <- NA
sc_merged$baraza.H77_binary <- (sc_merged$baraza.H77 == 1)
sc_merged$baraza.H78 <- as.numeric(as.character(sc_merged$baraza.H78))
sc_merged$baraza.H79 <- as.numeric(as.character(sc_merged$baraza.H79))
sc_merged$h216_binary <- (sc_merged$h216 == "yes")
sc_merged$sum_h217a_to_h217d <- (sc_merged$h217a + sc_merged$h217b + sc_merged$h217c + sc_merged$h217d)

#SUBSECTION H3: EDUCATION#
sc_merged$baraza.H88 <- as.numeric(as.character(sc_merged$baraza.H88))
sc_merged$baraza.H89 <- as.numeric(as.character(sc_merged$baraza.H89))
sc_merged$baraza.H90 <- as.numeric(as.character(sc_merged$baraza.H90))
sc_merged$baraza.H91 <- as.numeric(as.character(sc_merged$baraza.H91))
sc_merged$baraza.H93 <- as.numeric(as.character(sc_merged$baraza.H93))
sc_merged$baraza.H94 <- as.numeric(as.character(sc_merged$baraza.H94))
sc_merged$baraza.H95 <- as.numeric(as.character(sc_merged$baraza.H95))
sc_merged$baraza.H102 <- as.numeric(as.character(sc_merged$baraza.H102))
sc_merged$baraza.H103 <- as.numeric(as.character(sc_merged$baraza.H103))
sc_merged$baraza.H104 <- as.numeric(as.character(sc_merged$baraza.H104))
sc_merged$baraza.H105 <- as.numeric(as.character(sc_merged$baraza.H105))
sc_merged$baraza.H106 <- as.numeric(as.character(sc_merged$baraza.H106))
sc_merged$baraza.H107 <- as.numeric(as.character(sc_merged$baraza.H107))
sc_merged$baraza.H111 <- as.numeric(as.character(sc_merged$baraza.H111))
sc_merged$baraza.H112 <- as.numeric(as.character(sc_merged$baraza.H112))

sc_merged$h3111b[sc_merged$h3111b==0.1] <- 10
sc_merged$h3111a[sc_merged$h3111a==0.15000001] <- 15.000001
sc_merged$baraza.H109[sc_merged$baraza.H109==98] <- NA
sc_merged$baraza.H109_binary <- (sc_merged$baraza.H109 == 1)
sc_merged$h342_binary <- (sc_merged$h342 == "yes")
sc_merged$baraza.H110[sc_merged$baraza.H110==98] <- NA
sc_merged$baraza.H110_binary <- (sc_merged$baraza.H110 == 1)
sc_merged$h345[sc_merged$h345==""] <- NA
sc_merged$h345_binary <- (sc_merged$h345 == "yes")
sc_merged$sum_h346a_to_h346d <- (sc_merged$h346a + sc_merged$h346b + sc_merged$h346c + sc_merged$h346d)
sc_merged$baraza.H113_binary <- (sc_merged$baraza.H113 == 1)
sc_merged$h3511[sc_merged$h3511==""] <- NA
sc_merged$h3511_binary <- (sc_merged$h3511 == "yes")
sc_merged$baraza.H114[sc_merged$baraza.H114==98] <- NA
sc_merged$baraza.H114_binary <- (sc_merged$baraza.H114 == 1)
sc_merged$h361[sc_merged$h361==""] <- NA
sc_merged$h361_binary <- (sc_merged$h361 == "yes")

#SUBSECTION K: AGRICULTURE AND EXTENSION (NAADS; OPERATION WEALTH CREATION#
sc_merged$baraza.maleex.K1[sc_merged$baraza.maleex.K1==999] <- NA
sc_merged$baraza.maleex.K2[sc_merged$baraza.maleex.K2==999] <- NA
sc_merged$baraza.maleex.K4[sc_merged$baraza.maleex.K4==999] <- NA
sc_merged$baraza.maleex.K5[sc_merged$baraza.maleex.K5==999] <- NA
sc_merged$sum_maleex.K1_maleex.K2 <- (sc_merged$baraza.maleex.K1 + sc_merged$baraza.maleex.K2)
sc_merged$sum_maleex.K4_maleex.K5 <- (sc_merged$baraza.maleliv.K4 + sc_merged$baraza.maleliv.K5)
sc_merged$baraza.malec.K8[sc_merged$baraza.malec.K8==83] <- NA
#no analysis for baraza.H91 because 124 NA's in endline, 70 NA's in baseline

summary(sc_merged$h3233b)


########LOOPS########
#loop if NA cannot be interpreted as 0
outcomes <- c("baraza.km.D1","baraza.km.D2","baraza.km.D3","baraza.km.D4a","baraza.km.D4b","baraza.production.E1a_binary","baraza.health.E2a_binary","baraza.gender1.E3a_binary","baraza.works.E4a_binary","baraza.finance1.E5a_binary","baraza.E7_binary","baraza.meeting.F3","baraza.H1","baraza.H2","baraza.H2b","baraza.H12","baraza.H20","baraza.H29","baraza.H32","baraza.H65_binary","baraza.H69_binary","baraza.H86")
baseline_outcomes <- c("d12","d13","d14","d15a","d15b","e11a_binary","e11b_binary","e11c_binary","e11d_binary","e11e_binary","e13_binary","f14c","h11","h12","h15","h110","h1121c","h1125c","h1126c","h1171_binary","h121_binary","h3233a")
#df_ols <- array(NA,dim=c(6,3,length(outcomes)))
#df_ols <- array(NA,dim=c(7,3,length(outcomes)))
df_ols <- array(NA,dim=c(3,3,length(outcomes)))

sc_merged[outcomes[!duplicated(outcomes)]] <- lapply(sc_merged[outcomes[!duplicated(outcomes)]], function(x) replace(x, x == 999, NA) )
sc_merged[outcomes[!duplicated(outcomes)]] <- lapply(sc_merged[outcomes[!duplicated(outcomes)]], function(x) replace(x, x == "n/a", NA) )
sc_merged[baseline_outcomes[!duplicated(baseline_outcomes)]] <- lapply(sc_merged[baseline_outcomes[!duplicated(baseline_outcomes)]], function(x) replace(x, x == "", NA) )
sc_merged[baseline_outcomes[!duplicated(baseline_outcomes)]] <- lapply(sc_merged[baseline_outcomes[!duplicated(baseline_outcomes)]], function(x) replace(x, x == "N/A", NA) )

# sc_merged[outcomes] <- lapply(sc_merged[outcomes], function(x) replace(x, x == 999, NA) )
# sc_merged[outcomes] <- lapply(sc_merged[outcomes], function(x) replace(x, x == "n/a", NA) )
# sc_merged[baseline_outcomes] <- lapply(sc_merged[baseline_outcomes], function(x) replace(x, x == "", NA) )
# sc_merged[baseline_outcomes] <- lapply(sc_merged[baseline_outcomes], function(x) replace(x, x == "N/A", NA) )

### parallel computing for RI
cl <- makeCluster(detectCores(all.tests = FALSE, logical = TRUE))
registerDoParallel(cl)


for (i in 1:length(outcomes)) {
  print(i)
  ols <- lm(as.formula(paste(paste(outcomes[i],"information*deliberation+region.x",sep="~"),baseline_outcomes[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)

if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes, baseline_outcomes, subset(sc_merged, district_baraza == 0) , ctrls = "region.x", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <- RI_store$conf_2 
conf[3,4:5] <- RI_store$conf_3
res[2,5] <- RI_store$pval_2
res[3,5] <- RI_store$pval_3
}
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
if (RI_conf_switch) {
conf[6,4:5] <- RI_store$conf_1
res[6,5] <- RI_store$pval_1
}
  #df_ols[,1,i] <- c(res[6,1],res[6,2],res[6,5],res[6,6],conf[6,4],conf[6,5],nobs(ols))
  #df_ols[,1,i] <- c(res[6,1],res[6,2],res[6,5],conf[6,4],conf[6,5],nobs(ols))
  df_ols[,1,i] <- c(res[6,1],res[6,5],nobs(ols))
}









#loop if NA can be interpreted as 0
outcomes_NAcouldbe0 <- c("baraza.H3","baraza.H4","baraza.H6","baraza.H7","baraza.H9","baraza.H10","baraza.H13","baraza.H18","baraza.H19","baraza.H21","baraza.H22","baraza.H24","baraza.H25","baraza.H27","baraza.H28","baraza.H30","baraza.H31","baraza.H33","baraza.H34","baraza.H35","baraza.H36","baraza.H37","baraza.H38","baraza.H39","baraza.H40","baraza.H45","baraza.H46","baraza.H49","baraza.H50","baraza.H51","baraza.H52","baraza.H53","baraza.H54","baraza.H55","baraza.H56","baraza.H57","baraza.H58","baraza.H59","baraza.H60","baraza.H61","baraza.H62","baraza.H63","baraza.H64","baraza.H67","baraza.H67","baraza.H68","baraza.H70","baraza.H71","baraza.H72","baraza.H73","baraza.H73","baraza.H74","baraza.H75","baraza.H75","baraza.H77_binary","baraza.H78","baraza.H79","baraza.H67","baraza.H67","baraza.H80","baraza.H81","baraza.H82","baraza.H83","baraza.H84","baraza.H85","baraza.H87","baraza.H88","baraza.H89","baraza.H90","baraza.H92","baraza.H93","baraza.H94","baraza.H95","baraza.H96","baraza.H97","baraza.H98","baraza.H99","baraza.H100","baraza.H101","baraza.H102","baraza.H103","baraza.H104","baraza.H105","baraza.H108","baraza.H109_binary","baraza.H110_binary","baraza.H111","baraza.H112","baraza.H113_binary","baraza.H114_binary","baraza.H116","baraza.H116","baraza.maleex.K1","baraza.maleex.K2","sum_maleex.K1_maleex.K2","baraza.maleliv.K4","baraza.maleliv.K5","sum_maleex.K4_maleex.K5","baraza.K10")
baseline_outcomes_NAcouldbe0 <- c("h181a","h181b","sum_h182a_h183a","sum_h182b_h183b","h184a","h184b","h111","h1121a","h1121b","h1122a","h1122b","sum_h1123a_h1124a","sum_h1123b_h1124b","h1125a","h1125b","h1126a","h1126b","h1131a","h1131b","h1131c","h1131d","h1131e","h1131f","h1131g","h1131h","h1131m","h1131n","h1132a","h1132b","h1132c","h1132d","h1132e","h1132f","h1132g","h1132h","h1132i","h1132j","h1132k","h1132l","h1132m","h1132n","h1132o","h1132p","h119a","h119b","h1201","sum_h1221_to_h1224","h123","h213a","h213b","h213i","h213c","h213d","h213e","h216_binary","sum_h217a_to_h217d","h218","h119a","h119b","h3111b","h3111a","h321a","h322a","h323a","h3231a","h321b","h322b","h323b","h3231b","h321d","h322d","h321e","h322e","h331a","h331b","h331c","h331e","h331f","h331g","h332a","h332b","h332c","h332e","h341","h342_binary","h345_binary","sum_h346a_to_h346d","h356","h3511_binary","h361_binary","h386a","h386b","h42101a","h42101b","h422a","h42101_1","h42101_2","h422b","h441")
#df_ols_NAcouldbe0 <- array(NA,dim=c(6,10,length(outcomes_NAcouldbe0)))
df_ols_NAcouldbe0 <- array(NA,dim=c(3,9,length(outcomes_NAcouldbe0)))

sc_merged[outcomes_NAcouldbe0[!duplicated(outcomes_NAcouldbe0)]] <- lapply(sc_merged[outcomes_NAcouldbe0[!duplicated(outcomes_NAcouldbe0)]], function(x) replace(x, x == 999, NA) )
sc_merged[outcomes_NAcouldbe0[!duplicated(outcomes_NAcouldbe0)]] <- lapply(sc_merged[outcomes_NAcouldbe0[!duplicated(outcomes_NAcouldbe0)]], function(x) replace(x, x == "n/a", NA) )
sc_merged[baseline_outcomes_NAcouldbe0[!duplicated(baseline_outcomes_NAcouldbe0)]] <- lapply(sc_merged[baseline_outcomes_NAcouldbe0[!duplicated(baseline_outcomes_NAcouldbe0)]], function(x) replace(x, x == "", NA) )
sc_merged[baseline_outcomes_NAcouldbe0[!duplicated(baseline_outcomes_NAcouldbe0)]] <- lapply(sc_merged[baseline_outcomes_NAcouldbe0[!duplicated(baseline_outcomes_NAcouldbe0)]], function(x) replace(x, x == "N/A", NA) )

#sc_merged[outcomes_NAcouldbe0] <- lapply(sc_merged[outcomes_NAcouldbe0], function(x) replace(x, x == 999, NA) )
#sc_merged[outcomes_NAcouldbe0] <- lapply(sc_merged[outcomes_NAcouldbe0], function(x) replace(x, x == "n/a", NA) )
#sc_merged[baseline_outcomes_NAcouldbe0] <- lapply(sc_merged[baseline_outcomes_NAcouldbe0], function(x) replace(x, x == "", NA) )
#sc_merged[baseline_outcomes_NAcouldbe0] <- lapply(sc_merged[baseline_outcomes_NAcouldbe0], function(x) replace(x, x == "N/A", NA) )


for (i in 1:length(outcomes_NAcouldbe0)) {
  #A: NAs remain NAs
  print(i)
  ols <- lm(as.formula(paste(paste(outcomes_NAcouldbe0[i],"information*deliberation+region.x",sep="~"),baseline_outcomes_NAcouldbe0[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes_NAcouldbe0, baseline_outcomes_NAcouldbe0, subset(sc_merged, district_baraza == 0) , ctrls = "region.x", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <- RI_store$conf_2 
conf[3,4:5] <- RI_store$conf_3
res[2,5] <- RI_store$pval_2
res[3,5] <- RI_store$pval_3
}
  #df_ols_NAcouldbe0[,2,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5],nobs(ols))
  df_ols_NAcouldbe0[,2,i] <- c(res[2,1],res[2,5],nobs(ols))
  #df_ols_NAcouldbe0[,3,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5],nobs(ols))
  df_ols_NAcouldbe0[,3,i] <- c(res[3,1],res[3,5],nobs(ols))
  
  ols <- lm(as.formula(paste(paste(outcomes_NAcouldbe0[i],"information:deliberation+region.x",sep="~"),baseline_outcomes_NAcouldbe0[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
conf[6,4:5] <- RI_store$conf_1
res[6,5] <- RI_store$pval_1
}
  #df_ols_NAcouldbe0[,1,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5],nobs(ols))
  df_ols_NAcouldbe0[,1,i] <- c(res[6,1],res[6,5],nobs(ols))
  
  #B: recode NAs in baseline data as 0
  sc_merged[baseline_outcomes_NAcouldbe0[i]] <- replace(sc_merged[baseline_outcomes_NAcouldbe0[i]], is.na(sc_merged[baseline_outcomes_NAcouldbe0[i]]), 0)
  
  print(i)
  ols <- lm(as.formula(paste(paste(outcomes_NAcouldbe0[i],"information*deliberation+region.x",sep="~"),baseline_outcomes_NAcouldbe0[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes_NAcouldbe0, baseline_outcomes_NAcouldbe0, subset(sc_merged, district_baraza == 0) , ctrls = "region.x", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <- RI_store$conf_2 
conf[3,4:5] <- RI_store$conf_3
res[2,5] <- RI_store$pval_2
res[3,5] <- RI_store$pval_3
}
  #df_ols_NAcouldbe0[,5,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5],nobs(ols))
  df_ols_NAcouldbe0[,5,i] <- c(res[2,1],res[2,5],nobs(ols))
  #df_ols_NAcouldbe0[,6,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5],nobs(ols))
  df_ols_NAcouldbe0[,6,i] <- c(res[3,1],res[3,5],nobs(ols))
  ols <- lm(as.formula(paste(paste(outcomes_NAcouldbe0[i],"information:deliberation+region.x",sep="~"),baseline_outcomes_NAcouldbe0[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
conf[6,4:5] <- RI_store$conf_1
res[6,5] <- RI_store$pval_1
}
  #df_ols_NAcouldbe0[,4,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5],nobs(ols))
  df_ols_NAcouldbe0[,4,i] <- c(res[6,1],res[6,5],nobs(ols))
    
  #C: recode NAs in baseline and endline data as 0
  sc_merged[outcomes_NAcouldbe0[i]] <- replace(sc_merged[outcomes_NAcouldbe0[i]], is.na(sc_merged[outcomes_NAcouldbe0[i]]), 0)
  
  print(i)
  ols <- lm(as.formula(paste(paste(outcomes_NAcouldbe0[i],"information*deliberation+region.x",sep="~"),baseline_outcomes_NAcouldbe0[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
RI_store <- RI_conf_sc(i,outcomes_NAcouldbe0, baseline_outcomes_NAcouldbe0, subset(sc_merged, district_baraza == 0) , ctrls = "region.x", nr_repl = glob_repli, sig = glob_sig)
conf[2,4:5] <- RI_store$conf_2 
conf[3,4:5] <- RI_store$conf_3
res[2,5] <- RI_store$pval_2
res[3,5] <- RI_store$pval_3
}
  #df_ols_NAcouldbe0[,8,i] <- c(res[2,1],res[2,2],res[2,5], conf[2,4],conf[2,5],nobs(ols))
  df_ols_NAcouldbe0[,8,i] <- c(res[2,1],res[2,5],nobs(ols))
  #df_ols_NAcouldbe0[,9,i] <- c(res[3,1],res[3,2],res[3,5], conf[3,4],conf[3,5],nobs(ols))
  df_ols_NAcouldbe0[,9,i] <- c(res[3,1],res[3,5],nobs(ols))
  
  ols <- lm(as.formula(paste(paste(outcomes_NAcouldbe0[i],"information:deliberation+region.x",sep="~"),baseline_outcomes_NAcouldbe0[i],sep="+")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
if (RI_conf_switch) {
conf[6,4:5] <- RI_store$conf_1
res[6,5] <- RI_store$pval_1
}
  #df_ols_NAcouldbe0[,7,i] <- c(res[5,1],res[5,2],res[5,5], conf[5,4],conf[5,5],nobs(ols))
  df_ols_NAcouldbe0[,7,i] <- c(res[6,1],res[6,5],nobs(ols))
  
}

#loop if baseline variable is missing/ contains too many NAs
#baraza.H15, baraza.H16, baraza.H17: there is no equivalent in baseline data
#baraza.H23: h1122c is baseline variable but contains 212 NA's
#baraza.meeting.F1: f14a 106 NAs
#baraza.meeting.F2: f14b 112 NAs
#baraza.meeting.F4: f14d 105 NAs
#baraza.H8: sum_h182d_h183d 109 NAs
#baraza.H11: h184d 159 NAs
#baraza.H26: sum_h1123c_h1124c 81 NAs
#baraza.H41: no corresponding baseline variable
#baraza.H42 and baraza.H43 and baraza.H44: no corresponding baseline variable
#baraza.H47: h1131o 235 NAs
#baraza.H76: no corresponding baseline variable
#baraza.H106 and baraza.H107: no corresponding baseline variable
#baraza.maleex.K3, baraza.maleliv.K6, baraza.malec.K9: h421a, h421b, h421c have too many NAs

outcomes_nobaseline <- c("baraza.H15","baraza.H16","baraza.H17","baraza.H23","baraza.meeting.F1","baraza.meeting.F2","baraza.meeting.F4","baraza.H26","baraza.H76","baraza.H106","baraza.H107","baraza.malec.K7","baraza.malec.K8","baraza.maleex.K3","baraza.maleliv.K6","baraza.malec.K9")
#df_ols_nobaseline <- array(NA,dim=c(6,3,length(outcomes_nobaseline)))
#df_ols_nobaseline <- array(NA,dim=c(7,3,length(outcomes_nobaseline)))
df_ols_nobaseline <- array(NA,dim=c(3,3,length(outcomes_nobaseline)))

sc_merged[outcomes_nobaseline[!duplicated(outcomes_nobaseline)]] <- lapply(sc_merged[outcomes_nobaseline[!duplicated(outcomes_nobaseline)]], function(x) replace(x, x == 999, NA) )
sc_merged[outcomes_nobaseline[!duplicated(outcomes_nobaseline)]] <- lapply(sc_merged[outcomes_nobaseline[!duplicated(outcomes_nobaseline)]], function(x) replace(x, x == "n/a", NA) )

#sc_merged[outcomes_nobaseline] <- lapply(sc_merged[outcomes_nobaseline], function(x) replace(x, x == 999, NA) )
#sc_merged[outcomes_nobaseline] <- lapply(sc_merged[outcomes_nobaseline], function(x) replace(x, x == "n/a", NA) )

for (i in 1:length(outcomes_nobaseline)) {
  print(i)
  ols <- lm(as.formula(paste(outcomes_nobaseline[i],"information*deliberation+region.x",sep="~")), data=sc_merged[sc_merged$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols_nobaseline[,2,i] <- c(res[2,1],res[2,2],res[2,5],res[2,6],conf[2,4],conf[2,5],nobs(ols))
  #df_ols_nobaseline[,2,i] <- c(res[2,1],res[2,2],res[2,5],conf[2,4],conf[2,5],nobs(ols))
  df_ols_nobaseline[,2,i] <- c(res[2,1],res[2,5],nobs(ols))
  #df_ols_nobaseline[,3,i] <- c(res[3,1],res[3,2],res[3,5],res[3,6],conf[3,4],conf[3,5],nobs(ols))
  #df_ols_nobaseline[,3,i] <- c(res[3,1],res[3,2],res[3,5],conf[3,4],conf[3,5],nobs(ols))
  df_ols_nobaseline[,3,i] <- c(res[3,1],res[3,5],nobs(ols))
  ols <- lm(as.formula(paste(outcomes_nobaseline[i],"information:deliberation+region.x",sep="~")), data=sc_merged[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_merged$clusterID[sc_merged$district_baraza == 0 & (sc_merged$information == sc_merged$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  #df_ols_nobaseline[,1,i] <- c(res[6,1],res[6,2],res[6,5],res[6,6],conf[6,4],conf[6,5],nobs(ols))
  #df_ols_nobaseline[,1,i] <- c(res[6,1],res[6,2],res[6,5],conf[6,4],conf[6,5],nobs(ols))
  df_ols_nobaseline[,1,i] <- c(res[6,1],res[6,5],nobs(ols))
}
