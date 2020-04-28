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


#load sc endline and baseline data
sc_baseline <- read.csv(paste(path,"data/public/sc_level_baseline.csv", sep ="/"))
sc_baseline$actor <- NA
sc_baseline$actor[sc_baseline$designation %in% c("LC3_chair","Vice Chairperson LCIII")] <- "politician"
sc_baseline$actor[sc_baseline$designation %in% c("Subcounty Chief/Town Clerk","Parish Chief","Acting Subcounty Chief/Town Clerk","Deputy subcounty chief/Town Clerk","Community Development Officer")] <- "civil servant"

# sc_baseline$subcounty <- as.character(sc_baseline$subcounty)
# sc_baseline$subcounty[sc_baseline$subcounty == "SEMBABULE TC"] <- "SEMBABULE_TC"
# sc_baseline$subcounty[sc_baseline$subcounty ==  "NTUSI"] <- "NTUUSI"

#load treatment assignment
treats <- read.csv(paste(path,"data/public/treats.csv", sep ="/"))
#merge in treatments
sc_baseline <- merge(treats, sc_baseline, by.x=c("district","subcounty"), by.y=c("district","subcounty"))


# setdiff(sc_endline$subcounty,sc_baseline$subcounty)
# [1] "HARUGALI"     "KIGULYA": these two were apparenty not inteviewed during baseline? 
#     "NTUUSI"       "SEMBABULE_TC":: these two have spelling errors

# ########RECODING########
# #SECTION D: SUBCOUNTY'S BASIC INFORMATION#
# sc_merged$d15a[sc_merged$d15a==0.60000002] <- 60
# 
# #SECTION E: GOVERNMENT AND COMMUNITY BODIES IN THE SUBCOUNTY#
# sc_merged$baraza.production.E1a_binary <- (sc_merged$baraza.production.E1a == 1)
# sc_merged$baraza.health.E2a_binary <- (sc_merged$baraza.health.E2a == 1)
# sc_merged$baraza.gender1.E3a_binary <- (sc_merged$baraza.gender1.E3a == 1)
# sc_merged$baraza.works.E4a_binary <- (sc_merged$baraza.works.E4a == 1)
# sc_merged$baraza.finance1.E5a_binary <- (sc_merged$baraza.finance1.E5a == 1)
# sc_merged$baraza.E7[sc_merged$baraza.E7==96] <- NA
# sc_merged$baraza.E7_binary <- (sc_merged$baraza.E7 > 2) #median is 2
# sc_merged$e11a_binary <- (sc_merged$e11a == "Yes")
# sc_merged$e11b_binary <- (sc_merged$e11b == "Yes")
# sc_merged$e11c_binary <- (sc_merged$e11c == "Yes")
# sc_merged$e11d_binary <- (sc_merged$e11d == "Yes")
# sc_merged$e11e_binary <- (sc_merged$e11e == "Yes")
sc_baseline$e13_binary <- (sc_baseline$e13 %in% c("After every two months", "Monthly","Quarterly"))
# 
# #SECTION F: COMMUNITY MEETINGS HELD IN THE SUBCOUNTY#
# sc_merged$f14a[sc_merged$f14a==2015] <- NA
# 
# #SUBSECTION H1 - HEALTH#
# sc_merged$baraza.H3 <- as.numeric(as.character(sc_merged$baraza.H3))
# sc_merged$baraza.H4 <- as.numeric(as.character(sc_merged$baraza.H4))
# sc_merged$baraza.H5 <- as.numeric(as.character(sc_merged$baraza.H5))
# sc_merged$baraza.H6 <- as.numeric(as.character(sc_merged$baraza.H6))
# sc_merged$baraza.H7 <- as.numeric(as.character(sc_merged$baraza.H7))
# sc_merged$baraza.H8 <- as.numeric(as.character(sc_merged$baraza.H8))
# sc_merged$baraza.H9 <- as.numeric(as.character(sc_merged$baraza.H9))
# sc_merged$baraza.H10 <- as.numeric(as.character(sc_merged$baraza.H10))
# sc_merged$baraza.H11 <- as.numeric(as.character(sc_merged$baraza.H11))
# sc_merged$baraza.H13 <- as.numeric(as.character(sc_merged$baraza.H13))
# sc_merged$baraza.H15 <- as.numeric(as.character(sc_merged$baraza.H15))
# sc_merged$baraza.H16 <- as.numeric(as.character(sc_merged$baraza.H16))
# sc_merged$baraza.H17 <- as.numeric(as.character(sc_merged$baraza.H17))
# sc_merged$baraza.H18 <- as.numeric(as.character(sc_merged$baraza.H18))
# sc_merged$baraza.H19 <- as.numeric(as.character(sc_merged$baraza.H19))
# sc_merged$baraza.H20 <- as.numeric(as.character(sc_merged$baraza.H20))
# sc_merged$baraza.H21 <- as.numeric(as.character(sc_merged$baraza.H21))
# sc_merged$baraza.H22 <- as.numeric(as.character(sc_merged$baraza.H22))
# sc_merged$baraza.H23 <- as.numeric(as.character(sc_merged$baraza.H23))
# sc_merged$baraza.H24 <- as.numeric(as.character(sc_merged$baraza.H24))
# sc_merged$baraza.H25 <- as.numeric(as.character(sc_merged$baraza.H25))
# sc_merged$baraza.H26 <- as.numeric(as.character(sc_merged$baraza.H26))
# sc_merged$baraza.H27 <- as.numeric(as.character(sc_merged$baraza.H27))
# sc_merged$baraza.H28 <- as.numeric(as.character(sc_merged$baraza.H28))
# sc_merged$baraza.H29 <- as.numeric(as.character(sc_merged$baraza.H29))
# sc_merged$baraza.H30 <- as.numeric(as.character(sc_merged$baraza.H30))
# sc_merged$baraza.H31 <- as.numeric(as.character(sc_merged$baraza.H31))
# sc_merged$baraza.H32 <- as.numeric(as.character(sc_merged$baraza.H32))
# sc_merged$baraza.H33 <- as.numeric(as.character(sc_merged$baraza.H33))
# sc_merged$baraza.H34 <- as.numeric(as.character(sc_merged$baraza.H34))
# sc_merged$baraza.H35 <- as.numeric(as.character(sc_merged$baraza.H35))
# sc_merged$baraza.H36 <- as.numeric(as.character(sc_merged$baraza.H36))
# sc_merged$baraza.H37 <- as.numeric(as.character(sc_merged$baraza.H37))
# sc_merged$baraza.H38 <- as.numeric(as.character(sc_merged$baraza.H38))
# sc_merged$baraza.H39 <- as.numeric(as.character(sc_merged$baraza.H39))
# sc_merged$baraza.H40 <- as.numeric(as.character(sc_merged$baraza.H40))
# sc_merged$baraza.H41 <- as.numeric(as.character(sc_merged$baraza.H41))
# sc_merged$baraza.H42 <- as.numeric(as.character(sc_merged$baraza.H42))
# sc_merged$baraza.H43 <- as.numeric(as.character(sc_merged$baraza.H43))
# sc_merged$baraza.H44 <- as.numeric(as.character(sc_merged$baraza.H44))
# sc_merged$baraza.H45 <- as.numeric(as.character(sc_merged$baraza.H45))
# sc_merged$baraza.H46 <- as.numeric(as.character(sc_merged$baraza.H46))
# sc_merged$baraza.H47 <- as.numeric(as.character(sc_merged$baraza.H47))
# sc_merged$baraza.H49 <- as.numeric(as.character(sc_merged$baraza.H49))
# sc_merged$baraza.H50 <- as.numeric(as.character(sc_merged$baraza.H50))
# sc_merged$baraza.H51 <- as.numeric(as.character(sc_merged$baraza.H51))
# sc_merged$baraza.H52 <- as.numeric(as.character(sc_merged$baraza.H52))
# sc_merged$baraza.H53 <- as.numeric(as.character(sc_merged$baraza.H53))
# sc_merged$baraza.H54 <- as.numeric(as.character(sc_merged$baraza.H54))
# sc_merged$baraza.H55 <- as.numeric(as.character(sc_merged$baraza.H55))
# sc_merged$baraza.H56 <- as.numeric(as.character(sc_merged$baraza.H56))
# sc_merged$baraza.H57 <- as.numeric(as.character(sc_merged$baraza.H57))
# sc_merged$baraza.H58 <- as.numeric(as.character(sc_merged$baraza.H58))
# sc_merged$baraza.H59 <- as.numeric(as.character(sc_merged$baraza.H59))
# sc_merged$baraza.H60 <- as.numeric(as.character(sc_merged$baraza.H60))
# sc_merged$baraza.H61 <- as.numeric(as.character(sc_merged$baraza.H61))
# sc_merged$baraza.H62 <- as.numeric(as.character(sc_merged$baraza.H62))
# sc_merged$baraza.H63 <- as.numeric(as.character(sc_merged$baraza.H63))
# sc_merged$baraza.H64 <- as.numeric(as.character(sc_merged$baraza.H64))
# sc_merged$baraza.H68 <- as.numeric(as.character(sc_merged$baraza.H68))
# sc_merged$baraza.H70 <- as.numeric(as.character(sc_merged$baraza.H70))
# sc_merged$baraza.H71 <- as.numeric(as.character(sc_merged$baraza.H71))
# 
# #baraza.H5 has 115 NA's and is excluded from analysis
# sc_merged$baraza.H24[sc_merged$baraza.H24==49] <- NA
# sc_merged$baraza.H25[sc_merged$baraza.H25==49] <- NA
# sc_merged$h181a <- replace(sc_merged$h181a, is.na(sc_merged$h181a), 0)
# sc_merged$h181b <- replace(sc_merged$h181b, is.na(sc_merged$h181b), 0)
# sc_merged$h181d <- replace(sc_merged$h181d, is.na(sc_merged$h181d), 0)
sc_baseline$sum_h182a_h183a <- (sc_baseline$h182a + sc_baseline$h183a)
# sc_merged$sum_h182b_h183b <- (sc_merged$h182b + sc_merged$h183b)
# sc_merged$sum_h182d_h183d <- (sc_merged$h182d + sc_merged$h183d)
# sc_merged$sum_h1123a_h1124a <- (sc_merged$h1123a + sc_merged$h1124a)
# sc_merged$sum_h1123b_h1124b <- (sc_merged$h1123b + sc_merged$h1124b)
sc_baseline$sum_h1123c_h1124c <- (sc_baseline$h1123c + sc_baseline$h1124c)
# sc_merged$h1131f[sc_merged$h1131f=="Na"] <- NA
# sc_merged$h1131h[sc_merged$h1131h=="I1"] <- 1
# #baraza.H48 excluded from analysis because no HC2 has isolation room
# sc_merged$"h1132f" <- replace(sc_merged$"h1132f", is.na(sc_merged$"h1132f"), 0) #otherwise loop does not run because 59 1's and 177 NA's
# sc_merged$"h1132j" <- replace(sc_merged$"h1132j", is.na(sc_merged$"h1132j"), 0) #otherwise loop does not run
# sc_merged$"h1132k" <- replace(sc_merged$"h1132k", is.na(sc_merged$"h1132k"), 0) #otherwise loop does not run
# sc_merged$"h1132o" <- replace(sc_merged$"h1132o", is.na(sc_merged$"h1132o"), 0) #otherwise loop does not run
# sc_merged$"h1132p" <- replace(sc_merged$"h1132p", is.na(sc_merged$"h1132p"), 0) #otherwise loop does not run
# sc_merged$baraza.H65[sc_merged$baraza.H65==98] <- NA
# sc_merged$baraza.H65_binary <- (sc_merged$baraza.H65 == 1)
# sc_merged$h1171_binary <- (sc_merged$h1171 == "yes")
sc_baseline$h119a[sc_baseline$h119a==0.25] <- 25
sc_baseline$h119b[sc_baseline$h119b==0.2] <- 20
# sc_merged$baraza.H69[sc_merged$baraza.H69==98] <- NA
# sc_merged$baraza.H69_binary <- (sc_merged$baraza.H69 == 1)
# sc_merged$h121_binary <- (sc_merged$h121 == "yes")
# sc_merged$sum_h1221_to_h1224 <- (sc_merged$h1221 + sc_merged$h1222 + sc_merged$h1223 + sc_merged$h1224)
# #baraza.H8, baraza.H11, baraza.H41, baraza.H42, baraza.H43, baraza.H44, baraza.H47 excluded from analysis because too many NA's
# 
# #SUBSECTION: WATER INFRASTRUCTURE#
# sc_merged$baraza.H77[sc_merged$baraza.H77==98] <- NA
# sc_merged$baraza.H77_binary <- (sc_merged$baraza.H77 == 1)
# sc_merged$baraza.H78 <- as.numeric(as.character(sc_merged$baraza.H78))
# sc_merged$baraza.H79 <- as.numeric(as.character(sc_merged$baraza.H79))
# sc_merged$h216_binary <- (sc_merged$h216 == "yes")
# sc_merged$sum_h217a_to_h217d <- (sc_merged$h217a + sc_merged$h217b + sc_merged$h217c + sc_merged$h217d)
# 
# #SUBSECTION H3: EDUCATION#
# sc_merged$baraza.H88 <- as.numeric(as.character(sc_merged$baraza.H88))
# sc_merged$baraza.H89 <- as.numeric(as.character(sc_merged$baraza.H89))
# sc_merged$baraza.H90 <- as.numeric(as.character(sc_merged$baraza.H90))
# sc_merged$baraza.H91 <- as.numeric(as.character(sc_merged$baraza.H91))
# sc_merged$baraza.H93 <- as.numeric(as.character(sc_merged$baraza.H93))
# sc_merged$baraza.H94 <- as.numeric(as.character(sc_merged$baraza.H94))
# sc_merged$baraza.H95 <- as.numeric(as.character(sc_merged$baraza.H95))
# sc_merged$baraza.H102 <- as.numeric(as.character(sc_merged$baraza.H102))
# sc_merged$baraza.H103 <- as.numeric(as.character(sc_merged$baraza.H103))
# sc_merged$baraza.H104 <- as.numeric(as.character(sc_merged$baraza.H104))
# sc_merged$baraza.H105 <- as.numeric(as.character(sc_merged$baraza.H105))
# sc_merged$baraza.H106 <- as.numeric(as.character(sc_merged$baraza.H106))
# sc_merged$baraza.H107 <- as.numeric(as.character(sc_merged$baraza.H107))
# sc_merged$baraza.H111 <- as.numeric(as.character(sc_merged$baraza.H111))
# sc_merged$baraza.H112 <- as.numeric(as.character(sc_merged$baraza.H112))
# 
# sc_merged$h3111b[sc_merged$h3111b==0.1] <- 10
# sc_merged$h3111a[sc_merged$h3111a==0.15000001] <- 15.000001
# sc_merged$baraza.H109[sc_merged$baraza.H109==98] <- NA
# sc_merged$baraza.H109_binary <- (sc_merged$baraza.H109 == 1)
# sc_merged$h342_binary <- (sc_merged$h342 == "yes")
# sc_merged$baraza.H110[sc_merged$baraza.H110==98] <- NA
# sc_merged$baraza.H110_binary <- (sc_merged$baraza.H110 == 1)
# sc_merged$h345[sc_merged$h345==""] <- NA
# sc_merged$h345_binary <- (sc_merged$h345 == "yes")
# sc_merged$sum_h346a_to_h346d <- (sc_merged$h346a + sc_merged$h346b + sc_merged$h346c + sc_merged$h346d)
# sc_merged$baraza.H113_binary <- (sc_merged$baraza.H113 == 1)
# sc_merged$h3511[sc_merged$h3511==""] <- NA
# sc_merged$h3511_binary <- (sc_merged$h3511 == "yes")
# sc_merged$baraza.H114[sc_merged$baraza.H114==98] <- NA
# sc_merged$baraza.H114_binary <- (sc_merged$baraza.H114 == 1)
# sc_merged$h361[sc_merged$h361==""] <- NA
# sc_merged$h361_binary <- (sc_merged$h361 == "yes")
# 
# #SUBSECTION K: AGRICULTURE AND EXTENSION (NAADS; OPERATION WEALTH CREATION#
# sc_merged$baraza.maleex.K1[sc_merged$baraza.maleex.K1==999] <- NA
# sc_merged$baraza.maleex.K2[sc_merged$baraza.maleex.K2==999] <- NA
# sc_merged$baraza.maleliv.K4[sc_merged$baraza.maleliv.K4==999] <- NA
# sc_merged$baraza.maleliv.K5[sc_merged$baraza.maleliv.K5==999] <- NA
# sc_merged$sum_maleex.K1_maleex.K2 <- (sc_merged$baraza.maleex.K1 + sc_merged$baraza.maleex.K2)
# sc_merged$sum_maleliv.K4_maleliv.K5 <- (sc_merged$baraza.maleliv.K4 + sc_merged$baraza.maleliv.K5)
# sc_merged$baraza.malec.K8[sc_merged$baraza.malec.K8==83] <- NA
# #no analysis for baraza.H91 because 124 NA's in endline, 70 NA's in baseline
# sc_merged$improved_livestock_breedsh451_br[sc_merged$improved_livestock_breedsh451_br==0.1] <- 10
# sc_merged$improved_livestock_breedsh451_br[sc_merged$improved_livestock_breedsh451_br==0.5] <- 50
# sc_merged$baraza.K16[sc_merged$baraza.K16==98] <- NA
# sc_merged$baraza.K16_binary <- (sc_merged$baraza.K16 == 1)
# sc_merged$baraza.K17[sc_merged$baraza.K17==98] <- NA
# sc_merged$baraza.K17_binary <- (sc_merged$baraza.K17 == 1)
# sc_merged$h48s_binary <- (sc_merged$h48s == "yes")
# sc_merged$h491_binary <- (sc_merged$h491 == "yes")
# sc_merged$sum_seed1_seed2_seed3 <- (sc_merged$seed_1h4112_ns + sc_merged$seed_2h4112_nc + sc_merged$seed_3h4112_seed3)
# #sc_merged$improved_goat_noteh4112_goat excluded because 235 NA's, 1 answer
# #sc_merged$improved_breed_noteh4112_breed excluded because 228 NA's, 8 answers
# #sc_merged$improved_pig_noteh4112_pig excluded because 228 NA's, 8 answers
# sc_merged$baraza.K24[sc_merged$baraza.K24==98] <- NA
# sc_merged$baraza.K24_binary <- (sc_merged$baraza.K24 == 1)
# sc_merged$baraza.K24b <- as.numeric(as.character(sc_merged$baraza.K24b))
# sc_merged$baraza.K25 <- as.numeric(as.character(sc_merged$baraza.K25))
# sc_merged$h4115_binary <- (sc_merged$h4115 == "yes")
# sc_merged$h4116from_individuals[sc_merged$h4116from_individuals==3000] <- NA
# sc_merged$sum_from_farmersforum_ngos_individuals_other <- (sc_merged$h4116from_farmers_forum + sc_merged$h4116from_ngos + sc_merged$h4116from_individuals + sc_merged$h4116other_complaint)
# sc_merged$h4117[sc_merged$h4117==300] <- NA
# 
# #SECTION I:  PRIORITIES/ RANKINGS OF PROBLEMS#
# var_rank <- c("i_2i_2_1","i_2i_2_2","i_2i_2_3","i_2i_2_4","i_2i_2_5","i_2i_2_6","i_2i_2_7","i_2i_2_8","i_2i_2_9","i_2i_2_10","i_2i_2_11","i_2i_2_12","i_2i_2_13","i_2i_2_14")
# sc_merged[var_rank] <- lapply(sc_merged[var_rank],function(x)  as.character(unlist(x)) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "0_Not_applicable", NA) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "10_An_extreme_problem", 10) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "1_Not_a_problem", 1) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "2_2", 2) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "3_3", 3) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "4_4", 4) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "5_5", 5) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "6_6", 6) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "7_7", 7) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "8_8", 8) )
# sc_merged[var_rank] <- lapply(sc_merged[var_rank], function(x) replace(x, x == "9_9", 9) )
# sc_merged$i_2i_2_1 <- as.numeric(as.character(sc_merged$i_2i_2_1))
# sc_merged$i_2i_2_2 <- as.numeric(as.character(sc_merged$i_2i_2_2))
# sc_merged$i_2i_2_3 <- as.numeric(as.character(sc_merged$i_2i_2_3))
# sc_merged$i_2i_2_4 <- as.numeric(as.character(sc_merged$i_2i_2_4))
# sc_merged$i_2i_2_5 <- as.numeric(as.character(sc_merged$i_2i_2_5))
# sc_merged$i_2i_2_6 <- as.numeric(as.character(sc_merged$i_2i_2_6))
# sc_merged$i_2i_2_7 <- as.numeric(as.character(sc_merged$i_2i_2_7))
# sc_merged$i_2i_2_8 <- as.numeric(as.character(sc_merged$i_2i_2_8))
# sc_merged$i_2i_2_9 <- as.numeric(as.character(sc_merged$i_2i_2_9))
# sc_merged$i_2i_2_10 <- as.numeric(as.character(sc_merged$i_2i_2_10))
# sc_merged$i_2i_2_11 <- as.numeric(as.character(sc_merged$i_2i_2_11))
# sc_merged$i_2i_2_12 <- as.numeric(as.character(sc_merged$i_2i_2_12))
# sc_merged$i_2i_2_13 <- as.numeric(as.character(sc_merged$i_2i_2_13))
# sc_merged$i_2i_2_14 <- as.numeric(as.character(sc_merged$i_2i_2_14))
# # sc_merged$baraza.L1_binary <- (sc_merged$baraza.L1 > 5)
# # sc_merged$baraza.L2_binary <- (sc_merged$baraza.L2 > 5)
# # sc_merged$baraza.L3_binary <- (sc_merged$baraza.L3 > 5)
# # sc_merged$baraza.L4_binary <- (sc_merged$baraza.L4 > 5)
# # sc_merged$baraza.L5_binary <- (sc_merged$baraza.L5 > 5)
# # sc_merged$baraza.L6_binary <- (sc_merged$baraza.L6 > 5)
# # sc_merged$baraza.L7_binary <- (sc_merged$baraza.L7 > 5)
# # sc_merged$baraza.L8_binary <- (sc_merged$baraza.L8 > 5)
# # sc_merged$baraza.L9_binary <- (sc_merged$baraza.L9 > 5)
# # sc_merged$baraza.L10_binary <- (sc_merged$baraza.L10 > 5)
# # sc_merged$baraza.L11_binary <- (sc_merged$baraza.L11 > 5)
# # sc_merged$baraza.L12_binary <- (sc_merged$baraza.L12 > 5)
# # sc_merged$baraza.L13_binary <- (sc_merged$baraza.L13 > 5)
# # sc_merged$baraza.L14_binary <- (sc_merged$baraza.L14 > 5)
# # sc_merged$i_2i_2_1_binary <- (sc_merged$i_2i_2_1 > 5)
# # sc_merged$i_2i_2_2_binary <- (sc_merged$i_2i_2_2 > 5)
# # sc_merged$i_2i_2_3_binary <- (sc_merged$i_2i_2_3 > 5)
# # sc_merged$i_2i_2_4_binary <- (sc_merged$i_2i_2_4 > 5)
# # sc_merged$i_2i_2_5_binary <- (sc_merged$i_2i_2_5 > 5)
# # sc_merged$i_2i_2_6_binary <- (sc_merged$i_2i_2_6 > 5)
# # sc_merged$i_2i_2_7_binary <- (sc_merged$i_2i_2_7 > 5)
# # sc_merged$i_2i_2_8_binary <- (sc_merged$i_2i_2_8 > 5)
# # sc_merged$i_2i_2_9_binary <- (sc_merged$i_2i_2_9 > 5)
# # sc_merged$i_2i_2_10_binary <- (sc_merged$i_2i_2_10 > 5)
# # sc_merged$i_2i_2_11_binary <- (sc_merged$i_2i_2_11 > 5)
# # sc_merged$i_2i_2_12_binary <- (sc_merged$i_2i_2_12 > 5)
# # sc_merged$i_2i_2_13_binary <- (sc_merged$i_2i_2_13 > 5)
# # sc_merged$i_2i_2_14_binary <- (sc_merged$i_2i_2_14 > 5)

########LOOPS########
#loop if NA cannot be interpreted as 0
outcomes <- c("e13_binary","h119a","h119b","d13","d14","h42101a","improved_seedsh451_seeds","sum_h1123c_h1124c","sum_h182a_h183a","h322a","h322b")
df_ancova <- array("",dim=c(2,9,length(outcomes)))

sc_baseline[outcomes[!duplicated(outcomes)]] <- lapply(sc_baseline[outcomes[!duplicated(outcomes)]], function(x) replace(x, x == 999, NA) )
sc_baseline[outcomes[!duplicated(outcomes)]] <- lapply(sc_baseline[outcomes[!duplicated(outcomes)]], function(x) replace(x, x == "n/a", NA) )
sc_baseline <- sc_baseline %>%  mutate(clusterID = group_indices(., district, subcounty))

for (i in 1:length(outcomes)) {
  print(i)
  
  df_ancova[1,1,i] <- round(mean(as.matrix(sc_baseline[outcomes[i]]), na.rm=T),3)
  df_ancova[2,1,i] <- paste(paste("(",round(sd(as.matrix(sc_baseline[outcomes[i]]), na.rm=T),3),sep=""),")",sep="")
  
  ols <- lm(as.formula(paste(outcomes[i],"information*deliberation+region",sep="~")), data=sc_baseline[sc_baseline$district_baraza == 0,])
  vcov_cluster <- vcovCR(ols, cluster = sc_baseline$clusterID[sc_baseline$district_baraza == 0], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  
  #df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,5],res[2,6],conf[2,4],conf[2,5],nobs(ols))
  #df_ols[,2,i] <- c(res[2,1],res[2,2],res[2,5],conf[2,4],conf[2,5],nobs(ols))
  #df_ancova[,2,i] <- c(res[2,1],res[2,5],nobs(ols))
  df_ancova[,4,i] <- c(round(res[2,1],3),paste(paste("(",round(res[2,2],3), sep=""),")",sep=""))
  df_ancova[1,5,i] <- ifelse(res[2,5]<.01,"**",ifelse(res[2,5]<.05,"*",ifelse(res[2,5]<.1,"+","")))
  df_ancova[2,5,i] <- nobs(ols)
  #df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,5],res[3,6],conf[3,4],conf[3,5],nobs(ols))
  #df_ols[,3,i] <- c(res[3,1],res[3,2],res[3,5],conf[3,4],conf[3,5],nobs(ols))
  #df_ols[,3,i] <- c(res[3,1],res[3,5],nobs(ols))
  df_ancova[,6,i] <- c(round(res[3,1],3),paste(paste("(",round(res[3,2],3), sep=""),")",sep=""))
  df_ancova[1,7,i] <- ifelse(res[3,5]<.01,"**",ifelse(res[3,5]<.05,"*",ifelse(res[3,5]<.1,"+","")))
  df_ancova[2,7,i] <- nobs(ols)
  ols <- lm(as.formula(paste(outcomes[i],"information:deliberation+region",sep="~")), data=sc_baseline[sc_baseline$district_baraza == 0 & (sc_baseline$information == sc_baseline$deliberation),])
  vcov_cluster <- vcovCR(ols, cluster = sc_baseline$clusterID[sc_baseline$district_baraza == 0 & (sc_baseline$information == sc_baseline$deliberation)], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  if (RI_conf_switch) {
    conf[6,4:5] <- RI_store$conf_1
    res[6,5] <- RI_store$pval_1
  }
  #df_ols[,1,i] <- c(res[6,1],res[6,2],res[6,5],res[6,6],conf[6,4],conf[6,5],nobs(ols))
  #df_ols[,1,i] <- c(res[6,1],res[6,2],res[6,5],conf[6,4],conf[6,5],nobs(ols))
  df_ancova[,2,i] <- c(round(res[6,1],3),paste(paste("(",round(res[6,2],3), sep=""),")",sep=""))
  df_ancova[1,3,i] <- ifelse(res[6,5]<.01,"**",ifelse(res[6,5]<.05,"*",ifelse(res[6,5]<.1,"+","")))
  df_ancova[2,3,i] <- nobs(ols)
  
  
  
  #district vs. sc
  ols <- lm(as.formula(paste(outcomes[i],"district_baraza+region",sep="~")), data=sc_baseline[(sc_baseline$information == 1 & sc_baseline$deliberation==1) | sc_baseline$district_baraza == 1 ,])
  vcov_cluster <- vcovCR(ols, cluster = sc_baseline$clusterID2[(sc_baseline$information == 1 & sc_baseline$deliberation==1) | sc_baseline$district_baraza == 1 ], type = "CR0")
  res <- coef_test(ols, vcov_cluster)
  conf <- conf_int(ols, vcov_cluster)
  if (RI_conf_switch) {
    RI_store <- RI_conf_dist(i,outcomes, subset(sc_baseline, ((information == 1 & deliberation==1) | district_baraza == 1)) , ctrls = "region", nr_repl = glob_repli, sig = glob_sig)
    conf[2,4:5] <-  RI_store$conf
    res[2,5] <- RI_store$pval
  }
  #df_ols[,4,i] <- c(res[6,4],res[6,2],res[6,5],res[6,6],conf[6,4],conf[6,5],nobs(ols))
  #df_ols[,4,i] <- c(res[6,4],res[6,2],res[6,5],conf[6,4],conf[6,5],nobs(ols))
  #df_ols[,4,i] <- c(res[2,1],res[2,5],nobs(ols))
  df_ancova[,8,i] <- c(round(res[2,1],3),paste(paste("(",round(res[2,2],3), sep=""),")",sep=""))
  df_ancova[1,9,i] <- ifelse(res[2,5]<.01,"**",ifelse(res[2,5]<.05,"*",ifelse(res[2,5]<.1,"+","")))
  df_ancova[2,9,i] <- nobs(ols)
  
}
table_maker(df_ancova, c(11), "C:/Users/u0127963/Desktop/PhD/baraza/sc_report/table1.csv")
