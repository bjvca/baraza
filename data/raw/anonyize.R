rm(list=ls())
library(readstata13)
set.seed(12345) ### for gps offset
### anonymize endline data and merge in treatment assignment
#dta <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/baraza_enline_data_dummy.csv")
dta <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/clean.csv")


dta <- dta[,1:274]
dta <- dta[c(10,11,12,13,18,26:length(names(dta)))]

dta$baraza.labour.1..D1.1 <- NULL
dta$baraza.labour.2..D1.1 <- NULL
dta$baraza.labour.3..D1.1 <- NULL
dta$baraza.labour.4..D1.1 <- NULL
dta$baraza.labour.5..D1.1 <- NULL
dta$baraza.labour.6..D1.1 <- NULL
dta$baraza.labour.7..D1.1 <- NULL
dta$baraza.labour.8..D1.1 <- NULL
dta$baraza.labour.9..D1.1 <- NULL
dta$baraza.labour.10..D1.1 <- NULL
dta$baraza.labour.11..D1.1 <- NULL
dta$baraza.labour.12..D1.1 <- NULL
dta$baraza.labour.13..D1.1 <- NULL
dta$baraza.labour.14..D1.1 <- NULL
dta$baraza.labour.15..D1.1 <- NULL

dta <- dta[!duplicated(dta$hh_id),] #this should not be needed if data is cleaned
### merge in treatment status

treats <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/questionnaire/final_list_5.csv")
endline <- merge(treats, dta, by.x=c("district","subcounty"), by.y=c("district","sub"))
## rename hh_id to hhid
names(endline)[names(endline) == 'hh_id'] <- 'hhid'
 
### save as a publically accessible endline dataset
write.csv(endline,"/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/endline.csv",row.names=FALSE)

###### also prepare baseline data
library(readstata13)
hh_level <- read.dta13("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/baseline/quant/Household data/BARAZA_cleaned datasets/BARAZA HOUSEHOLD MAIN DATA FILE.dta")
get_health <- read.dta13("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/baseline/quant/Household data/Parent file_hhd data.dta")[,c("key","feverd21_fever","delivery_birthd21_delivery_birth")]
hh_level <- merge(hh_level, get_health,by="key")
#merge in indicator of demo visit
hh_level <- merge(hh_level, read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/baseline/quant/Household data/working/baseline_hh_char.csv")[c("hhsize", "femhead", "agehead", "farmsize", "hhid", "treat")], by="hhid")
hh_level$has_phone <-  hh_level$a31!= 0

inputs <- read.dta13("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/baseline/quant/Household data/BARAZA_cleaned datasets/BARAZA INPUTS DATA FILE.dta")
seed <- subset(inputs, inputs == "Improved crop seed" & b121=="Yes")
seed$b122 <- seed$b122 %in% c("Extension worker","Commmunity based facilitator")
hh_level <- merge(hh_level,seed[c("key","b121")], by="key", all.x=T)
hh_level <- merge(hh_level,seed[c("key","b122")], by="key", all.x=T)
names(hh_level)[names(hh_level) == 'b121'] <- 'used_seed'
names(hh_level)[names(hh_level) == 'b122'] <- 'seed_NARO'

fertilizer <- subset(inputs, inputs == "Fertilizers (DAP, Urea, NPK, etc.)" & b121=="Yes")
hh_level <- merge(hh_level,fertilizer[c("key","b121")], by="key", all.x=T) 
names(hh_level)[names(hh_level) == 'b121'] <- 'used_fert'

livestock_tech <- subset(inputs, inputs %in% c("Improved livestock breeds","Artificial Insemination","Veterinary drugs","Livestock feed") & b121=="Yes")
hh_level <- merge(hh_level,livestock_tech[c("key","b121")], by="key", all.x=T) 
names(hh_level)[names(hh_level) == 'b121'] <- 'used_livestock_tech'

chemicals <- subset(inputs, inputs %in% c("Pesticide","Herbicide") & b121=="Yes")
hh_level <- merge(hh_level,chemicals[c("key","b121")], by="key", all.x=T) 
names(hh_level)[names(hh_level) == 'b121'] <- 'used_chem'
hh_level$used_fert[is.na(hh_level$used_fert)] <- "No"
hh_level$used_seed[is.na(hh_level$used_seed)] <- "No"
hh_level$seed_NARO[is.na(hh_level$seed_NARO)] <- FALSE
hh_level$used_chem[is.na(hh_level$used_chem)] <- "No"
hh_level$used_livestock_tech[is.na(hh_level$used_livestock_tech)] <- "No"

get_sickdays <- read.dta13("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/baseline/quant/Household data/BARAZA_cleaned datasets/BARAZA HOUSEHOLD IMPACT DATA FILE.dta")
sickdays_hh <- aggregate(cbind(get_sickdays[c("d111","d112","d113")]), by=list(get_sickdays$key), FUN= sum, na.rm=T)
names(sickdays_hh) <- c("key","tot_sick","not_work","not_school")
hh_level <- merge(hh_level,sickdays_hh, by="key", all.x=T)

get_participation <- read.dta13("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/baseline/quant/Household data/BARAZA_cleaned datasets/BARAZA ELECTION PARTICIPATION DATA FILE.dta")

hh_level <-  merge(hh_level,reshape(get_participation[,1:3], idvar = "key", timevar = "sno", direction = "wide"), by="key", all.x=T)

get_contributions <- read.dta13("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/baseline/quant/Household data/BARAZA_cleaned datasets/BARAZA PROJECTS CONTRIBUTIONS DATA FILE.dta")

contributions <- reshape(get_contributions[,c(1:2,4)], idvar = "key", timevar = "project", direction = "wide")

contributions$cschoolk <- contributions$g21.School %in% c("In-kind","Both 1 & 2")
contributions$cschoolk[is.na(contributions$cschoolk)] <- NA
contributions$cschoolc <- contributions$g21.School %in% c("Money","Both 1 & 2")
contributions$cschoolc[is.na(contributions$cschoolc)] <- NA

contributions$chealthk <- contributions$"g21.Health centre" %in% c("In-kind","Both 1 & 2")
contributions$chealthk[is.na(contributions$chealthk)] <- NA
contributions$chealthc <- contributions$"g21.Health centre" %in% c("Money","Both 1 & 2")
contributions$chealthc[is.na(contributions$chealthc)] <- NA

contributions$croadk <- contributions$"g21.Road/bridge" %in% c("In-kind","Both 1 & 2")
contributions$croadk[is.na(contributions$croadk)] <- NA
contributions$croadc <- contributions$"g21.Road/bridge" %in% c("Money","Both 1 & 2")
contributions$croadc[is.na(contributions$croadc)] <- NA

contributions$cwaterk <- contributions$"g21.Drinking water facility" %in% c("In-kind","Both 1 & 2")
contributions$cwaterk[is.na(contributions$cwaterk)] <- NA
contributions$cwaterc <- contributions$"g21.Drinking water facility" %in% c("Money","Both 1 & 2")
contributions$cwaterc[is.na(contributions$cwaterc)] <- NA

contributions$cdamk <- contributions$"g21.Dam/irrigation facility" %in% c("In-kind","Both 1 & 2")
contributions$cdamk[is.na(contributions$cdamk)] <- NA
contributions$cdamc <- contributions$"g21.Dam/irrigation facility" %in% c("Money","Both 1 & 2")
contributions$cdamc[is.na(contributions$cdamc)] <- NA


contributions$cbuildk <- contributions$"g21.Other building/structure" %in% c("In-kind","Both 1 & 2")
contributions$cbuildk[is.na(contributions$cbuildk)] <- NA
contributions$cbuildc <- contributions$"g21.Other building/structure" %in% c("Money","Both 1 & 2")
contributions$cbuildc[is.na(contributions$cbuildc)] <- NA
contributions <- contributions[,c(1,9:20)]


hh_level <-  merge(hh_level,contributions, by="key", all.x=T)

####
hh_level <- hh_level[c(2,12:15,18,19, 23:718, 724:758)]

## we need an offset for the gps coordinates
hh_level$a26a <- hh_level$a26a +rnorm(dim(hh_level)[1],0,.05)
#for the longitudes (b)
hh_level$a26b <- hh_level$a26b + rnorm(dim(hh_level)[1],0,.025)
hh_level <- hh_level[!duplicated(hh_level$hhid),]


write.csv(hh_level,"/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/baseline.csv",row.names=FALSE)

###anonymize sc level data


sc_dta <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/sc_level.csv")
sc_dta <- sc_dta[c(9:13,15,18:229)]
#create a sc_id - just used row number as ID
sc_dta$scID <- as.numeric(rownames(sc_dta))
write.csv(sc_dta,"/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/sc_level_endline.csv",row.names=FALSE)

### also baseline sc level data

sc_base <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/baseline/quant/Subcounty data/main.csv")

sc_base <- sc_base[c(7:10,15,17:20,22:575)]

write.csv(sc_base,"/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/sc_level_baseline.csv",row.names=FALSE)

#exit

