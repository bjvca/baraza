rm(list=ls())
library(readstata13)
set.seed(12345) ### for gps offset
### anonymize endline data and merge in treatment assignment
#dta <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/baraza_enline_data_dummy.csv")
dta <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/clean.csv")


dta <- dta[,1:256]
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
#dta$baraza.labour.13..D1.1 <- NULL
#dta$baraza.labour.14..D1.1 <- NULL
#dta$baraza.labour.15..D1.1 <- NULL

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
hh_level <- merge(hh_level,seed[c("key","b121")], by="key", all.x=T)
names(hh_level)[names(hh_level) == 'b121'] <- 'used_seed'

fertilizer <- subset(inputs, inputs == "Fertilizers (DAP, Urea, NPK, etc.)" & b121=="Yes")
hh_level <- merge(hh_level,fertilizer[c("key","b121")], by="key", all.x=T) 
names(hh_level)[names(hh_level) == 'b121'] <- 'used_fert'
hh_level$used_fert[is.na(hh_level$used_fert)] <- "No"
hh_level$used_seed[is.na(hh_level$used_seed)] <- "No"

get_sickdays <- read.dta13("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/baseline/quant/Household data/BARAZA_cleaned datasets/BARAZA HOUSEHOLD IMPACT DATA FILE.dta")
sickdays_hh <- aggregate(cbind(get_sickdays[c("d111","d112","d113")]), by=list(get_sickdays$key), FUN= sum, na.rm=T)
names(sickdays_hh) <- c("key","tot_sick","not_work","not_school")
hh_level <- merge(hh_level,sickdays_hh, by="key", all.x=T)

hh_level <- hh_level[c(2,12:15,18,19, 23:718, 724:736)]

## we need an offset for the gps coordinates
hh_level$a26a <- hh_level$a26a +rnorm(dim(hh_level)[1],0,.05)
#for the longitudes (b)
hh_level$a26b <- hh_level$a26b + rnorm(dim(hh_level)[1],0,.025)
hh_level <- hh_level[!duplicated(hh_level$hhid),]


write.csv(hh_level,"/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/public/baseline.csv",row.names=FALSE)

#exit

