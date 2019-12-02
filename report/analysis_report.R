library(ggplot2)
library(MatchIt) 
set.seed(12345)

endline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/endline_dummy.csv")[c("hhid","information","deliberation","ka21_region")]

###randomly generate adoption 
endline$access_extension <- rbinom(n=dim(endline)[1], size=1, prob=0.1) 

##simple treatment-control difference

ols <- lm(access_extension~information*deliberation+ka21_region, data=endline) 
df <-data.frame( rbind(summary(ols)$coef[2,1],summary(ols)$coef[2,2],summary(ols)$coef[2,4]), rbind(summary(ols)$coef[3,1],summary(ols)$coef[3,2],summary(ols)$coef[3,4]), rbind(summary(ols)$coef[7,1],summary(ols)$coef[7,2],summary(ols)$coef[7,4]))


##ancova
## merge in baseline
baseline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/endline_dummy.csv")[c("hhid","b21")]

dta <- merge(endline, baseline, by="hhid")
ols <- lm(access_extension~information*deliberation+ka21_region+b21 , data=dta) 
df <-data.frame( rbind(summary(ols)$coef[2,1],summary(ols)$coef[2,2],summary(ols)$coef[2,4]), rbind(summary(ols)$coef[3,1],summary(ols)$coef[3,2],summary(ols)$coef[3,4]), rbind(summary(ols)$coef[8,1],summary(ols)$coef[8,2],summary(ols)$coef[8,4]))


## dif-in-dif
baseline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/endline_dummy.csv")[c("hhid","b21","information","deliberation","ka21_region")]
baseline$time <- 0
names(baseline)[names(baseline) == "b21"] <- "access_extension"
endline$time <- 1
endline <- endline[names(baseline)]
baseline$access_extension <- as.numeric((baseline$access_extension == "Yes"))

dta <- rbind(endline, baseline)

ols <- lm(access_extension~information*deliberation*time+ka21_region, data=dta)
df <-data.frame( rbind(summary(ols)$coef[9,1],summary(ols)$coef[9,2],summary(ols)$coef[9,4]), rbind(summary(ols)$coef[10,1],summary(ols)$coef[10,2],summary(ols)$coef[10,4]), rbind(summary(ols)$coef[11,1],summary(ols)$coef[11,2],summary(ols)$coef[11,4]))

###matched dif-in-dif
baseline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/endline_dummy.csv")[c("hhid","information","deliberation","ka21_region","b21")]
###merge in baseline hhchar (these were constructed using initial.do at the time of the preparation of baseline balance table)
hh_char <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/data/baseline_hh_char.csv") 
baseline <- merge(baseline,read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/report/data/baseline_hh_char.csv") , by="hhid")
baseline$farmsize <- log(baseline$farmsize)
baseline$ironroof <- as.numeric(baseline$a512 =="Corrugated iron sheets")
baseline$improved_wall <- as.numeric(baseline$a513 %in% c("Mud_bricks_burnt_bricks","Concrete_blocks") )
baseline$farmsize[is.infinite(baseline$farmsize)] <- NA
baseline$has_phone <-  baseline$a31!= 0
baseline$head_sec <- as.numeric(baseline$a36) > 15
baseline <- baseline[complete.cases(baseline[c("hhsize","femhead","agehead","farmsize","ironroof","improved_wall","has_phone","head_sec","a26a","a26b")]),] 
nearest.match <- matchit(formula = information ~ hhsize + femhead + agehead + farmsize + ironroof + improved_wall + has_phone +head_sec+a26a+a26b,  data =baseline,method = "nearest",distance = "logit")
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)

### this is the matched data -  we now need to extract these ids from the endline and stack base and endline to do a dif-in-dif

matched.endline <- endline[endline$hhid %in% matched.baseline$hhid,]
matched.baseline$time <- 0
matched.endline$time <- 1
names(matched.baseline)[names(matched.baseline) == "b21"] <- "access_extension"
matched.baseline <- matched.baseline[,names(matched.endline)]

matched.baseline$access_extension <- as.numeric((matched.baseline$access_extension == "Yes"))

matched.dta <- rbind(matched.endline, matched.baseline[,names(matched.endline)])

ols <- lm(access_extension~information*deliberation*time+ka21_region, data=matched.dta)
df <-data.frame( rbind(summary(ols)$coef[9,1],summary(ols)$coef[9,2],summary(ols)$coef[9,4]))

### now for delib
nearest.match <- matchit(formula = deliberation ~ hhsize + femhead + agehead + farmsize + ironroof + improved_wall + has_phone +head_sec+a26a+a26b,  data =baseline,method = "nearest",distance = "logit")
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)

### this is the matched data -  we now need to extract these ids from the endline and stack base and endline to do a dif-in-dif

matched.endline <- endline[endline$hhid %in% matched.baseline$hhid,]
matched.baseline$time <- 0
matched.endline$time <- 1
names(matched.baseline)[names(matched.baseline) == "b21"] <- "access_extension"
matched.baseline <- matched.baseline[,names(matched.endline)]

matched.baseline$access_extension <- as.numeric((matched.baseline$access_extension == "Yes"))

matched.dta <- rbind(matched.endline, matched.baseline[,names(matched.endline)])

ols <- lm(access_extension~information*deliberation*time+ka21_region, data=matched.dta)
df <-data.frame(df, rbind(summary(ols)$coef[10,1],summary(ols)$coef[10,2],summary(ols)$coef[10,4]))

### now for delib
nearest.match <- matchit(formula = deliberation*information ~ hhsize + femhead + agehead + farmsize + ironroof + improved_wall + has_phone +head_sec+a26a+a26b,  data =baseline,method = "nearest",distance = "logit")
summary(nearest.match)
#plot(nearest.match)
#plot(nearest.match, type="hist")
#plot(nearest.match, type="jitter")

matched.baseline <- match.data(nearest.match)

### this is the matched data -  we now need to extract these ids from the endline and stack base and endline to do a dif-in-dif

matched.endline <- endline[endline$hhid %in% matched.baseline$hhid,]
matched.baseline$time <- 0
matched.endline$time <- 1
names(matched.baseline)[names(matched.baseline) == "b21"] <- "access_extension"
matched.baseline <- matched.baseline[,names(matched.endline)]

matched.baseline$access_extension <- as.numeric((matched.baseline$access_extension == "Yes"))

matched.dta <- rbind(matched.endline, matched.baseline[,names(matched.endline)])

ols <- lm(access_extension~information*deliberation*time+ka21_region, data=matched.dta)
df <-data.frame(df, rbind(summary(ols)$coef[11,1],summary(ols)$coef[11,2],summary(ols)$coef[11,4]))


