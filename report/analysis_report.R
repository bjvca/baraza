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
baseline <- read.csv("/home/bjvca/Dropbox (IFPRI)/baraza/Impact Evaluation Surveys/endline/data/raw/endline_dummy.csv")[c("hhid","b21","a6","b11","information","deliberation","ka21_region")]
nearest.match <- matchit(formula = information ~ a6,  data =baseline,method = "nearest",distance = "logit",discard="control")

data.matched <- match.data(nearest.match)
imbalance(group=data.matched$information, data=data.matched)

ols <- lm((b21=="Yes") ~ information*deliberation +ka21_region,data = data.matched)
df <-data.frame( rbind(summary(ols)$coef[2,1],summary(ols)$coef[2,2],summary(ols)$coef[2,4]), rbind(summary(ols)$coef[3,1],summary(ols)$coef[3,2],summary(ols)$coef[3,4]), rbind(summary(ols)$coef[7,1],summary(ols)$coef[7,2],summary(ols)$coef[7,4]))
