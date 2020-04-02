rm(list=ls())
#remove and rm can be used to remove objects
#list: a character vector naming objects to be removed
#ls shows what data sets and functions a user has defined

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

### this loads sub-county endline data
sc_endline <- read.csv(paste(path,"data/public/sc_level_endline.csv", sep ="/"))
### this loads the treatment assignment
treats <- read.csv(paste(path,"questionnaire/final_list_5.csv", sep ="/"))

### this merges in treatments
sc_endline <- merge(treats, sc_endline, by.x=c("district","subcounty"), by.y=c("district","sub"))

#drop sc that received the district treatment
sc_endline_no_dist <- subset(sc_endline, district_baraza == 0)

## make an indicator that pools I, D and ID - | means OR

sc_endline_no_dist$ind_treat <- (sc_endline_no_dist$information == 1 | sc_endline_no_dist$deliberation == 1)  

#tabulate 
table(sc_endline_no_dist$ind_treat)

#baraza.km.D1 is total km of tarmac roads in SC (see questionnaire)
#first look at this variable
sc_endline_no_dist$baraza.km.D1
# see those 999s? That is a missings code, you need to get rid of it because it will screw up your means (frequent beginer mistake for data analysts, the do not look at their data)
#other ways to detect outliers etc
hist(sc_endline_no_dist$baraza.km.D1) 
summary(sc_endline_no_dist$baraza.km.D1)
boxplot(sc_endline_no_dist$baraza.km.D1)

# in any case, we need to replace the 999s by NA (code for not available)
## this is how you assign something: if D1 == 999, assign NA
sc_endline_no_dist$baraza.km.D1[ sc_endline_no_dist$baraza.km.D1 == 999 ] <- NA
#look at data again, not means etc make sense

### now you can look at impact
### means by treatment group
tapply(sc_endline_no_dist$baraza.km.D1, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
### ughh, those that received a baraza seem to have on average only 3km, while those that did not receive a baraza had 6km
#is this difference significant? Lets do a t-test
t.test(sc_endline_no_dist$baraza.km.D1~sc_endline_no_dist$ind_treat, var.equal=T)
## run a regression
summary(lm(baraza.km.D1~ind_treat, data = sc_endline_no_dist))

###Caroline's first R day:
#####C1#####
sc_endline_no_dist$baraza.visit_times.C1
hist(sc_endline_no_dist$baraza.visit_times.C1)
summary(sc_endline_no_dist$baraza.visit_times.C1)
#is 48 an outlier?
boxplot(sc_endline_no_dist$baraza.visit_times.C1)
#means by treatment group
tapply(sc_endline_no_dist$baraza.visit_times.C1, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#na.rm=T is Excluding Missing Values from Analyses
#difference significant? t-test
t.test(sc_endline_no_dist$baraza.visit_times.C1~sc_endline_no_dist$ind_treat, var.equal=T)
#regression
summary(lm(baraza.visit_times.C1~ind_treat, data = sc_endline_no_dist))
#lm means linear model

#####C2#####
sc_endline_no_dist$baraza.visit_times.C2
#replace 999s by NA
sc_endline_no_dist$baraza.visit_times.C2[ sc_endline_no_dist$baraza.visit_times.C2 == 999 ] <- NA
hist(sc_endline_no_dist$baraza.visit_times.C2)
summary(sc_endline_no_dist$baraza.visit_times.C2)
boxplot(sc_endline_no_dist$baraza.visit_times.C2)
#means by treatment group
tapply(sc_endline_no_dist$baraza.visit_times.C2, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#na.rm=T is Excluding Missing Values from Analyses
#difference significant? t-test
t.test(sc_endline_no_dist$baraza.visit_times.C2~sc_endline_no_dist$ind_treat, var.equal=T)
#regression
summary(lm(baraza.visit_times.C2~ind_treat, data = sc_endline_no_dist))
#lm means linear model

###C3###
sc_endline_no_dist$baraza.visit_times.C3
#replace 999s by NA
sc_endline_no_dist$baraza.visit_times.C3[ sc_endline_no_dist$baraza.visit_times.C3 == 999 ] <- NA
#means by treatment group
tapply(sc_endline_no_dist$baraza.visit_times.C3, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#regression
summary(lm(baraza.visit_times.C3~ind_treat, data = sc_endline_no_dist))

###C4###
sc_endline_no_dist$baraza.visit_times.C4
#means by treatment group
tapply(sc_endline_no_dist$baraza.visit_times.C4, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#regression
summary(lm(baraza.visit_times.C4~ind_treat, data = sc_endline_no_dist))

###C5###
sc_endline_no_dist$baraza.visit_times.C5
#replace 999s by NA
sc_endline_no_dist$baraza.visit_times.C5[ sc_endline_no_dist$baraza.visit_times.C5 == 999 ] <- NA
#means by treatment group
tapply(sc_endline_no_dist$baraza.visit_times.C5, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#regression
summary(lm(baraza.visit_times.C5~ind_treat, data = sc_endline_no_dist))

###C6###
sc_endline_no_dist$baraza.visit_times.C6
#replace 999s by NA
sc_endline_no_dist$baraza.visit_times.C6[ sc_endline_no_dist$baraza.visit_times.C6 == 999 ] <- NA
#means by treatment group
tapply(sc_endline_no_dist$baraza.visit_times.C6, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#regression
summary(lm(baraza.visit_times.C6~ind_treat, data = sc_endline_no_dist))

###C7###
sc_endline_no_dist$baraza.visit_times.C7
#replace 999s by NA
sc_endline_no_dist$baraza.visit_times.C7[ sc_endline_no_dist$baraza.visit_times.C7 == 999 ] <- NA
#means by treatment group
tapply(sc_endline_no_dist$baraza.visit_times.C7, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#regression
summary(lm(baraza.visit_times.C7~ind_treat, data = sc_endline_no_dist))

###C8###
sc_endline_no_dist$baraza.visit_times.C8
#means by treatment group
tapply(sc_endline_no_dist$baraza.visit_times.C8, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#regression
summary(lm(baraza.visit_times.C8~ind_treat, data = sc_endline_no_dist))

###C9###
sc_endline_no_dist$baraza.visit_times.C9
#means by treatment group
tapply(sc_endline_no_dist$baraza.visit_times.C9, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#regression
summary(lm(baraza.visit_times.C9~ind_treat, data = sc_endline_no_dist))

###D1###
sc_endline_no_dist$baraza.km.D1
#replace 999s by NA
sc_endline_no_dist$baraza.km.D1[ sc_endline_no_dist$baraza.km.D1 == 999 ] <- NA
#means by treatment group
tapply(sc_endline_no_dist$baraza.km.D1, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
#regression
summary(lm(baraza.km.D1~ind_treat, data = sc_endline_no_dist))

###D2-4b###
sc_endline_no_dist$baraza.km.D2
sc_endline_no_dist$baraza.km.D2[ sc_endline_no_dist$baraza.km.D2 == 999 ] <- NA
sc_endline_no_dist$baraza.km.D3
sc_endline_no_dist$baraza.km.D3[ sc_endline_no_dist$baraza.km.D3 == 999 ] <- NA
sc_endline_no_dist$baraza.km.D4a
sc_endline_no_dist$baraza.km.D4a[ sc_endline_no_dist$baraza.km.D4a == 999 ] <- NA
sc_endline_no_dist$baraza.km.D4b
sc_endline_no_dist$baraza.km.D4b[ sc_endline_no_dist$baraza.km.D4b == 999 ] <- NA
#regression
summary(lm(baraza.km.D2~ind_treat, data = sc_endline_no_dist))
summary(lm(baraza.km.D3~ind_treat, data = sc_endline_no_dist))
summary(lm(baraza.km.D4a~ind_treat, data = sc_endline_no_dist))
tapply(sc_endline_no_dist$baraza.km.D4a, sc_endline_no_dist$ind_treat, FUN=mean, na.rm=T)
t.test(sc_endline_no_dist$baraza.km.D4a~sc_endline_no_dist$ind_treat, var.equal=T)

###E1a-e7###
sc_endline_no_dist$baraza.health.E2a
summary(lm(baraza.health.E2a~ind_treat, data = sc_endline_no_dist))
sc_endline_no_dist$baraza.health.E2a
summary(lm(baraza.health.E2a~ind_treat, data = sc_endline_no_dist))
sc_endline_no_dist$baraza.works.E4a
summary(lm(baraza.works.E4a~ind_treat, data = sc_endline_no_dist))
sc_endline_no_dist$baraza.works.E4a
summary(lm(baraza.works.E4a~ind_treat, data = sc_endline_no_dist))
sc_endline_no_dist$baraza.E7
summary(lm(baraza.E7~ind_treat, data = sc_endline_no_dist))
#how to measure difference date and today
sc_endline_no_dist$baraza.finance1.E5a
summary(lm(baraza.finance1.E5a~ind_treat, data = sc_endline_no_dist))

########31/03/20
#SUBSECTION H1 - HEALTH
sc_endline_no_dist$baraza.H1[sc_endline_no_dist$baraza.H1==999] <- NA
sc_endline_no_dist$baraza.H1 <- as.numeric(as.character(sc_endline_no_dist$baraza.H1))
plot(density(sc_endline_no_dist$baraza.H1, na.rm=T))
ols <- lm(baraza.H1~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H2[sc_endline_no_dist$baraza.H2==999] <- NA
sc_endline_no_dist$baraza.H2 <- as.numeric(as.character(sc_endline_no_dist$baraza.H2))
plot(density(sc_endline_no_dist$baraza.H2, na.rm=T))
ols <- lm(baraza.H2~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H5[sc_endline_no_dist$baraza.H5==999] <- NA
sc_endline_no_dist$baraza.H5 <- as.numeric(as.character(sc_endline_no_dist$baraza.H5))
plot(density(sc_endline_no_dist$baraza.H5, na.rm=T))
ols <- lm(baraza.H5~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H6[sc_endline_no_dist$baraza.H6==999] <- NA
sc_endline_no_dist$baraza.H6 <- as.numeric(as.character(sc_endline_no_dist$baraza.H6))
plot(density(sc_endline_no_dist$baraza.H6, na.rm=T))
ols <- lm(baraza.H6~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H7[sc_endline_no_dist$baraza.H7==999] <- NA
sc_endline_no_dist$baraza.H7 <- as.numeric(as.character(sc_endline_no_dist$baraza.H7))
plot(density(sc_endline_no_dist$baraza.H7, na.rm=T))
ols <- lm(baraza.H7~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#does NOT look like the gauss curve
sc_endline_no_dist$baraza.H10[sc_endline_no_dist$baraza.H10==999] <- NA
sc_endline_no_dist$baraza.H10 <- as.numeric(as.character(sc_endline_no_dist$baraza.H10))
plot(density(sc_endline_no_dist$baraza.H10, na.rm=T))
ols <- lm(baraza.H10~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#outlier?
sc_endline_no_dist$baraza.H24[sc_endline_no_dist$baraza.H24==999] <- NA
sc_endline_no_dist$baraza.H24 <- as.numeric(as.character(sc_endline_no_dist$baraza.H24))
plot(density(sc_endline_no_dist$baraza.H24, na.rm=T))
ols <- lm(baraza.H24~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#outlier?
sc_endline_no_dist$baraza.H25[sc_endline_no_dist$baraza.H25==999] <- NA
sc_endline_no_dist$baraza.H25 <- as.numeric(as.character(sc_endline_no_dist$baraza.H25))
plot(density(sc_endline_no_dist$baraza.H25, na.rm=T))
ols <- lm(baraza.H25~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#does NOT look like the gauss curve
sc_endline_no_dist$baraza.H30[sc_endline_no_dist$baraza.H30==999] <- NA
sc_endline_no_dist$baraza.H30 <- as.numeric(as.character(sc_endline_no_dist$baraza.H30))
plot(density(sc_endline_no_dist$baraza.H30, na.rm=T))
ols <- lm(baraza.H30~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#does NOT look like the gauss curve
sc_endline_no_dist$baraza.H33[sc_endline_no_dist$baraza.H33==999] <- NA
sc_endline_no_dist$baraza.H33 <- as.numeric(as.character(sc_endline_no_dist$baraza.H33))
plot(density(sc_endline_no_dist$baraza.H33, na.rm=T))
ols <- lm(baraza.H33~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#does NOT look like the gauss curve
sc_endline_no_dist$baraza.H37[sc_endline_no_dist$baraza.H37==999] <- NA
sc_endline_no_dist$baraza.H37 <- as.numeric(as.character(sc_endline_no_dist$baraza.H37))
plot(density(sc_endline_no_dist$baraza.H37, na.rm=T))
ols <- lm(baraza.H37~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#outlier?
sc_endline_no_dist$baraza.H43[sc_endline_no_dist$baraza.H43==999] <- NA
sc_endline_no_dist$baraza.H43 <- as.numeric(as.character(sc_endline_no_dist$baraza.H43))
plot(density(sc_endline_no_dist$baraza.H43, na.rm=T))
ols <- lm(baraza.H43~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#does NOT look like the gauss curve
sc_endline_no_dist$baraza.H44[sc_endline_no_dist$baraza.H44==999] <- NA
sc_endline_no_dist$baraza.H44 <- as.numeric(as.character(sc_endline_no_dist$baraza.H44))
plot(density(sc_endline_no_dist$baraza.H44, na.rm=T))
ols <- lm(baraza.H44~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H53[sc_endline_no_dist$baraza.H53==999] <- NA
sc_endline_no_dist$baraza.H53 <- as.numeric(as.character(sc_endline_no_dist$baraza.H53))
plot(density(sc_endline_no_dist$baraza.H53, na.rm=T))
ols <- lm(baraza.H53~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H54[sc_endline_no_dist$baraza.H54==999] <- NA
sc_endline_no_dist$baraza.H54 <- as.numeric(as.character(sc_endline_no_dist$baraza.H54))
plot(density(sc_endline_no_dist$baraza.H54, na.rm=T))
ols <- lm(baraza.H54~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H55[sc_endline_no_dist$baraza.H55==999] <- NA
sc_endline_no_dist$baraza.H55 <- as.numeric(as.character(sc_endline_no_dist$baraza.H55))
plot(density(sc_endline_no_dist$baraza.H55, na.rm=T))
ols <- lm(baraza.H55~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H57[sc_endline_no_dist$baraza.H57==999] <- NA
sc_endline_no_dist$baraza.H57 <- as.numeric(as.character(sc_endline_no_dist$baraza.H57))
plot(density(sc_endline_no_dist$baraza.H57, na.rm=T))
ols <- lm(baraza.H57~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H58[sc_endline_no_dist$baraza.H58==999] <- NA
sc_endline_no_dist$baraza.H58 <- as.numeric(as.character(sc_endline_no_dist$baraza.H58))
plot(density(sc_endline_no_dist$baraza.H58, na.rm=T))
ols <- lm(baraza.H58~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H62[sc_endline_no_dist$baraza.H62==999] <- NA
sc_endline_no_dist$baraza.H62 <- as.numeric(as.character(sc_endline_no_dist$baraza.H62))
plot(density(sc_endline_no_dist$baraza.H62, na.rm=T))
ols <- lm(baraza.H62~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H62[sc_endline_no_dist$baraza.H62==999] <- NA
sc_endline_no_dist$baraza.H62 <- as.numeric(as.character(sc_endline_no_dist$baraza.H62))
plot(density(sc_endline_no_dist$baraza.H62, na.rm=T))
ols <- lm(baraza.H62~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#yes/no question
#sc_endline_no_dist$baraza.H65[sc_endline_no_dist$baraza.H65==98] <- NA
#sc_endline_no_dist$baraza.H65[sc_endline_no_dist$baraza.H65==999] <- NA
#sc_endline_no_dist$baraza.H65 <- as.numeric(as.character(sc_endline_no_dist$baraza.H65))
#plot(density(sc_endline_no_dist$baraza.H65, na.rm=T))
#ols <- lm(baraza.H65~information*deliberation+region, data = sc_endline_no_dist)
#vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
#res <- coef_test(ols, vcov_cluster)
#print(res)

#OUTLIER
sc_endline_no_dist$baraza.H70[sc_endline_no_dist$baraza.H70==999] <- NA
sc_endline_no_dist$baraza.H70 <- as.numeric(as.character(sc_endline_no_dist$baraza.H70))
plot(density(sc_endline_no_dist$baraza.H70, na.rm=T))
ols <- lm(baraza.H70~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#SUBSECTION : WATER INFRASTRUCTURE
sc_endline_no_dist$baraza.H72[sc_endline_no_dist$baraza.H72==999] <- NA
sc_endline_no_dist$baraza.H72 <- as.numeric(as.character(sc_endline_no_dist$baraza.H72))
plot(density(sc_endline_no_dist$baraza.H72, na.rm=T))
ols <- lm(baraza.H72~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#SUBSECTION H3: EDUCATION
sc_endline_no_dist$baraza.H80[sc_endline_no_dist$baraza.H80==999] <- NA
sc_endline_no_dist$baraza.H80 <- as.numeric(as.character(sc_endline_no_dist$baraza.H80))
plot(density(sc_endline_no_dist$baraza.H80, na.rm=T))
ols <- lm(baraza.H80~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H81[sc_endline_no_dist$baraza.H81==999] <- NA
sc_endline_no_dist$baraza.H81 <- as.numeric(as.character(sc_endline_no_dist$baraza.H81))
plot(density(sc_endline_no_dist$baraza.H81, na.rm=T))
ols <- lm(baraza.H81~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#outlier?
sc_endline_no_dist$baraza.H109[sc_endline_no_dist$baraza.H109==999] <- NA
sc_endline_no_dist$baraza.H109 <- as.numeric(as.character(sc_endline_no_dist$baraza.H109))
plot(density(sc_endline_no_dist$baraza.H109, na.rm=T))
ols <- lm(baraza.H109~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

#yes/no question
#sc_endline_no_dist$baraza.H109[sc_endline_no_dist$baraza.H109==999] <- NA
#sc_endline_no_dist$baraza.H109 <- as.numeric(as.character(sc_endline_no_dist$baraza.H109))
#plot(density(sc_endline_no_dist$baraza.H109, na.rm=T))
#ols <- lm(baraza.H109~information*deliberation+region, data = sc_endline_no_dist)
#vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
#res <- coef_test(ols, vcov_cluster)
#print(res)

#yes/no question
#sc_endline_no_dist$baraza.H110[sc_endline_no_dist$baraza.H110==999] <- NA
#sc_endline_no_dist$baraza.H110 <- as.numeric(as.character(sc_endline_no_dist$baraza.H110))
#plot(density(sc_endline_no_dist$baraza.H110, na.rm=T))
#ols <- lm(baraza.H110~information*deliberation+region, data = sc_endline_no_dist)
#vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
#res <- coef_test(ols, vcov_cluster)
#print(res)

#####01/04/20
#####SECTION D: SUBCOUNTY'S BASIC INFORMATION
sc_endline_no_dist$baraza.km.D2[sc_endline_no_dist$baraza.km.D2==999] <- NA
sc_endline_no_dist$baraza.km.D2
summary(sc_endline_no_dist$baraza.km.D2)
plot(density(sc_endline_no_dist$baraza.km.D2, na.rm=T))
ols <- lm(baraza.km.D2~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
sum(is.na (sc_endline_no_dist$baraza.km.D2))

sc_endline_no_dist$baraza.km.D4b[sc_endline_no_dist$baraza.km.D4b==999] <- NA
sc_endline_no_dist$baraza.km.D4b
summary(sc_endline_no_dist$baraza.km.D4b)
plot(density(sc_endline_no_dist$baraza.km.D4b, na.rm=T))
ols <- lm(baraza.km.D4b~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
sum(is.na (sc_endline_no_dist$baraza.km.D4b))

#####SECTION E: GOVERNMENT AND COMMUNITY BODIES IN THE SUBCOUNTY
#/

#####SECTION F:  COMMUNITY MEETINGS HELD IN THE SUBCOUNTY
#/

#####SUBSECTION H1 - HEALTH
sc_endline_no_dist$baraza.H1[sc_endline_no_dist$baraza.H1==999] <- NA
sc_endline_no_dist$baraza.H1 <- as.numeric(as.character(sc_endline_no_dist$baraza.H1))
plot(density(sc_endline_no_dist$baraza.H1, na.rm=T))
ols <- lm(baraza.H1~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H2[sc_endline_no_dist$baraza.H2==999] <- NA
sc_endline_no_dist$baraza.H2 <- as.numeric(as.character(sc_endline_no_dist$baraza.H2))
plot(density(sc_endline_no_dist$baraza.H2, na.rm=T))
ols <- lm(baraza.H2~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H5[sc_endline_no_dist$baraza.H5==999] <- NA
sc_endline_no_dist$baraza.H5 <- as.numeric(as.character(sc_endline_no_dist$baraza.H5))
sc_endline_no_dist$baraza.H5
summary(sc_endline_no_dist$baraza.H5)
table(sc_endline_no_dist$baraza.H5)
plot(density(sc_endline_no_dist$baraza.H5, na.rm=T))
#median is 1
sc_endline_no_dist$baraza.H5.binary <- (sc_endline_no_dist$baraza.H5 > 1)
ols <- lm(baraza.H5.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H5.binary)

sc_endline_no_dist$baraza.H6[sc_endline_no_dist$baraza.H6==999] <- NA
sc_endline_no_dist$baraza.H6 <- as.numeric(as.character(sc_endline_no_dist$baraza.H6))
plot(density(sc_endline_no_dist$baraza.H6, na.rm=T))
ols <- lm(baraza.H6~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H7[sc_endline_no_dist$baraza.H7==999] <- NA
sc_endline_no_dist$baraza.H7 <- as.numeric(as.character(sc_endline_no_dist$baraza.H7))
plot(density(sc_endline_no_dist$baraza.H7, na.rm=T))
ols <- lm(baraza.H7~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H10[sc_endline_no_dist$baraza.H10==999] <- NA
sc_endline_no_dist$baraza.H10 <- as.numeric(as.character(sc_endline_no_dist$baraza.H10))
sc_endline_no_dist$baraza.H10
summary(sc_endline_no_dist$baraza.H10)
table(sc_endline_no_dist$baraza.H10)
plot(density(sc_endline_no_dist$baraza.H10, na.rm=T))
#median is 0
sc_endline_no_dist$baraza.H10.binary <- (sc_endline_no_dist$baraza.H10 > 0)
ols <- lm(baraza.H10.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H10.binary)

sc_endline_no_dist$baraza.H12[sc_endline_no_dist$baraza.H12==999] <- NA
sc_endline_no_dist$baraza.H12 <- as.numeric(as.character(sc_endline_no_dist$baraza.H12))
sc_endline_no_dist$baraza.H12
summary(sc_endline_no_dist$baraza.H12)
table(sc_endline_no_dist$baraza.H12)
plot(density(sc_endline_no_dist$baraza.H12, na.rm=T))
#median is 1
sc_endline_no_dist$baraza.H12.binary <- (sc_endline_no_dist$baraza.H12 > 1)
ols <- lm(baraza.H12.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H12.binary)

sc_endline_no_dist$baraza.H24[sc_endline_no_dist$baraza.H24==999] <- NA
sc_endline_no_dist$baraza.H24 <- as.numeric(as.character(sc_endline_no_dist$baraza.H24))
plot(density(sc_endline_no_dist$baraza.H24, na.rm=T))
ols <- lm(baraza.H24~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H25[sc_endline_no_dist$baraza.H25==999] <- NA
sc_endline_no_dist$baraza.H25 <- as.numeric(as.character(sc_endline_no_dist$baraza.H25))
plot(density(sc_endline_no_dist$baraza.H25, na.rm=T))
ols <- lm(baraza.H25~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H29[sc_endline_no_dist$baraza.H29==999] <- NA
sc_endline_no_dist$baraza.H29 <- as.numeric(as.character(sc_endline_no_dist$baraza.H29))
sc_endline_no_dist$baraza.H29
summary(sc_endline_no_dist$baraza.H29)
table(sc_endline_no_dist$baraza.H29)
plot(density(sc_endline_no_dist$baraza.H29, na.rm=T))
#median is 2
sc_endline_no_dist$baraza.H29.binary <- (sc_endline_no_dist$baraza.H29 > 2)
ols <- lm(baraza.H29.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H29.binary)

sc_endline_no_dist$baraza.H30[sc_endline_no_dist$baraza.H30==999] <- NA
sc_endline_no_dist$baraza.H30 <- as.numeric(as.character(sc_endline_no_dist$baraza.H30))
sc_endline_no_dist$baraza.H30
summary(sc_endline_no_dist$baraza.H30)
table(sc_endline_no_dist$baraza.H30)
plot(density(sc_endline_no_dist$baraza.H30, na.rm=T))
#median is 1
sc_endline_no_dist$baraza.H30.binary <- (sc_endline_no_dist$baraza.H30 > 1)
ols <- lm(baraza.H30.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H30.binary)

sc_endline_no_dist$baraza.H37[sc_endline_no_dist$baraza.H37==999] <- NA
sc_endline_no_dist$baraza.H37 <- as.numeric(as.character(sc_endline_no_dist$baraza.H37))
sc_endline_no_dist$baraza.H37
summary(sc_endline_no_dist$baraza.H37)
table(sc_endline_no_dist$baraza.H37)
plot(density(sc_endline_no_dist$baraza.H37, na.rm=T))
#median is 1
sc_endline_no_dist$baraza.H37.binary <- (sc_endline_no_dist$baraza.H37 > 1)
ols <- lm(baraza.H37.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H37.binary)

sc_endline_no_dist$baraza.H39[sc_endline_no_dist$baraza.H39==999] <- NA
sc_endline_no_dist$baraza.H39 <- as.numeric(as.character(sc_endline_no_dist$baraza.H39))
sc_endline_no_dist$baraza.H39
summary(sc_endline_no_dist$baraza.H39)
table(sc_endline_no_dist$baraza.H39)
plot(density(sc_endline_no_dist$baraza.H39, na.rm=T))
#median is 1
sc_endline_no_dist$baraza.H39.binary <- (sc_endline_no_dist$baraza.H39 > 1)
ols <- lm(baraza.H39.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H39.binary)

sc_endline_no_dist$baraza.H40[sc_endline_no_dist$baraza.H40==999] <- NA
sc_endline_no_dist$baraza.H40 <- as.numeric(as.character(sc_endline_no_dist$baraza.H40))
sc_endline_no_dist$baraza.H40
summary(sc_endline_no_dist$baraza.H40)
table(sc_endline_no_dist$baraza.H40)
plot(density(sc_endline_no_dist$baraza.H40, na.rm=T))
#median is 1
sc_endline_no_dist$baraza.H40.binary <- (sc_endline_no_dist$baraza.H40 > 1)
ols <- lm(baraza.H40.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H40.binary)

sc_endline_no_dist$baraza.H41[sc_endline_no_dist$baraza.H41==999] <- NA
sc_endline_no_dist$baraza.H41 <- as.numeric(as.character(sc_endline_no_dist$baraza.H41))
sc_endline_no_dist$baraza.H41
summary(sc_endline_no_dist$baraza.H41)
table(sc_endline_no_dist$baraza.H41)
plot(density(sc_endline_no_dist$baraza.H41, na.rm=T))
#median is 1
sc_endline_no_dist$baraza.H41.binary <- (sc_endline_no_dist$baraza.H41 > 1)
ols <- lm(baraza.H41.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H41.binary)

sc_endline_no_dist$baraza.H43[sc_endline_no_dist$baraza.H43==999] <- NA
sc_endline_no_dist$baraza.H43 <- as.numeric(as.character(sc_endline_no_dist$baraza.H43))
sc_endline_no_dist$baraza.H43
summary(sc_endline_no_dist$baraza.H43)
table(sc_endline_no_dist$baraza.H43)
plot(density(sc_endline_no_dist$baraza.H43, na.rm=T))
#median is 0
sc_endline_no_dist$baraza.H43.binary <- (sc_endline_no_dist$baraza.H43 > 0)
ols <- lm(baraza.H43.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H43.binary)

sc_endline_no_dist$baraza.H44[sc_endline_no_dist$baraza.H44==999] <- NA
sc_endline_no_dist$baraza.H44 <- as.numeric(as.character(sc_endline_no_dist$baraza.H44))
summary(sc_endline_no_dist$baraza.H44)
table(sc_endline_no_dist$baraza.H44)
#median is 0
sc_endline_no_dist$baraza.H44.binary <- (sc_endline_no_dist$baraza.H44 > 0)
ols <- lm(baraza.H44.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H44.binary)

sc_endline_no_dist$baraza.H46[sc_endline_no_dist$baraza.H46==999] <- NA
sc_endline_no_dist$baraza.H46 <- as.numeric(as.character(sc_endline_no_dist$baraza.H46))
summary(sc_endline_no_dist$baraza.H46)
table(sc_endline_no_dist$baraza.H46)
#median is 1
sc_endline_no_dist$baraza.H46.binary <- (sc_endline_no_dist$baraza.H46 > 1)
ols <- lm(baraza.H46.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H46.binary)

sc_endline_no_dist$baraza.H47[sc_endline_no_dist$baraza.H47==999] <- NA
sc_endline_no_dist$baraza.H47 <- as.numeric(as.character(sc_endline_no_dist$baraza.H47))
summary(sc_endline_no_dist$baraza.H47)
table(sc_endline_no_dist$baraza.H47)
#median is 0
sc_endline_no_dist$baraza.H47.binary <- (sc_endline_no_dist$baraza.H47 > 0)
ols <- lm(baraza.H47.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H47.binary)

sc_endline_no_dist$baraza.H53[sc_endline_no_dist$baraza.H53==999] <- NA
sc_endline_no_dist$baraza.H53 <- as.numeric(as.character(sc_endline_no_dist$baraza.H53))
summary(sc_endline_no_dist$baraza.H53)
table(sc_endline_no_dist$baraza.H53)
#median is 1
sc_endline_no_dist$baraza.H53.binary <- (sc_endline_no_dist$baraza.H53 > 1)
ols <- lm(baraza.H53.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H53.binary)

sc_endline_no_dist$baraza.H54[sc_endline_no_dist$baraza.H54==999] <- NA
sc_endline_no_dist$baraza.H54 <- as.numeric(as.character(sc_endline_no_dist$baraza.H54))
summary(sc_endline_no_dist$baraza.H54)
table(sc_endline_no_dist$baraza.H54)
#median is 1
sc_endline_no_dist$baraza.H54.binary <- (sc_endline_no_dist$baraza.H54 > 1)
ols <- lm(baraza.H54.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H54.binary)

sc_endline_no_dist$baraza.H55[sc_endline_no_dist$baraza.H55==999] <- NA
sc_endline_no_dist$baraza.H55 <- as.numeric(as.character(sc_endline_no_dist$baraza.H55))
summary(sc_endline_no_dist$baraza.H55)
table(sc_endline_no_dist$baraza.H55)
#median is 1
sc_endline_no_dist$baraza.H55.binary <- (sc_endline_no_dist$baraza.H55 > 1)
ols <- lm(baraza.H55.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H55.binary)

sc_endline_no_dist$baraza.H57[sc_endline_no_dist$baraza.H57==999] <- NA
sc_endline_no_dist$baraza.H57 <- as.numeric(as.character(sc_endline_no_dist$baraza.H57))
summary(sc_endline_no_dist$baraza.H57)
table(sc_endline_no_dist$baraza.H57)
#median is 1
sc_endline_no_dist$baraza.H57.binary <- (sc_endline_no_dist$baraza.H57 > 1)
ols <- lm(baraza.H57.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H57.binary)

sc_endline_no_dist$baraza.H58[sc_endline_no_dist$baraza.H58==999] <- NA
sc_endline_no_dist$baraza.H58 <- as.numeric(as.character(sc_endline_no_dist$baraza.H58))
summary(sc_endline_no_dist$baraza.H58)
table(sc_endline_no_dist$baraza.H58)
#median is 1
sc_endline_no_dist$baraza.H58.binary <- (sc_endline_no_dist$baraza.H58 > 1)
ols <- lm(baraza.H58.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H58.binary)

sc_endline_no_dist$baraza.H69[sc_endline_no_dist$baraza.H69==999] <- NA
sc_endline_no_dist$baraza.H69[sc_endline_no_dist$baraza.H69==98] <- NA
sc_endline_no_dist$baraza.H69 <- as.numeric(as.character(sc_endline_no_dist$baraza.H69))
summary(sc_endline_no_dist$baraza.H69)
table(sc_endline_no_dist$baraza.H69)
#median is 1
sc_endline_no_dist$baraza.H69.binary <- (sc_endline_no_dist$baraza.H69 > 1)
ols <- lm(baraza.H69.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H69.binary)

#####SUBSECTION : WATER INFRASTRUCTURE
sc_endline_no_dist$baraza.H72[sc_endline_no_dist$baraza.H72==999] <- NA
sc_endline_no_dist$baraza.H72 <- as.numeric(as.character(sc_endline_no_dist$baraza.H72))
plot(density(sc_endline_no_dist$baraza.H72, na.rm=T))
ols <- lm(baraza.H72~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H77[sc_endline_no_dist$baraza.H77==999] <- NA
sc_endline_no_dist$baraza.H77[sc_endline_no_dist$baraza.H77==98] <- NA
sc_endline_no_dist$baraza.H77 <- as.numeric(as.character(sc_endline_no_dist$baraza.H77))
summary(sc_endline_no_dist$baraza.H77)
table(sc_endline_no_dist$baraza.H77)
#median is 1
sc_endline_no_dist$baraza.H77.binary <- (sc_endline_no_dist$baraza.H77 > 1)
ols <- lm(baraza.H77.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H77.binary)

#####SUBSECTION H3: EDUCATION
sc_endline_no_dist$baraza.H80[sc_endline_no_dist$baraza.H80==999] <- NA
sc_endline_no_dist$baraza.H80 <- as.numeric(as.character(sc_endline_no_dist$baraza.H80))
plot(density(sc_endline_no_dist$baraza.H80, na.rm=T))
ols <- lm(baraza.H80~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H81[sc_endline_no_dist$baraza.H81==999] <- NA
sc_endline_no_dist$baraza.H81 <- as.numeric(as.character(sc_endline_no_dist$baraza.H81))
plot(density(sc_endline_no_dist$baraza.H81, na.rm=T))
ols <- lm(baraza.H81~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H87[sc_endline_no_dist$baraza.H87==999] <- NA
sc_endline_no_dist$baraza.H87 <- as.numeric(as.character(sc_endline_no_dist$baraza.H87))
summary(sc_endline_no_dist$baraza.H87)
table(sc_endline_no_dist$baraza.H87)
#median is 1
sc_endline_no_dist$baraza.H87.binary <- (sc_endline_no_dist$baraza.H87 > 1)
ols <- lm(baraza.H87.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H87.binary)

sc_endline_no_dist$baraza.H91[sc_endline_no_dist$baraza.H91==999] <- NA
sc_endline_no_dist$baraza.H91 <- as.numeric(as.character(sc_endline_no_dist$baraza.H91))
plot(density(sc_endline_no_dist$baraza.H91, na.rm=T))
ols <- lm(baraza.H91~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)

sc_endline_no_dist$baraza.H103[sc_endline_no_dist$baraza.H103==999] <- NA
sc_endline_no_dist$baraza.H103 <- as.numeric(as.character(sc_endline_no_dist$baraza.H103))
summary(sc_endline_no_dist$baraza.H103)
table(sc_endline_no_dist$baraza.H103)
#median is 1
sc_endline_no_dist$baraza.H103.binary <- (sc_endline_no_dist$baraza.H103 > 1)
ols <- lm(baraza.H103.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H103.binary)


sc_endline_no_dist$baraza.H109[sc_endline_no_dist$baraza.H109==999] <- NA
sc_endline_no_dist$baraza.H109 <- as.numeric(as.character(sc_endline_no_dist$baraza.H109))
summary(sc_endline_no_dist$baraza.H109)
table(sc_endline_no_dist$baraza.H109)
#median is 1
sc_endline_no_dist$baraza.H109.binary <- (sc_endline_no_dist$baraza.H109 > 1)
ols <- lm(baraza.H109.binary~information*deliberation+region, data = sc_endline_no_dist)
vcov_cluster <- vcovCR(ols, cluster = sc_endline_no_dist$subcounty, type = "CR0")
res <- coef_test(ols, vcov_cluster)
print(res)
table(sc_endline_no_dist$baraza.H109.binary)