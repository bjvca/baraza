rm(list=ls())
#remove and rm can be used to remove objects
#list: a character vector naming objects to be removed
#ls shows what data sets and functions a user has defined

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
sc_endline_no_dist$baraza.production.E1a
summary(lm(baraza.production.E1a~ind_treat, data = sc_endline_no_dist))
sc_endline_no_dist$baraza.health.E2a
summary(lm(baraza.health.E2a~ind_treat, data = sc_endline_no_dist))
sc_endline_no_dist$baraza.gender1.E3a
summary(lm(baraza.gender1.E3a~ind_treat, data = sc_endline_no_dist))
sc_endline_no_dist$baraza.works.E4a
summary(lm(baraza.works.E4a~ind_treat, data = sc_endline_no_dist))
sc_endline_no_dist$baraza.finance1.E5a
summary(lm(baraza.finance1.E5a~ind_treat, data = sc_endline_no_dist))
#how to measure difference date and today
sc_endline_no_dist$baraza.finance1.E5a
summary(lm(baraza.finance1.E5a~ind_treat, data = sc_endline_no_dist))

########Caroline's first R day, after call at 4pm
#sc, information and deliberation
ols <- lm( baraza.H1~ information * deliberation + region, data = sc_endline[sc_endline$district_baraza == 0,] )
summary(ols)
coef(ols)