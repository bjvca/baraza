rm(list=ls())

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

