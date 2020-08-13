rm(list=ls())

### motivating graph
library(ggplot2)
siglev <-  1.96

### this is executed in the /report subdirectory, need to ..
path <- strsplit(getwd(), "/report")[[1]]


dta <- read.csv(paste(path,"data/public/bagezza.csv", sep ="/"), stringsAsFactors = TRUE)
dta <- subset(dta, subcounty != "HAPUYO")

dta$treat <- dta$subcounty =="BAGEZZA"

means <- tapply(dta$dist_source,dta$treat, FUN=mean, na.rm=T)
sds <- tapply(dta$dist_source,dta$treat, FUN=sd, na.rm=T)
sel1_1 <- data.frame(tapply(dta$dist_source,dta$treat, FUN=mean, na.rm=T))
names(sel1_1) <- "mean"
sel1_1$group <- rownames(sel1_1)
sel1_1$up <- means +  siglev*(sds/sqrt(dim(dta)[1]))
sel1_1$down <- means -  siglev*(sds/sqrt(dim(dta)[1]))
sel1_1$time <- 1

means <- tapply(dta$baraza.end.dist,dta$treat, FUN=mean, na.rm=T)
sds <- tapply(dta$baraza.end.dist,dta$treat, FUN=sd, na.rm=T)
sel1_2 <- data.frame(tapply(dta$baraza.end.dist,dta$treat, FUN=mean, na.rm=T))
names(sel1_2) <- "mean"
sel1_2$group <- rownames(sel1_2)
sel1_2$up <- means +   siglev*(sds/sqrt(dim(dta)[1]))
sel1_2$down <- means -  siglev*(sds/sqrt(dim(dta)[1]))
sel1_2$time <- 2

means <- tapply(dta$i1,dta$treat, FUN=mean)
sds <- tapply(dta$i1,dta$treat, FUN=sd)
sel2_1 <- data.frame(tapply(dta$i1,dta$treat, FUN=mean))
names(sel2_1) <- "mean"
sel2_1$group <- rownames(sel2_1)
sel2_1$up <- means +  siglev*(sds/sqrt(dim(dta)[1]))
sel2_1$down <- means -   siglev*(sds/sqrt(dim(dta)[1]))
sel2_1$time <- 1

means <- tapply(dta$baraza.H1,dta$treat, FUN=mean)
sds <- tapply(dta$baraza.H1,dta$treat, FUN=sd)
sel2_2 <- data.frame(tapply(dta$baraza.H1,dta$treat, FUN=mean))
names(sel2_2) <- "mean"
sel2_2$group <- rownames(sel2_2)
sel2_2$up <- means +  siglev*(sds/sqrt(dim(dta)[1]))
sel2_2$down <- means -  siglev*(sds/sqrt(dim(dta)[1]))
sel2_2$time <- 2

### now for store

means <- tapply(dta$i2,dta$treat, FUN=mean)
sds <- tapply(dta$i2,dta$treat, FUN=sd)
store1_1 <- data.frame(tapply(dta$i2,dta$treat, FUN=mean))
names(store1_1) <- "mean"
store1_1$group <- rownames(store1_1)
store1_1$up <- means +  siglev*(sds/sqrt(dim(dta)[1]))
store1_1$down <- means -  siglev*(sds/sqrt(dim(dta)[1]))
store1_1$time <- 1

means <- tapply(dta$baraza.H2,dta$treat, FUN=mean)
sds <- tapply(dta$baraza.H2,dta$treat, FUN=sd)
store1_2 <- data.frame(tapply(dta$baraza.H2,dta$treat, FUN=mean))
names(store1_2) <- "mean"
store1_2$group <- rownames(store1_2)
store1_2$up <- means +  siglev*(sds/sqrt(dim(dta)[1]))
store1_2$down <- means -  siglev*(sds/sqrt(dim(dta)[1]))
store1_2$time <- 2

means <- tapply(dta$base_unprotected,dta$treat, FUN=mean)
store2_1 <- data.frame(tapply(dta$base_unprotected,dta$treat, FUN=mean))
names(store2_1) <- "mean"
store2_1$group <- rownames(store2_1)
store2_1$up <- means +  siglev*sqrt(means*(1-means)/dim(dta)[1])
store2_1$down <- means -  siglev*sqrt(means*(1-means)/dim(dta)[1])
store2_1$time <- 1

means <- tapply(dta$unprotected,dta$treat, FUN=mean)
store2_2 <- data.frame(tapply(dta$unprotected,dta$treat, FUN=mean))
names(store2_2) <- "mean"
store2_2$group <- rownames(store2_2)
store2_2$up <- means +  siglev*sqrt(means*(1-means)/dim(dta)[1])
store2_2$down <- means -  siglev*sqrt(means*(1-means)/dim(dta)[1])
store2_2$time <- 2





sel1 <- rbind(sel1_1, sel1_2)
sel2 <- rbind(sel2_1, sel2_2)

store1 <- rbind(store1_1, store1_2)
store2 <- rbind(store2_1, store2_2)

sel1$time <- as.factor(sel1$time)
sel1$group <- as.factor(sel1$group)
levels(sel1$group) <- c("Ctrl","Treat")

sel2$time <- as.factor(sel2$time)
sel2$group <- as.factor(sel2$group)
levels(sel2$group) <- c("Ctrl","Treat")

store1$time <- as.factor(store1$time)
store1$group <- as.factor(store1$group)
levels(store1$group) <- c("Ctrl","Treat")

store2$time <- as.factor(store2$time)
store2$group <- as.factor(store2$group)
levels(store2$group) <- c("Ctrl","Treat")


p1 <- ggplot(sel1, aes(x=time, y=mean, group=group)) + 
  geom_line(aes(linetype=group)) +
  geom_point()+
  geom_pointrange(aes(ymin=up, ymax=down))+ ylab("kilometers") + ggtitle("Distance to water source", subtitle = NULL)+ theme(plot.title = element_text(hjust = 0.5))


p2 <- ggplot(sel2, aes(x=time, y=mean, group=group)) + 
 geom_line(aes(linetype=group)) +
  geom_point()+
  geom_pointrange(aes(ymin=up, ymax=down))+ ylab("scale (1-10)")+ ggtitle("Access to water is a problem", subtitle = NULL)+ theme(plot.title = element_text(hjust = 0.5))


p3 <- ggplot(store1, aes(x=time, y=mean, group=group)) + 
 geom_line(aes(linetype=group)) +
  geom_point()+
  geom_pointrange(aes(ymin=up, ymax=down))+ ylab("scale (1-10)")+ ggtitle("Drinking water is usually dirty", subtitle = NULL) + theme(plot.title = element_text(hjust = 0.5))

p4 <- ggplot(store2, aes(x=time, y=mean, group=group)) + 
 geom_line(aes(linetype=group)) +
  geom_point()+
  geom_pointrange(aes(ymin=up, ymax=down))+ ylab("proportion of households")+ ggtitle("Use unprotected source", subtitle = NULL) + theme(plot.title = element_text(hjust = 0.5))

require(gridExtra)
png(paste(path,"report/figure/bagezza.png", sep="/"), units="px", height=3200, width= 4200, res=600)
grid.arrange(p1, p2, p3, p4, ncol=2)
dev.off()
