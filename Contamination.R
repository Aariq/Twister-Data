setwd("/Users/scottericr/Google Drive/R Code/Twister Data")
source("/Users/scottericr/Google Drive/R Code/myfunctions.R")

library(ggplot2)
library(tidyr)
library(dplyr)

####Read in Data####
contam<-read.csv("contamination.csv")
contam$RPA_adj<-contam$RPA*1000
glimpse(contam)
contam
####Figures####
contam1<-contam
contam1$Distance<-as.factor(contam1$Distance)
data<-group_by(contam1,Distance)
data<-summarise(data,mean=mean(RPA),sd=sd(RPA))
data$err.max=data$mean+data$sd
data$err.min=data$mean-data$sd
data
str(data)
plot(contam1$Distance, contam$RPA_adj)

#all plants combined
ggplot(data, aes(as.factor(Distance), mean,ymin=err.min, ymax=err.max)) + 
  geom_bar(stat="identity",fill="dark blue") +
  geom_errorbar(width=.4) +
  labs(x="Distance from plant (cm)", y="Relative peak area") + mytheme +
  theme(axis.title.y=element_text(vjust=1.4),axis.title.x=element_text(vjust=0))

#only left plants
contamleft<-contam[contam$Plant=="left",]
glimpse(contamleft)
dataleft<-group_by(contamleft,Distance)
dataleft<-summarise(dataleft,mean=mean(RPA),sd=sd(RPA))
dataleft$err.max=dataleft$mean+dataleft$sd
dataleft$err.min=dataleft$mean-dataleft$sd
dataleft


ggplot(dataleft, aes(as.factor(Distance), mean,ymin=err.min, ymax=err.max)) + 
  geom_bar(stat="identity",fill="dark blue") +
  geom_errorbar(width=.4) +
  labs(x="Distance from plant (cm)", y="Relative peak area") + mytheme +
  theme(axis.title.y=element_text(vjust=1.4),axis.title.x=element_text(vjust=0))

#only right plants
contamright<-contam[contam$Plant=="right",]
glimpse(contamright)
dataright<-group_by(contamright,Distance)
dataright<-summarise(dataright,mean=mean(RPA),sd=sd(RPA))
dataright$err.max=dataright$mean+dataright$sd
dataright$err.min=dataright$mean-dataright$sd
dataright

contamright[contamright$Distance==150,]
ggplot(dataright, aes(as.factor(Distance), mean,ymin=err.min, ymax=err.max)) + 
  geom_bar(stat="identity",fill="dark blue") +
  geom_errorbar(width=.4) +
  labs(x="Distance from plant (cm)", y="Relative peak area") + mytheme +
  theme(axis.title.y=element_text(vjust=1.4),axis.title.x=element_text(vjust=0))

####Analysis####

contam
glimpse(contam)
mean(contam$RPA)

#regression
m0<-lm(RPA~1,data=contam)
m1<-lm(RPA~Distance, data=contam)
summary(m1) #slope not sig. diff from 0
confint(m1)
plot(contam$Distance, contam$RPA)
points(contam$Distance, predict(m1, type = "response"), type = "l", col = "purple3" , lwd = 3)

#regression with chemical identity as random effect.  
#Not sure that this is appropriate since chemical identiy and distance confounded

library(lme4)
library(car)
m1a<-lmer(RPA_adj~Distance + (1|Chem), data=contam)
summary(m1a)
confint(m1a) #slope overlaps zero, still no effect

#anova
m2<-lm(RPA~Distance, data=contam1)
summary(m2)
anova(m2) #slight effect of distance

m2a<-lm(RPA~Chem*Distance*Plant, data=contam1)
anova(m2a)
step(m2a)

m2b<-lm(RPA~Chem,data=contam1)
anova(m2b)

m2c<-lm(RPA~Distance,data=contam1)
anova(m2c)

#only plants on the left
m3<-lm(RPA~as.factor(Distance), data=contamleft)
summary(m3)
anova(m3)

m4<-lm(RPA~as.factor(Distance), data=contamright)
summary(m4)
anova(m4)

