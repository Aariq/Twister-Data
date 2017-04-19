setwd("/Users/scottericr/Google Drive/R Code/Twister Data")
source("/Users/scottericr/Google Drive/R Code/myfunctions.R")

library(ggplot2)
library(tidyr)
timing<-read.csv("timing.csv",header=T)
timing
trunc_time<-timing[1:5,]

ggplot(data=timing, aes(x=Time,y=Vanillin)) + geom_point(size=5,color="#993300") +
  geom_line(size=1, color="#993300")+
  labs(title="Vanillin", x="Time (hours)", y="Relative Peak Area")+
  mytheme

ggplot(data=timing, aes(x=Time,y=Methyl.Salicylate)) + geom_point(size=5,color="#3339CC") +
  geom_line(size=1, color="#3339CC")+
  labs(title="Methyl Salicylate", x="Time (hours)", y="Relative Peak Area")+
  mytheme

ggplot(data=timing, aes(x=Time,y=Hexanal)) + geom_point(size=5,color="#33FF00") +
  geom_line(size=1, color="#33FF00")+
  labs(title="Hexanal", x="Time (hours)", y="Relative Peak Area")+
  mytheme


new_time<-gather(trunc_time, Compound, area, -Time)
new_time$Compound[6:10]<-rep("Methyl Salicylate",5) #just makes the label prettier for the legend later

ggplot(data=new_time, aes(x=Time, y=area)) +geom_point(size=5, aes(color=Compound)) +
  geom_line(size=1, aes(color=Compound)) + labs(x="Time (hours)", y="Relative Peak Area")+
  mytheme +theme(legend.position=c(0.23,0.81))+
  theme(legend.text=element_text(size=18), legend.title=element_text(size=18))

new_time2<-gather(timing, Compound, area, -Time)
new_time2$Compound[7:12]<-rep("Methyl Salicylate",6)
new_time2$Compound[19:24]<-rep("a-Pinene",6)
data<-new_time2[1:18,]

ggplot(data=data, aes(x=Time, y=area)) +geom_point(size=5, aes(color=Compound)) +
  geom_line(size=1, aes(color=Compound)) + labs(x="Time (hours)", y="Relative Peak Area")+
  mytheme +theme(legend.position=c(0.23,0.81))+
  theme(legend.text=element_text(size=18), legend.title=element_text(size=18))

