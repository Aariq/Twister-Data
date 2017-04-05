setwd("/Users/scottericr/Google Drive/R Code/Twister Data")
#setwd("/Users/escott03/Google Drive/R Code/Twister Data")
library(reshape2)
library(dplyr)
library(ggplot2)
library(car)
library(scales)
library(MASS)
library(RColorBrewer)

source("/Users/scottericr/Google Drive/R Code/myfunctions.R")

#####Read in data#######
twister = read.csv("Field.csv")
glimpse(twister)
twister
twist2 = recast(twister, variable~Compound, measure.var=c(2:12)) #this transposes the data so each column is a compound
glimpse(twist2)
twist3<-twist2
twist3[is.na(twist3)]=0 #missing cells represent undetected compounds
colnames(twist3)[1]<-"samples"
glimpse(twist3)
summary(twist3$samples)
treatment<-c(rep("E. obliqua",3),rep("MeJA",4),rep("control",4))

###################################################################################################


####Exploratory stats####
#how many compounds detected total?
dim(twist2)-1 #149

#how many compounds detected in each treatment?
library(venn)
removednas<-twister
removednas[is.na(removednas)]=0 #missing cells represent undetected compounds

x<-transmute(removednas,Control=con1+con2+con3+con4, MeJA=meja1+meja2+meja3+meja4, Geometrid=geo1+geo2+geo3)
x[x>0]=1
x
venn(x, cexil=1.5, cexsn=1.1)
sum(x$Control)
sum(x$Geometrid)
sum(x$MeJA)
#unique to Geometrid: 
twister$Compound[c(5,70,114,134)]

#unique to MeJA:
twister[123,]

#-0.0378833907  3.147150e-02  2.403748e-02

####PCA####
tea.pca2<-prcomp(twist3[2:length(twist3)],scale.=TRUE, center = TRUE)
summary(tea.pca2)
tea.pca2$rotation
screeplot(tea.pca2, type="lines")

####PCA with Pareto scaling####
pareto_scale<- function(x) {
  apply(x, 2, function(y) (y - mean(y))/sqrt(sd(y)))
}

pareto.data<-pareto_scale(twist3[2:length(twist3)])

pareto.pca<-prcomp(pareto.data,scale=FALSE, center=FALSE)
summary(pareto.pca)
screeplot(pareto.pca, type="lines")


######## Try with FactoMineR#######
library(FactoMineR)
tea.pca3<-PCA(twist3[2:length(twist3)], scale.unit = TRUE)
tea.pca3
summary(tea.pca3)
tea.pca3$var$cor
tea.pca3$var$contrib
dimdesc(tea.pca3,axes = 1:2, proba = 0.05)
plot(tea.pca3, select = "contrib 10")
###I think basically the same results?#

####Plotting PCA####
#Treatment<-c(rep("Control",4),rep("MeJA",4),rep("Geometrid",3))
labels = expression("Control", italic("E. obliqua"), "MeJA")
ggplot(as.data.frame(tea.pca2$x), aes(PC1, PC2)) + 
  geom_point(aes(color = treatment), size = I(4)) + 
  scale_color_brewer(palette=("Set1"),labels = labels, name = NULL) +
  newphytol.theme + 
  theme(legend.position=c(0.13,0.12),
  legend.text=element_text(size=12),
  legend.background = element_rect(color = "black")) +
  ylab("PC2 (16.86%)") + xlab("PC1 (42.97%)")


ggplot(as.data.frame(pareto.pca$x), aes(PC1, PC2)) +
  geom_point(aes(color = treatment), size = I(5)) +
  mytheme 
  

####Output correlations####
twist.scaled<-as.data.frame(scale(twist3[2:length(twist3)]))
#Cors.2<-PC.cors(tea.pca2,twist3[2:length(twist3)],2)
Cors.2<-PC.cors(tea.pca2, twist.scaled, 2)
Cors.2

PC1.ordered<-format(Cors.2[order(Cors.2$PC1.corr.coefs),],digits=2,scientific=FALSE)
PC1.ordered
write.csv(PC1.ordered, file="PC1correlations.csv")

PC2.ordered<-format(Cors.2[order(Cors.2$PC2.corr.coefs),],digits=2,scientific=FALSE)
PC2.ordered
write.csv(PC2.ordered, file="PC2correlations.csv")

#example showing why isovaleric acid is correlated with PC2, but not PC1 
#even though it is absent from controls and controls are others are separated along PC1
cor.test(twist3[,"Isovaleric acid"],tea.pca2$x[,1])
cor.test(twist3[,"Isovaleric acid"],tea.pca2$x[,2])
opar<-par(mfrow=c(1,2))
plot(twist3[,"Isovaleric acid"]~tea.pca2$x[,1], xlab="PC1 values", ylab="Isovaleric acid peak area")
plot(twist3[,"Isovaleric acid"]~tea.pca2$x[,2], xlab="PC2 values", ylab="Isovaleric acid peak area")
par(opar)


####Anova####
m0<-lm(tea.pca2$x[,1]~treatment)
summary(m0)
Anova(m0)

m1<-lm(tea.pca2$x[,2]~treatment)
summary(m1)
Anova(m1)

m2<-lm(tea.pca2$x[,3]~treatment)
summary(m2)
Anova(m2)

#Calculate expected probabilities
Fams.count<-count(twister$Family)
p<-Fams.count$freq/sum(Fams.count$freq)
p

PC1.fams<-Cors.2$Family[Cors.2$PC1.p.values<=0.01]
PC1.fams.count<-merge(Fams.count$x,count(PC1.fams),all.y=T,all.x=T)
PC1.fams.count[is.na(PC1.fams.count)]<-0

chisq.test(PC1.fams.count$freq,p=p)

###PC1 not enriched in any particular family of compounds

PC2.fams<-Cors.2$Family[Cors.2$PC2.p.values<=0.01]
PC2.fams.count<-merge(Fams.count$x,count(PC2.fams),all.y=T,all.x=T)
PC2.fams.count[is.na(PC2.fams.count)]<-0

chisq.test(PC2.fams.count$freq,p=p)

###PC2 not enriched in any particular family of compounds)


#####################

#Do some univariate stats on treatment~no.compounds or total peak area or total relative peak area
#presence/absence of compounds with binomial glm?

##calculating peak area
head(twist3)
#add a column with sums across rows
data<-twist3

data$total.area=rowSums(twist3[2:length(twist3)])
head(data)
data$treatment<-as.factor(c(rep("Control",4),rep("MeJA",4),rep("Geometrid",3))) 

tea.glm<-glm(total.area~treatment,data=data,family=gaussian())

summary(tea.glm)
Anova(tea.glm)
tea.null<-glm(total.area~1,data=data,family=gaussian())
library(lmtest)
lrtest(tea.null,tea.glm)


tea.glm2<-glm(total.area~-1+treatment,data=data,family=gaussian())
coef(tea.glm2)


scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scientific_format()(x)))
}

ggplot(aes(treatment, total.area),data=data)+geom_boxplot()+ylab("Total Peak Area")+xlab(NULL)+mytheme + scale_y_continuous(labels=scientific_10)

#######total peak area does vary by treatment.  highest in MeJA####


#####number of compounds
just.data<-data[2:(length(data)-2)]
count.data<-as.data.frame(ifelse(just.data>0,1,0))
count.data$total.count<-rowSums(count.data)
count.data$treatment<-data$treatment
head(count.data)

count.Control<-count.data[count.data$treatment=="Control"]
count.Control

count.MeJA<-count.data[count.data$treatment=="MeJA"]
count.MeJA


tea.glm3<-glm(total.count~treatment,data=count.data,family=poisson(link="log"))
summary(tea.glm3)
Anova(tea.glm3)

tea.glm4<-glm(total.count~-1+treatment,data=count.data,family=poisson(link="log"))
exp(coef(tea.glm4))

ggplot(aes(treatment, total.count),data=count.data)+geom_boxplot()+ylab("Number of Compounds Detected")+xlab(NULL)+mytheme

####no effect of treatment on number of compounds, 
####but a trend toward more compounds in MeJA treated plants

counts.above<-as.data.frame(ifelse(just.data>10000,1,0))
counts.above$total.count<-rowSums(counts.above)
counts.above$treatment<-data$treatment
head(counts.above)

tea.glm5<-glm(total.count~treatment,data=counts.above,family=poisson)
summary(tea.glm5)
Anova(tea.glm5)
tea.glm6<-glm(total.count~-1+treatment,data=counts.above,family=poisson)
exp(coef(tea.glm6))
#####It doesn't change the results if I look at compounds present above some threshold



library(mixOmics)
dim(twist3)
X<-twist3[2:155]
Y<-Treatment

#uses only 50 variables per axis
ncomp=4
result<-splsda(X, Y, ncomp=ncomp, keepX=rep(50,ncomp))
result
plotcols<-c(brewer.pal(3,"Dark2"))
plotIndiv(result, comp=1:2, ind.names=F,add.legend = T,
          plot.ellipse=T,cex=4,col.per.group = plotcols, 
          pch=c(15:17), style="ggplot2", 
          main="", X.label="PLS1", Y.label = "PLS2")


plotVar(result, comp = 1:2, cex=3,var.names=T,cutoff=0.75,style = "ggplot2")
error<-perf(result,validation= "loo",method.predict="all",near.zero.var = F)
plot(error)

#without variable selection there is less separation and #9 is by itself
#result.pls<-splsda(X,Y,ncomp=3)
#plotIndiv(result.pls, comp=1:2, ind.names=F,add.legend = T)

test<-PLSDA.cors(result,twist3[2:155],2) #from MyFUnctions

test2<-format(test[order(test[,3]),],digits=2,scientific=FALSE)
test2

plotVar(result, comp = 1:2, cex=3,var.names=T,cutoff=0.95,style = "ggplot2")

#p<-#a vector of length X (number of compounds) for coloring by chemical family
#test with dummy variables
#p=rep(1:4, length.out=length(X))
#plotVar(result, comp = 1:2, cex=3,var.names=T,cutoff=0.75,style = "ggplot2",col=list(p))


######the perf() function keeps getting numerical errors, so I maybe have to do manual
######cross validation using predict() with a training set:
  i <- 1
samp <- sample(1:3, nrow(X), replace = TRUE) # Creation of a list of the same size
# as X containing 1, 2 or 3

test <- which(samp == i) # Search which column in samp has a value of 1
train <- setdiff(1:nrow(X), test) # Keeping the column that are not in test

## For sPLS-DA
splsda.train <- splsda(X[train, ], Y[train], ncomp = 3, keepX = c(50, 50, 50))
test.predict <- predict(splsda.train, X[test,], method = "max.dist")
Prediction <-levels(Y)[test.predict$class$max.dist[,1]]
cbind(Y = as.character(Y[test]), Prediction)






## Set number of components to 2
ncomp <- 2

## Total number of selected compunds on all ncomp dimensions
keepX <- round(c(seq(5, 45, 5), seq(50, 500, 50))/ncomp)
error <- matrix(NA, nrow = length(keepX), ncol = 3) 

for (i in 1:length(keepX)) {
  error[i, ] <- valid(X, Y, ncomp = 3, keepX = rep(keepX[i], ncomp),
                      pred.method = "max.dist", method = "splsda",
                      validation = "Mfold", M = 10)
}                   


## Plot the error obtained for each dimension
matplot(error, type = 'l', axes = FALSE, xlab = 'number of selected genes',
        ylab = 'error rate', col = c("black", "red", "blue"), lwd = 2, lty = 1)
axis(1, c(1:length(keepX)), labels = keepX)
axis(2)
legend(6, 0.45, lty = 1, legend = c('dim 1', 'dim 1:2', 'dim 1:3'),
       horiz = TRUE, cex = 0.9, col = c("black", "red", "blue"), lwd = 2)



