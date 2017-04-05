setwd("/Users/scottericr/Google Drive/R Code/Twister Data")
library(reshape2)
library(tidyr)
library(dplyr)
library(venn)
library(ggplot2)
source("/Users/scottericr/Google Drive/R Code/myfunctions.R")


#####Read in data#######
#contam<-read.csv("DHS Comparison.csv")
comp.data <- read.csv("DHS Comparison updated.csv") #Nicole added back two compounds she previously removed as contaminants.  They're not actually contaminants.
head(comp.data)
comp.data[is.na(comp.data)]=0  #NA's represent 0's
glimpse(comp.data)

#create treatment variables
vars <- colsplit(comp.data$Treatment, "\\.", c("Method", "Treatment", "Rep"))
vars

dim(comp.data)
data <- cbind(vars$Method, vars$Treatment, comp.data[,2:291])
colnames(data)[1:2] <- c("Method","Treatment")
glimpse(data[,1:15])

#find any compounds that are constants
counts<-head(summarise_each(data, funs(n_distinct)))
glimpse(counts)
rm <- which(counts <= 1)
rm
data <- subset(data, select = -rm)
#########################

####Venn Diagram####
x<-as.data.frame(t(data[3:length(data)]))
colnames(x) <- contam$Treatment
head(x)
x.logical<-transmute(x, Tenax = Tenax.Cont.1 + Tenax.Cont.2 + Tenax.Cont.3 + 
                Tenax.MeJA.1 + Tenax.MeJA.2 + Tenax.MeJA.3, 
                Twister = Twister.Cont.1 + Twister.Cont.2 + Twister.Cont.3 + 
                Twister.MeJA.1 + Twister.MeJA.2 + Twister.MeJA.3,
                PDMS = PDMS.Cont.1 + PDMS.Cont.2 + PDMS.Cont.3 +
                PDMS.MeJA.1 + PDMS.MeJA.2 + PDMS.MeJA.3,
                Strip = Strip.Cont.1 + Strip.Cont.2 + Strip.Cont.3 + 
                Strip.MeJA.1 + Strip.MeJA.2)
x.logical[x.logical>0] = 1
tail(x.logical)
head(x)
venn(x.logical, cexil = 1, cexsn = .85)

##Total##
colSums(x.logical)

#Find the compounds unique to each method
Tenax.unique <- subset(x, x.logical$Twister == 0 & x.logical$PDMS == 0 & x.logical$Strip == 0)
dim(Tenax.unique)


Twister.unique <- subset(x, x.logical$Tenax == 0 & x.logical$PDMS == 0 & x.logical$Strip == 0)
dim(Twister.unique)

PDMS.unique <- subset(x, x.logical$Tenax == 0 & x.logical$Twister == 0 & x.logical$Strip == 0)
dim(PDMS.unique)

Strip.unique <- subset(x, x.logical$Tenax == 0 & x.logical$Twister == 0 & x.logical$PDMS == 0)
dim(Strip.unique)

#in BOTH DHS methods but not in either DCSE method
DHS.unique.wrong <- subset(x, x.logical$Tenax == 1 & x.logical$PDMS == 1 & x.logical$Twister == 0 & x.logical$Strip == 0)
dim(DHS.unique.wrong)


DCSE.unique.wrong <- subset(x, x.logical$Tenax == 0 & x.logical$PDMS == 0 & x.logical$Twister == 1 & x.logical$Strip == 1)
dim(DCSE.unique.wrong)

#what you should do is add a column identifying all samples as DCSE or DHS, then do as above.  Compare with changing & to | with dim() to check

glimpse(x)
head(x)
glimpse(x.logical)
#Shit.  Gotta go back to before x.logical was created.
##By headspace vs. direct contact
x.method<-transmute(x, DHS = Tenax.Cont.1 + Tenax.Cont.2 + Tenax.Cont.3 + 
                       Tenax.MeJA.1 + Tenax.MeJA.2 + Tenax.MeJA.3 + 
                      PDMS.Cont.1 + PDMS.Cont.2 + PDMS.Cont.3 +
                      PDMS.MeJA.1 + PDMS.MeJA.2 + PDMS.MeJA.3, 
                     DCSE = Twister.Cont.1 + Twister.Cont.2 + Twister.Cont.3 + 
                       Twister.MeJA.1 + Twister.MeJA.2 + Twister.MeJA.3 + 
                      Strip.Cont.1 + Strip.Cont.2 + Strip.Cont.3 + 
                       Strip.MeJA.1 + Strip.MeJA.2)
x.method[x.method>0] = 1
tail(x.method)
venn(x.method, cexil = 1, cexsn = .85)

DHS.unique <- subset(x, x.method$DHS != 0 & x.method$DCSE == 0)
dim(DHS.unique)

DCSE.unique <- subset(x, x.method$DCSE !=0 & x.method$DHS == 0)
dim(DCSE.unique)

##By treatment
x.treatment<-transmute(x, Control = Tenax.Cont.1 + Tenax.Cont.2 + Tenax.Cont.3 + 
                       Twister.Cont.1 + Twister.Cont.2 + Twister.Cont.3 + 
                       PDMS.Cont.1 + PDMS.Cont.2 + PDMS.Cont.3 +
                       Strip.Cont.1 + Strip.Cont.2 + Strip.Cont.3,
                      MeJA = Tenax.MeJA.1 + Tenax.MeJA.2 + Tenax.MeJA.3 + 
                       Twister.MeJA.1 + Twister.MeJA.2 + Twister.MeJA.3 +
                       PDMS.MeJA.1 + PDMS.MeJA.2 + PDMS.MeJA.3 +
                       Strip.MeJA.1 + Strip.MeJA.2)

head(x.treatment)
x.trt.logical <- x.treatment
x.trt.logical[x.trt.logical>0] = 1
colSums(x.trt.logical)

Control.unique <- subset(x, x.treatment$Control != 0 & x.treatment$MeJA == 0)
dim(Control.unique)

MeJA.unique <- subset(x, x.treatment$Control == 0 & x.treatment$MeJA != 0)
dim(MeJA.unique)

U.compounds <- list(Tenax = rownames(Tenax.unique), Twister = rownames(Twister.unique), 
                    PDMS = rownames(PDMS.unique), Strip = rownames(Strip.unique),
                    DHS = rownames(DHS.unique), DCSE = rownames(DCSE.unique),
                    Control = rownames(Control.unique), MeJA = rownames(MeJA.unique))
#maybe try to figure out how to make this list prettier
#writes to file
sink("Unique Compounds.txt")
U.compounds
sink()

####PCA####
data.scaled <- scale(data[3:length(data)])
all.pca <- prcomp(data.scaled, scale. = FALSE, center = FALSE)
summary(all.pca)
#6 or 7 PCs to capture ~80% of variance
screeplot(all.pca,type="l")

Method <- data$Method
Treatment <- data$Treatment

#I can't figure out how else to change the text in the legend

shape.labels = c("Control", "MeJA")
col.labels = expression("Tenax DHS", "PDMS DHS", "Strip DCSE", paste(Twister^"®"," DCSE"))
cols<-c("PDMS" = "#FF5C00", "Tenax" = "#FFC000", "Strip" = "#00A779", "Twister" = "#1D1AB2")
ggplot(as.data.frame(all.pca$x), aes(PC1,PC2)) + 
  geom_point(aes (color=Method, shape=Treatment),size=I(4)) + 
  newphytol.theme + 
  scale_shape(labels = shape.labels, name = "Treatment") +
  theme(legend.position=c(0.85, 0.75))+
  theme(legend.text=element_text(size = 10), legend.title=element_text(size = 10),
        legend.background = element_rect(color = "black")) +
  #scale_color_brewer(palette="Set1") + 
  scale_color_manual(values = cols, breaks = c("Tenax", "PDMS", "Strip", "Twister"),
                     labels = col.labels) + ylim(-10,20) +
  ylab("PC2 (19.49%)") + xlab("PC1 (26.46%)")


pca.correlations<-PC.cors(all.pca, as.data.frame(data.scaled), 2)
pca.correlations


write.csv(pca.correlations, "DHS comparison correlations.csv")

#try with pareto scaling, which gives less importance to noise####
p.scaled<-pareto_scale(data[3:length(data)]) #from myfunctions.R
pareto.pca <- prcomp(p.scaled, scale. = FALSE, center = FALSE)
summary(pareto.pca)

ggplot(as.data.frame(pareto.pca$x), aes(PC1,PC2)) + 
  geom_point(aes (color=Method, shape=Treatment),size=I(5)) + 
  mytheme+ 
  theme(legend.position=c("bottom"))+
  theme(legend.text=element_text(size=18), legend.title=element_text(size=18))+
  scale_color_brewer(palette="Dark2") #+ ylab("PC2 (19.47%)") + xlab("PC1 (26.55%)")

#unsurprisingly, separation between methods is less good with pareto scaling, 
#since the differences between methods is probably due to low concentration compounds and pareto scaling doesn't give as much importance to small peaks
###########

##subset of just Tenax DHS and Twister DCSE####
data.sub1 <- subset(data, Method == "Tenax" | Method == "Twister")
dim(data.sub1)
glimpse(data.sub1)
#look for constants and remove
counts.sub1<-head(summarise_each(data.sub1, funs(n_distinct)))
glimpse(counts.sub1)
rm <- which(counts.sub1 <= 1)
rm
data.sub1 <- subset(data.sub1, select = -rm)

sub1.pca <- prcomp(data.sub1[3:length(data.sub1)], scale. = TRUE, center = TRUE)
summary(sub1.pca)

ggplot(as.data.frame(sub1.pca$x), aes(PC1,PC2)) + 
  geom_point(aes (color=data.sub1$Method, shape=data.sub1$Treatment),size=I(5)) + 
  mytheme+ 
  theme(legend.position=c("bottom"))+
  theme(legend.text=element_text(size=18), legend.title=element_text(size=18))+
  scale_color_brewer(palette="Dark2") #+ ylab("PC2 (19.47%)") + xlab("PC1 (26.55%)")
##better separation with Twister DCSE

###Subset with just control plants####
data.sub2 <- subset(data, Treatment == "Cont")
counts.sub2<-head(summarise_each(data.sub2, funs(n_distinct)))
glimpse(counts.sub2)
rm <- which(counts.sub2 <= 1)
rm
data.sub2 <- subset(data.sub2, select = -rm)

sub2.pca <- prcomp(data.sub2[3:length(data.sub2)], scale. = TRUE, center = TRUE)

ggplot(as.data.frame(sub2.pca$x), aes(PC1,PC2)) + 
  geom_point(aes (color=data.sub2$Method),size=I(5)) + 
  mytheme+ 
  theme(legend.position=c("bottom"))+
  theme(legend.text=element_text(size=18), legend.title=element_text(size=18))+
  scale_color_brewer(palette="Dark2") #+ ylab("PC2 (19.47%)") + xlab("PC1 (26.55%)")

###subset with just MeJA plants####
data.sub3 <- subset(data, Treatment == "MeJA")
counts.sub3<-head(summarise_each(data.sub3, funs(n_distinct)))
glimpse(counts.sub3)
rm <- which(counts.sub3 <= 1)
rm
data.sub3 <- subset(data.sub3, select = -rm)

sub3.pca <- prcomp(data.sub3[3:length(data.sub3)], scale. = TRUE, center = TRUE)

ggplot(as.data.frame(sub3.pca$x), aes(PC1,PC2)) + 
  geom_point(aes (color=data.sub3$Method),size=I(5)) + 
  mytheme+ 
  theme(legend.position=c("bottom"))+
  theme(legend.text=element_text(size=18), legend.title=element_text(size=18))+
  scale_color_brewer(palette="Dark2") #+ ylab("PC2 (19.47%)") + xlab("PC1 (26.55%)")

###subsets for each method####
data.Tenax <- subset(data, Method == "Tenax")
counts.Tenax<-head(summarise_each(data.Tenax, funs(n_distinct)))
glimpse(counts.Tenax)
rm <- which(counts.Tenax <= 1)
rm
data.Tenax <- subset(data.Tenax, select = -rm)

data.PDMS <- subset(data, Method == "PDMS")
counts.PDMS<-head(summarise_each(data.PDMS, funs(n_distinct)))
glimpse(counts.PDMS)
rm <- which(counts.PDMS <= 1)
rm
data.PDMS <- subset(data.PDMS, select = -rm)

data.Strip <- subset(data, Method == "Strip")
counts.Strip <- head(summarise_each(data.Strip, funs(n_distinct)))
glimpse(counts.Strip)
rm <- which(counts.Strip <= 1)
rm
data.Strip <- subset(data.Strip, select = -rm)

data.Twister <- subset(data, Method == "Twister")
counts.Twister <- head(summarise_each(data.Twister, funs(n_distinct)))
glimpse(counts.Twister)
rm <- which(counts.Twister <= 1)
rm
data.Twister <- subset(data.Twister, select = -rm)

Tenax.pca <- prcomp(data.Tenax[3:length(data.Tenax)], scale. = TRUE, center = TRUE)
PDMS.pca <- prcomp(data.PDMS[3:length(data.PDMS)], scale. = TRUE, center = TRUE)
Strip.pca <- prcomp(data.Strip[3:length(data.Strip)], scale. = TRUE, center = TRUE)
Twister.pca <- prcomp(data.Twister[3:length(data.Twister)], scale. = TRUE, center = TRUE)

Tenax.plot <- ggplot(as.data.frame(Tenax.pca$x), aes(PC1,PC2)) + 
  geom_point(aes (shape=data.Tenax$Treatment),size=I(5)) + 
  mytheme + theme(legend.position=c("bottom"))

PDMS.plot <- ggplot(as.data.frame(PDMS.pca$x), aes(PC1,PC2)) + 
  geom_point(aes (shape=data.PDMS$Treatment),size=I(5)) + 
  mytheme + theme(legend.position=c("bottom"))

Strip.plot <- ggplot(as.data.frame(Strip.pca$x), aes(PC1,PC2)) + 
  geom_point(aes (shape=data.Strip$Treatment),size=I(5)) + 
  mytheme + theme(legend.position=c("bottom"))

Twister.plot <- ggplot(as.data.frame(Twister.pca$x), aes(PC1,PC2)) + 
  geom_point(aes (shape=data.Twister$Treatment),size=I(5)) + 
  mytheme + theme(legend.position=c("bottom"))

opar<-par(mfrow=c(2,2))
Tenax.plot
PDMS.plot
Strip.plot
Twister.plot
par(opar)

###only compounds unique to DCSE
#subset by method
data.DCSE<- t(DCSE.unique)
data.DCSE <- subset(data.DCSE, vars$Method == "Twister")
#subset by list of compounds
data.DCSE
dim(data.DCSE)

scaled.DCSE<- scale(data.DCSE)
unique.PCA <- prcomp(scaled.DCSE, scale. = FALSE)
summary(unique.PCA)
screeplot(unique.PCA,type="l")

data.DCSE
#Method <- c(rep("Twister",6),rep("Strip",5))
#Treatment <- c(rep("Control", 3), rep("MeJA", 3), rep("Control",3), rep("MeJA", 2))
Treatment <- c(rep("Control", 3), rep("MeJA", 3))
ggplot(as.data.frame(unique.PCA$x), aes(PC1,PC2)) + 
  geom_point(aes(color=Treatment),size=I(4)) + 
  newphytol.theme #+ 
  #scale_shape(labels = shape.labels, name = "Treatment") +
  #theme(legend.position=c(0.85, 0.75))+
  #theme(legend.text=element_text(size = 10), legend.title=element_text(size = 10),
  #      legend.background = element_rect(color = "black")) +
  #scale_color_brewer(palette="Set1") + 
  #scale_color_manual(values = cols, breaks = c("Tenax", "PDMS", "Strip", "Twister"),
   #                  labels = col.labels) +
  #ylab("PC2 (19.47%)") + xlab("PC1 (26.55%)")

##Next steps:
#  -Decide on making a table with 14 rows or including this PCA (leaning toward table)
#  -Figure out appropriate correction for multiple testing
#  -Do t-tests on all 14 unique-to-DCSE compounds to look for important differences we are missing with DHS
#  -Figure out if those differentially abundant metabolites have known biological importance
#  -redo with read_csv() from the readr package (doesn't fuck with column names)



## t-Tests for DCSE unique compounds
glimpse(data.DCSE)
control.data <- data.DCSE[Treatment == "Control",]
MeJA.data <- data.DCSE[Treatment == "MeJA",]
t.test(control.data[,14],MeJA.data[,14])
#sample size too small to detect many differences :-(
head(control.data)

comparisons = dim(control.data)[2]
t.table <- data_frame(Compound = colnames(control.data), stat = NA, p.value = NA, Control.means = NA, MeJA.means = NA)
t.table
for(i in 1:comparisons){
  temp <- t.test(control.data[ ,i], MeJA.data[ ,i])
  t.table$p.value[i] <- temp$p.value
  t.table$stat[i] <- temp$statistic
  t.table$Control.means[i] <- colMeans(control.data)[i]
  t.table$MeJA.means[i] <- colMeans(MeJA.data)[i]
}
t.table$p.adj <- p.adjust(t.table$p.value, method = "bonferroni")
library(knitr)
t.table[t.table$p.value <= 0.05,] %>%
  kable(digits = 4)

## t-Tests for DHS unique compounds
data.DHS<- t(DHS.unique)
glimpse(data.DHS)
control.DHS <- data.DHS[Treatment == "Control",]
MeJA.DHS <- data.DHS[Treatment == "MeJA",]

comparisons = dim(control.DHS)[2]
t.table <- data_frame(Compound = colnames(control.DHS), p.value = NA, Control.means = NA, MeJA.means = NA)
t.table
for(i in 1:comparisons){
  temp <- t.test(control.DHS[ ,i], MeJA.DHS[ ,i])
  t.table$p.value[i] <- temp$p.value
  t.table$Control.means[i] <- colMeans(control.data)[i]
  t.table$MeJA.means[i] <- colMeans(MeJA.data)[i]
}
t.table$p.adj <- p.adjust(t.table$p.value, method = "fdr")
library(knitr)
t.table[t.table$p.value <= 0.2,] %>%
  kable(digits = 4)
