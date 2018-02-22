# Nosema Data Analysis and Modeling
# P. Alexander Burnham
# March 26, 2016


#-------------------------------------------------------------------
#Preliminaries

# Clear memory of characters
ls()
rm(list=ls())
ls()

setwd("~/Desktop/RScripts/NosemaData_2")

# read in data
NosemaDF <- read.table("Nosema_Data_R.csv",header=TRUE,row.names=1,sep=",",stringsAsFactors=FALSE)

#attatch column names to column data
attach(NosemaDF)

#-------------------------------------------------------------------
#histogram of Nosema Count

hist(NosemaAverage,
     breaks=30, 
     col="steelblue",
     ylim=c(0,250),
     xlim=c(0,1300),
     xlab="Nosema Count per Bee",
     main = "Histogram of Nosema Count"
     )

library(RColorBrewer)
colors<-colorRampPalette(brewer.pal(9,"Blues"))(19)
colors<-rev(colors)

#Log transform data


LogNosema <- log(NosemaAverage)
hist(LogNosema, col=colors, ylim=c(0,20), xlim=c(-2,8))

#No0Nosema<-NosemaAverage[NosemaAverage !=0]

NosemaTable<-table(NosemaAverage)
barplot(height = NosemaTable, 
        ylim = c(0,250),
        ylab = "Frequency",
        xlab = "Counts",
        col=colors
)


#-------------------------------------------------------------------
# barplot showing percent prevalence by Species

SpeciesTab <- table(Species,NosemaPA)
SpeciesTab <- SpeciesTab[-3,]
SpeciesTab <- SpeciesTab[-3,]
SpeciesTab <- SpeciesTab[-5,]
Infected<-SpeciesTab[,2]
Uninfected<-SpeciesTab[,1]
x<-Infected/(Infected+Uninfected)
mean(x)
Species1 <- c("bimac.", "bor.", "imp.", "tern.", "vag.")

library(RColorBrewer)
colors<-colorRampPalette(brewer.pal(9,"Blues"))(5)

plot1 <- barplot(height=SpeciesTab[,2]/SpeciesTab[,1],
                 names.arg=Species1,
                 xlab="Bombus spp.",
                 ylab="Probability of Infection",
                 ylim=c(0,0.6),
                 #main=paste("Prevalence of Nosema by Bombus spp."),
                 col.main="blue",
                 font.lab=2,
                 las=1,
                 lwd=2,
                 col = colors,
                 border="black"
                 #add=TRUE
)

n <- SpeciesTab[,2]+SpeciesTab[,1]
p <- SpeciesTab[,2]/SpeciesTab[,1]
SE <- (p*(1-p))/n
SD=SE*sqrt(n)


# plot SD bars on plot:
arrows(plot1, p-SD*2, plot1, p+SD*2, lwd=2, angle=90, code=3)
#segments(plot1, p-SE*2, plot1, p+SE*2, lwd=2)
#arrows(plot1, p-SE*2, plot1, p+SE*2, lwd=2, angle=90, code=3)

#-----------------------------------------------------------------
#Contingecy table chi-sq

SpeciesTab <- table(Species,NosemaPA)
SpeciesTab
SpeciesTab <- SpeciesTab[-3,]
SpeciesTab <- SpeciesTab[-3,]
SpeciesTab <- SpeciesTab[-5,]
Infected<-SpeciesTab[,2]
Uninfected<-SpeciesTab[,1]
data <- rbind(Infected,Uninfected)
rownames(data) <- c("Infected", "Uninfected")
colnames(data) <- c("bimac.", "bor.", "imp.", "tern.", "vag.")
print(data)

colors<-colorRampPalette(brewer.pal(9,"Blues"))(5)

print(chisq.test(data))

as.matrix(data)

fisher.test(data)

mosaicplot(data,col=colors,shade=FALSE,
           main="Prevalence of Nosema by Species",
           ylab = "Nosema Spp.")


#-------------------------------------------------------------------
# barplot showing percent prevalence by Caste

colors1<-colorRampPalette(brewer.pal(9,"Blues"))(3)
CasteTab <- table(Caste,NosemaPA)
CasteTab
Caste1 <- c("Male", "Queen", "Worker")

plot1 <- barplot(height=CasteTab[,2]/CasteTab[,1],
                 names.arg=Caste1,
                 xlab="Caste",
                 ylab="Probability of Infection",
                 ylim=c(0,0.6),
                 #main=paste("Prevalence of Nosema by Caste"),
                 col.main="steelblue",
                 font.lab=2,
                 las=1,
                 lwd=2,
                 col = colors1,
                 border="black"
                 # add=TRUE
)

n <- CasteTab[,1]+CasteTab[,1]
p <- CasteTab[,2]/CasteTab[,1]
SE <- (p*(1-p))/n
SE  
SD=SE*sqrt(n)
SD

# plot SD bars on plot:
arrows(plot1, p-SD*2, plot1, p+SD*2, lwd=2, angle=90, code=3)
#segments(plot1, p-SE*2, plot1, p+SE*2, lwd=2)
#arrows(plot1, p-SE*2, plot1, p+SE*2, lwd=2, angle=90, code=3)

#-----------------------------------------------------------------
# Chisq on Caste

CasteTab <- table(Caste,NosemaPA)
CasteTab
Infected<-CasteTab[,2]
Uninfected<-CasteTab[,1]
data <- rbind(Infected,Uninfected)
rownames(data) <- c("Infected", "Uninfected")
colnames(data) <- c("Males", "Queens", "Workers")
print(data)

colors<-colorRampPalette(brewer.pal(9,"Blues"))(3)

print(chisq.test(data))
mosaicplot(data,col=colors,
           shade=FALSE,
           main="Prevalence of Nosema by Caste",
           ylab = "Caste")

#-------------------------------------------------------------------
# barplot showing percent prevalence by Site

colors2<-colorRampPalette(brewer.pal(9,"Blues"))(15)

SiteTab <- table(Site,NosemaPA)
SiteTab <-SiteTab[-16,]
SiteTab <-SiteTab[-16,]
Site1 <- c("100", "Amy", "CBF", "DanG", "GMC", "HORT", "LIB", "LOW", "MRW", "NORD", "OVC", "PC", "PEL","PVC", "Robin")

plot1 <- barplot(height=SiteTab[,2]/SiteTab[,1],
                 names.arg=Site1,
                 #xlab="Site",
                 ylab="Probability of Infection",
                 ylim=c(0,1),
                 #main=paste("Prevalence of Nosema by Site"),
                 col.main="steelblue",
                 font.lab=2,
                 las=2,
                 lwd=2,
                 col = colors2,
                 border="black"
                 #add=TRUE
)

#-------------------------------------------------------------------

#Contingecy table chi-sq site

SiteTab <- table(Site,NosemaPA)
SiteTab <-SiteTab[-16,]
SiteTab <-SiteTab[-16,]

Infected<-SiteTab[,2]
Uninfected<-SiteTab[,1]
data <- rbind(Infected,Uninfected)
rownames(data) <- c("Infected", "Uninfected")
colnames(data) <- c("100", "Amy", "CBF", "DanG", "GMC", "HORT", "LIB", "LOW", "MRW", "NORD", "OVC", "PC", "PEL","PVC", "Robin")
print(data)

colors<-colorRampPalette(brewer.pal(9,"Blues"))(15)

print(chisq.test(data))
mosaicplot(data,col=colors,shade=FALSE)




#-----------------------------------------------------------------
# Nosmea prevalence by hbPA

hbTab <- table(hbPA,NosemaPA)
hbTab
Infected<-hbTab[,2]
Uninfected<-hbTab[,1]
x<-Infected/(Infected+Uninfected)
mean(x)
HB1 <- c("Near", "Far")

library(RColorBrewer)
colors4<-colorRampPalette(brewer.pal(9,"Blues"))(2)

plot1 <- barplot(height=hbTab[,2]/hbTab[,1],
                 names.arg=HB1,
                 xlab="Proximity to Honeybee Colonies",
                 ylab="Probability of Infection",
                 ylim=c(0,0.6),
                 #main=paste("Prevalence of Nosema by Bombus spp."),
                 col.main="blue",
                 font.lab=2,
                 las=1,
                 lwd=2,
                 col = colors4,
                 border="black"
                 #add=TRUE
)

n <- hbTab[,2]+hbTab[,1]
p <- hbTab[,2]/hbTab[,1]
SE <- (p*(1-p))/n
SE  
SD=SE*sqrt(n)
SD

# plot SD bars on plot:
arrows(plot1, p-SD*2, plot1, p+SD*2, lwd=2, angle=90, code=3)

#-----------------------------------------------------------------
# Chisq on HB proximity
hbTab <- table(hbPA,NosemaPA)
hbTab
Infected<-hbTab[,2]
Uninfected<-hbTab[,1]
data <- rbind(Infected,Uninfected)
rownames(data) <- c("Infected", "Uninfected")
colnames(data) <- c("Near", "Far")
print(data)

colors<-colorRampPalette(brewer.pal(9,"Blues"))(2)

print(chisq.test(data))
mosaicplot(data,col=colors,shade=FALSE, main="Prevalence of Nosema by Honeybee Proximity",
           ylab = "Honeybee Proximity")





















#-----------------------------------------------------------------
# Nosema prob. of infection by HB proximity and species


x<-aggregate(NosemaPA~Date, data=NosemaDF, mean)
x 

barplot(height=x$NosemaPA)

NosemaAverage

x<-which(NosemaAverage == 0)
length(x)

#-----------------------------------------------------------------
#Mean nosema by species and caste

x<-aggregate(NosemaAverage~Species+Caste, data=NosemaDF, mean)
x<-data.frame(x)
x

#-----------------------------------------------------------------
# Nosmea count by date

dateplotmat <- cbind(NosemaDF$Date,NosemaDF$NosemaAverage)
head(dateplotmat)

dateplotmat <- dateplotmat[-221,] 
dateplotmat <- dateplotmat[-224,]
dateplotDF <- data.frame(dateplotmat)

plot(dateplotDF$X1,dateplotDF$X2)

LineBF <- lm(dateplotDF$X2~dateplotDF$X1, data=dateplotDF)
line<-abline(LineBF, col = "blue")
summary(LineBF)


#-----------------------------------------------------------------
# ANOSIM (vegan) non parametric statistical test HB

library(vegan)

anoHB <- data.frame(cbind(hbPA, NosemaAverage))
attach(anoHB) 
anoHB <- anoHB[-239,]
anoHB

z <- anosim(dat=anoHB,grouping=hbPA, distance="binomial",permutations=999)
summary(z) 
plot(z) 
z


#-----------------------------------------------------------------
# ANOSIM (vegan) non parametric statistical test caste

library(vegan)

anoCaste <- data.frame(cbind(Caste, NosemaAverage))
attach(anoCaste) 
Caste <- anoCaste[-239,]
anoCaste

z <- anosim(dat=anoCaste,grouping=Caste, distance="binomial",permutations=999)
summary(z) 
plot(z) 
z

#-----------------------------------------------------------------
# ADONIS (vegan) non parametric statistical test HB/Caste

library(vegan)

anoHB <- data.frame(cbind(factor(hbPA),factor(Caste),NosemaAverage))
attach(anoHB) 
anoHB <- anoHB[-239,]
anoHB

z<-adonis(formula=anoHB~V1+V2,data=anoHB, permutations=999)
summary(z)
plot(z)
z

#-----------------------------------------------------------------
# K-W non parmaetric test Caste

head(NosemaAverage)
head(Caste)
kwCaste <- data.frame(cbind(as.factor(Caste), NosemaAverage))
attach(kwCaste) 
kwCaste <- kwCaste[-239,]

z <- kruskal.test(x=NosemaAverage ,g=V1 ,data=kwCaste, formula=NosemaAverage~V1 )
summary(z)
z

kwCaste

boxplot(NosemaAverage~V1,data=kwCaste)

#-----------------------------------------------------------------
# K-W non parmaetric test HB

head(NosemaAverage)
head(hbPA)
kwHB <- data.frame(cbind(as.factor(hbPA), NosemaAverage))
attach(kwHB) 
kwHB <- kwHB[-239,]

z <- kruskal.test(x=NosemaAverage ,g=V1 ,data=kwHB, formula=NosemaAverage~V1 )
summary(z)
z

kwHB

#-----------------------------------------------------------------
# K-W non parmaetric test Species

head(NosemaAverage)
head(Species)
kwSP <- data.frame(cbind(Species, NosemaAverage))
attach(kwSP) 
kwSP <- kwSP[-61,]
kwSP <- kwSP[-167,]
kwSP <- kwSP[-175,]
kwSP <- kwSP[-204,]
kwSP <- kwSP[-235,]
kwSP
z <- kruskal.test(x=NosemaAverage ,g=Species ,data=kwSP, formula=NosemaAverage~Species )
summary(z)
z




#-----------------------------------------------------------------
# rank nosema data mean

rankN <- rank(NosemaAverage)
hbPA

rankDF<-data.frame(cbind(rankN,hbPA))

aggregate(rankN, by=list(hbPA), FUN=mean)


aggregate(formula=rankN~hbPA, data=rankDF, FUN = function(x) c(M=mean(x), SD=sd(x)))

library(dplyr)
rankDF<-data.frame(cbind(rankN,hbPA))
group_by(rankDF, hbPA) %>%
  summarise(mean=mean(rankN), sd=sd(rankN))
rankDF$rankN




#----------------------------------------------------------------
NosemaAverage>44



mastertable <- table(NosemaPA)
x <- data.frame(mastertable)
x
Infected <- x[2,2]
Uninfected <- x[1,2]
x<-Infected/(Infected+Uninfected)
Infected 
Uninfected
barplot(height=x,)

8/55






x <- NosemaDF[which(NosemaDF$NosemaAverage>30),]

x
