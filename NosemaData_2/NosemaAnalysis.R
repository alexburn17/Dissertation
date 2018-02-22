#########################################################################################
# Nosema Data Analysis 
# P. Alexander Burnham
# November 7, 2017

#########################################################################################

# Clear memory of characters
ls()
rm(list=ls())

#Preliminaries
library(RColorBrewer)
library(plyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(car)
library(MASS)
library(scales)

# set working directory 
setwd("~/Dissertation/NosemaData_2")

# read in data
NosemaDF <- read.table("~/BurnhamAlexPrivate/BombusNosemaSurvey2014/Nosema_Data_R.csv",
                       header=TRUE,
                       sep=",",
                       stringsAsFactors=FALSE)

##########################################################################################

# new variable for log nosema:
NosemaDF$LogNosema <- log(NosemaDF$NosemaAverage + 1)
NosemaDF$LogSucrose <- log(NosemaDF$Sucrose + 1)

NosemaDF<-NosemaDF[!NosemaDF$Species==("Gris. "),]
NosemaDF<-NosemaDF[!NosemaDF$Species==("Terr."),]
NosemaDF<-NosemaDF[!NosemaDF$Species==("Ferv. "),]

##########################################################################################
# Nosema Prev by Species

# summary stats for plotting purposes:
VirusSummary <- ddply(NosemaDF, c("Species"), summarise, 
                      n = length(NosemaPA),
                      mean = mean(NosemaPA, na.rm = TRUE),
                      sd = sd(NosemaPA, na.rm = TRUE),
                      se = sd / sqrt(n))


# color pallette for graphics:
colors<-colorRampPalette(brewer.pal(9,"Blues"))(5)
colors<-rev(colors)

ggplot(VirusSummary, aes(x=Species, y=mean, fill=colors)) + geom_bar(stat="identity", col="black", position=position_dodge()) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.4, position=position_dodge(.9)) + scale_fill_manual(values=colors)  + labs(x="Species", y = "Nosema Prevalence") + theme_minimal(base_size = 17) + theme(legend.position=c(3, 3)) + coord_cartesian(ylim = c(0, .5)) + scale_y_continuous(labels = scales::percent)

##########################################################################################
# Nosema Prev by Caste

# summary stats for plotting purposes:
VirusSummary1 <- ddply(NosemaDF, c("Caste"), summarise, 
                      n = length(NosemaPA),
                      mean = mean(NosemaPA, na.rm = TRUE),
                      sd = sd(NosemaPA, na.rm = TRUE),
                      se = sd / sqrt(n))


# color pallette for graphics:
colors1<-colorRampPalette(brewer.pal(9,"Blues"))(3)
colors1<-rev(colors1)

# graphics for caste
ggplot(VirusSummary1, aes(x=Caste, y=mean, fill=colors1)) + geom_bar(stat="identity", col="black", position=position_dodge()) + geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.4, position=position_dodge(.9)) + scale_fill_manual(values=colors1)  + labs(x="Caste", y = "Nosema Prevalence") + theme_minimal(base_size = 17) + theme(legend.position=c(3, 3)) + coord_cartesian(ylim = c(0, .5)) + scale_y_continuous(labels = scales::percent)

##############################################
# statistical models:

# prevalence:
NosModPrev <- glmer(data=NosemaDF, NosemaPA~Caste+Species+Sucrose+hbPA + (1|Site), family = binomial(link = "logit"))
Anova(NosModPrev)

# load:
NosModLoad <- lmer(data=NosemaDF, LogNosema~Caste + Species + Sucrose + (1|Site))
Anova(NosModLoad)



# load of postive:
NosemaDFno0 <- NosemaDF[!NosemaDF$LogNosema==(0),]

NosModLoadno0  <- glmer(data=NosemaDFno0, LogNosema~Caste+Species+LogSucrose + (1|Site), family = Gamma)
Anova(NosModLoadno0)

r <- ddply(NosemaDFno0, c("Caste"), summarise, 
                       n = length(LogNosema),
                       mean = mean(LogNosema, na.rm = TRUE),
                       sd = sd(LogNosema, na.rm = TRUE),
                       se = sd / sqrt(n))


z <- ddply(NosemaDFno0, c("Species"), summarise, 
           n = length(LogNosema),
           mean = mean(LogNosema, na.rm = TRUE),
           sd = sd(LogNosema, na.rm = TRUE),
           se = sd / sqrt(n))



stepAIC(NosModLoadno0, direction = "both")







