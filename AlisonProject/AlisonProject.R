###########################################################################################
# Data Analysis for Pollen Study
# Alison Brody and P. Alexander Burnham
# June 1, 2017
###########################################################################################

# Preliminaries:
# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory: 
setwd("~/Desktop/AlisonProject")

# Read in Data:
data <- read.table("Pollen.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

###########################################################################################

# subset data to just include controls for question 1
dataControl<-data[data$StalkTrt==("C"),]
# removing black (not enough in data set)
dataControl<-dataControl[!dataControl$Color==("black"),]
dataControl<-dataControl[!dataControl$Color==("pink"),]
dataControl<-dataControl[!dataControl$Color==("red"),]


dataBoth <- data[!data$Color==("black"),]
dataBoth <- dataBoth[!dataBoth$Color==("pink"),]
dataBoth <- dataBoth[!dataBoth$Color==("red"),]

###########################################################################################
library(lme4)
library(car)

#DWV prevalence
Fullmod <- lmer(data=dataControl, formula = SumSeedWeight ~ PollenRemoval * Sex + (1|Hood/Color))

Fullmod <- lmer(data=dataControl, formula = ArcsinePerset ~ (PollenRemoval * Color) * Sex + (1|Hood))

summary(Fullmod)
Anova(Fullmod)




Fullmod1 <- lmer(data=dataBoth, formula = SumSeedWeight ~ StalkTrt + (1|Hood) + (1|Color))

summary(Fullmod)
Anova(Fullmod1)


###########################################################################################
# data visualization
library(plyr)
library(ggplot2)

sum <- ddply(dataControl, c("Color", "PollenRemoval", "Sex"), summarise, 
      n = length(ArcsinePerset),
      mean = mean(ArcsinePerset, na.rm=TRUE),
      sd = sd(ArcsinePerset, na.rm = TRUE),
      se = sd / sqrt(n))
sum <- sum[-13,]

colors <- c("slategray3", "dodgerblue4")

plot1 <- ggplot(sum, aes(x=Color, y=mean, fill=PollenRemoval)) + 
  geom_bar(stat="identity", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x="Week", y = "Arcsine of % Fruit Set") + theme_minimal(base_size = 17) + coord_cartesian(ylim = c(0, 1.5)) + scale_fill_manual(values=colors, 
     name="Pollen") + facet_wrap( ~ Sex)
plot1

