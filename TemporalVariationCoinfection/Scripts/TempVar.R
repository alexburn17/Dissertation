###########################################################################################
# Data Analysis for Temporal Variation Coinfection Study (Chapter 1)
# P. Alexander Burnham
# July 27, 2017
###########################################################################################

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/Dissertation/TemporalVariationCoinfection/Data")

# read in data:
TempVar <- read.csv("TempVarData.csv",header=TRUE,sep=",",stringsAsFactors=FALSE)

###############################################################################################
################################### PROGRAM BODY ##############################################
###############################################################################################

# source my functions
source("~/Dissertation/Scripts/BurnhamFunctions.R")
library(plyr)
library(ggplot2)
library(dplyr)
library(lme4)
library(car)
library(gridExtra)
library(grid)
library(cowplot)

###############################################################################################
# Variable scaling and manipulation:

# create binary variable for Nosema:
TempVar$NosemaBinary <- ifelse(TempVar$NosemaLoad == 0, 0, 1)

# create binary variable for Varroa:
TempVar$VarroaBinary <- ifelse(TempVar$VarroaLoad == 0, 0, 1)

# create scaled (0 to 1) variable for Brood Pattern:
TempVar$BroodPatternScaled <- TempVar$BroodPattern * 0.1 

# create scaled (0 to 1) variable for Frames of Brood:
TempVar$FOBnorm <- ((TempVar$FOB) - min(TempVar$FOB, na.rm = TRUE))/(max(TempVar$FOB, na.rm = TRUE)-min(TempVar$FOB, na.rm = TRUE))

###############################################################################################
# create temporal variation data frame:


###########################################################################
# function name: DFmaker
# description: creates a dataframe from sampling event variable and variable name
# parameters: 
# data = data frame
# Var1 = sampleing event
# Var2 = continous variable 
# name = name of varaible
# repNum = number of reps in each section for the data frame
# returns a dataframe with columns for sampling event, variable and variable name
###########################################################################

DFmaker <- function(data = TempVar,
                    Var1 = "SamplingEvent", 
                    Var2 = "DWVbinary", 
                    repNum = 80,
                    name = "DWV"){
  
  x <- cbind(select(data, Var1, Var2), rep(name, repNum))
  
  names(x)[2] <- "variable"
  names(x)[3] <- "variableName"
  
  return(x)
}

###########################################################################
# END OF FUNCITON
###########################################################################

# create data frames for each variable and merge them
DW <- DFmaker(name = "DWV", Var2 = "DWVbinary")

BQ <- DFmaker(name = "BQCV", Var2 = "BQCVbinary")

NO <- DFmaker(name = "Nosema", Var2 = "NosemaBinary")

VA <- DFmaker(name = "Varroa", Var2 = "VarroaBinary")

BP <- DFmaker(name = "Brood Pattern", Var2 = "BroodPatternScaled")

FB <- DFmaker(name = "Frames of Bees", Var2 = "FOBnorm")

# create TempDat dataframe:
TempDat <- rbind(DW, BQ, NO, VA, BP, FB)



###############################################################################################
# Plotting Temporal Variation data:

Temporal <- ddply(TempDat, c("variableName", "SamplingEvent"), summarise, 
                  n = length(variable),
                  mean = mean(variable, na.rm=TRUE),
                  sd = sd(variable, na.rm = TRUE),
                  se = sd / sqrt(n))

# split data frame into brood variables and pathogen variables
Brood <- Temporal[13:18,]
Path <- Temporal[1:12,]

# plotting pathogens through time 
p1 <- ggplot(data = Path, 
       aes(x = SamplingEvent, 
           y = mean, 
           color = variableName)
) + geom_point(size=4) + labs(x = NULL, y = "Prevalance") + coord_cartesian(ylim = c(0, 1), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(size=1) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.85, .24), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), axis.title.y=element_text(margin=margin(0,20,0,0))) + labs(color="Pathogen:") + scale_x_continuous(breaks=c(1,2,3))

# plotting brood measures through time
p2 <- ggplot( data = Brood, 
       aes(x = SamplingEvent, 
           y = mean, 
           group = variableName)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Rel. Intensity") + coord_cartesian(ylim = c(0, 1), xlim = c(1,3)) + geom_line(aes(linetype=variableName), size=1) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.22, .87), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Population:") + scale_x_continuous(breaks=c(1,2,3))


# use cowlpot package to combine the two figures:
plot_grid(p1, p2, align = "v", nrow = 2, rel_heights = c(5/8, 3/8))

###############################################################################################

Full <- lmer(data = TempVar, formula = logBQCV ~ VarroaLoad * SamplingEvent + (1|labID) + (SamplingEvent|Yard), REML=F)


summary(Full)

Null <- lmer(data = TempVar, formula = logBQCV ~ SamplingEvent + (1|labID) + (SamplingEvent|Yard), REML=F)

anova(Full, Null, test="LRT")



Full1 <- lmer(data = TempVar, formula = logDWV ~ VarroaBinary * SamplingEvent + (1|labID) + (SamplingEvent|Yard), REML=F)


Null1 <- lmer(data = TempVar, formula = logDWV ~ 1 * SamplingEvent + (1|labID) + (SamplingEvent|Yard), REML=F)



anova(Full1, Null1, test="LRT")


Full1 <- glmer(data = TempVar, formula = DWVbinary ~ VarroaLoad * SamplingEvent + (1|labID) + (SamplingEvent|Yard), family = binomial(link = "logit"))

Null1 <- glmer(data = TempVar, formula = DWVbinary ~ 1 * SamplingEvent + (1|labID) + (SamplingEvent|Yard), family = binomial(link = "logit"))

anova(Full1, Null1, test="LRT")





