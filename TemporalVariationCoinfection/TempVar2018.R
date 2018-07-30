###########################################################################################
# Temporal Variation and Coinfection Analysis Script 2018:
# Written by: P. Alexander Burnham
# July 28, 2018
###########################################################################################

# Preliminaries: 

# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory:
setwd("~/Documents/GitHub/Dissertation/TemporalVariationCoinfection")

# load packages and functions:
source("BurnhamFunctionsTempVar.R")
library("ggplot2")
library("dplyr")
library("lme4")
library("car")
library("plyr")

# load in data
TempVar <- read.csv("TempVar.csv", header=TRUE, stringsAsFactors=FALSE)
dilutions <- read.csv("dilutions_TempVar.csv", header=TRUE, stringsAsFactors=FALSE)
Temp <- read.csv("TemporaryData.csv", header=TRUE, stringsAsFactors=FALSE)
FieldDat <- read.csv("tempVarFieldDat.csv", header=TRUE, stringsAsFactors=FALSE)

# merge together Temp and Field Data
EcoDat <- merge(Temp, FieldDat, by = c("Date", "Site", "FieldID"), all.x = TRUE)

####################################################################################################
# Program Body:

# preliminary cleaning -> removes duplicate rows and control data
TempVarClean <- PrelimClean(data = TempVar)

# merge data sets to inlude dilution data for Normalization:
TempVarClean <- merge(TempVarClean, dilutions, by = "ID", all.x = TRUE)

# normalize data set Viral Load
TempVarClean <- VirusNorm(number_bees = 1, data = TempVarClean)

# normalize viral laod by actin
TempVarClean <- actinNormal(data = TempVarClean)

# make binary variable and use threashld of ct for limit of detection: 
TempVarClean <- CT_Threash(data = TempVarClean)

# merge eco data by lab data by ID
TempVarClean <- merge(TempVarClean, EcoDat, by = "ID", all.x = TRUE)

# log base 10 transform genome copies
TempVarClean$logVL <- log10(TempVarClean$NormGenomeCopy + 1)

# split data set DWV and BQCV
TempVarCleanSplit <- split(TempVarClean, TempVarClean$target_name)
BQCV <- TempVarCleanSplit$BQCV
DWV <- TempVarCleanSplit$DWV

# create data set that has binary and log viral loads for each variable 
DWVtemp <- select(DWV, ID, virusBINY, logVL)
names(DWVtemp) <- c("ID", "binaryDWV", "logDWV")
modelDat <- merge(BQCV, DWVtemp, by = "ID")


mod <- lmer(data=modelDat, formula = logVL ~ binaryDWV * sampling_event + Caste + (1|Site) + (1|Host_Plant))
Anova(mod)



# family = binomial(link = "logit")