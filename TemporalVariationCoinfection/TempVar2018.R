###########################################################################################
# Temporal Variation and Coinfection Analysis Script 2018:
# Written by: P. Alexander Burnham
# July 28, 2018
###########################################################################################

# Preliminaries: 
# Clear memory of characters:
ls()
rm(list=ls())

# Set Working Directory 
setwd("~/Documents/GitHub/Dissertation/TemporalVariationCoinfection")

library("ggplot2")
library("dplyr")
library("lme4")
library("car")
library("plyr")

# load in data
TempVar <- read.csv("TempVar.csv", header=TRUE, stringsAsFactors=FALSE)







