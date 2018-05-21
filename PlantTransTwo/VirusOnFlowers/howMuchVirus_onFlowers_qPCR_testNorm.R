###########################################################################################
# Figuring out how much virus to put on flowers for plant trans experiment part II
# P. Alexander Burnham
# May 16, 2018
##########################################################################################


# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/Dissertation/PlantTransTwo/VirusOnFlowers")



##############################################################################
##############################################################################
# PLATE 1 ####################################################################
##############################################################################
##############################################################################


###########################################################################################
# Read in Virus Data:
PlantVirTest <- read.table("20180518_PlantTrans2_HowMuchVirusOnPlants_Run11.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

###########################################################################
# function name: PrelimClean
# description: removes unneeded PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=data){
  
  # take only columns that we want:
  library(dplyr)
  data <- dplyr::select(data, ID, target_name, Ct_mean, Ct_sd, quantity_mean)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="No Sample"),]
  
  # remove Gblock rows from dataframe:
  data<-data[!(data$ID=="Gblock"),]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="NTC"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################

# cean virus data set:
PlantVirTest <- PrelimClean(data=PlantVirTest)
PlantVirTest$dil.factor <- rep(1, 5)


###########################################################################
# function name: VirusNorm
# description: normalizes virus data with dilutions and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 1, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 200
  eluteRNA <- 50
  GITCperbee <- 200
  cDNA_eff <- 0.1
  rxn_vol <- 3
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  
  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################


PlantVirTest <- VirusNorm(number_bees = 1, PlantVirTest)

format(17137954170, scientific=TRUE) # stock solution Raw
format(44714450, scientific=TRUE) # stock solution after RNA extraction







##############################################################################
##############################################################################
# PLATE 2 ####################################################################
##############################################################################
##############################################################################







# Read in Virus Data:
PlantVirTest2 <- read.table("20180519_PlantTrans2_HowMuchVirusOnPlants_Run12.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

###########################################################################
# function name: PrelimClean
# description: removes unneeded PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=data){
  
  # take only columns that we want:
  library(dplyr)
  data <- dplyr::select(data, ID, target_name, Ct_mean, Ct_sd, quantity_mean)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="No Sample"),]
  
  # remove Gblock rows from dataframe:
  data<-data[!(data$ID=="G-Block"),]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="NTC"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################

# cean virus data set:
PlantVirTest2 <- PrelimClean(data=PlantVirTest2)
PlantVirTest2$dil.factor <- c(1, 1, 7.7, 3, 1, 1, 3.4, 4.4)


###########################################################################
# function name: VirusNorm
# description: normalizes virus data with dilutions and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 1, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 350
  eluteRNA <- 50
  GITCperbee <- 600
  cDNA_eff <- 0.1
  rxn_vol <- 3
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  
  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################

PlantVirTest2$quantity_mean <- as.numeric(PlantVirTest2$quantity_mean)
PlantVirTest2 <- VirusNorm(number_bees = 1, PlantVirTest2)
















##############################################################################
##############################################################################
# PLATE Single BQCV RUN ######################################################
##############################################################################
##############################################################################


# Read in Virus Data:
PlantVirTest3 <- read.table("SingleFlowerBQCVResultsPos.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

###########################################################################
# function name: PrelimClean
# description: removes unneeded PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=data){
  
  # take only columns that we want:
  library(dplyr)
  data <- dplyr::select(data, ID, target_name, Ct_mean, Ct_sd, quantity_mean, dil.factor)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="No Sample"),]
  
  # remove Gblock rows from dataframe:
  data<-data[!(data$ID=="G-Block"),]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="NTC"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################

# cean virus data set:
PlantVirTest3 <- PrelimClean(data=PlantVirTest3)



###########################################################################
# function name: VirusNorm
# description: normalizes virus data with dilutions and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 1, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 350
  eluteRNA <- 50
  GITCperbee <- 600
  cDNA_eff <- 0.1
  rxn_vol <- 3
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  
  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################

PlantVirTest3 <- VirusNorm(number_bees = 1, PlantVirTest3)







##############################################################################
##############################################################################
# PLATE 3 ####################################################################
##############################################################################
##############################################################################








# Read in Virus Data:
PlantVirTest4 <- read.table("DataFromVirusPlants_qPCR_plate3.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

###########################################################################
# function name: PrelimClean
# description: removes unneeded PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=data){
  
  # take only columns that we want:
  library(dplyr)
  data <- dplyr::select(data, ID, target_name, Ct_mean, Ct_sd, quantity_mean)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="No Sample"),]
  
  # remove Gblock rows from dataframe:
  data<-data[!(data$ID=="G-Block"),]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="NTC"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################

# cean virus data set:
PlantVirTest4 <- PrelimClean(data=PlantVirTest4)
PlantVirTest4$dil.factor <- c(15, 15.4, 14, 15.2)


###########################################################################
# function name: VirusNorm
# description: normalizes virus data with dilutions and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 1, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 350
  eluteRNA <- 50
  GITCperbee <- 600
  cDNA_eff <- 0.1
  rxn_vol <- 3
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  
  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * data$dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################


PlantVirTest4 <- VirusNorm(number_bees = 1, PlantVirTest4)

