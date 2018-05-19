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

###########################################################################################
# Read in Virus Data:
PlantVir <- read.table("PlantqPCRResults.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

###########################################################################
# function name: PrelimClean
# description: removes unneeded PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=ImidVirus){
  
  # take only columns that we want:
  library(dplyr)
  data <- select(data, ID, target_name, Ct_mean, Ct_sd, quantity_mean)
  
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
PlantVir <- PrelimClean(PlantVir)

# merge eco data and dilution factors
dataFrame <- merge(PlantDil, PlantEco, by = "ID", all.x = TRUE)

# merge virus and eco and dilution data:
dataFrame2 <- merge(PlantVir, dataFrame, by = "ID")

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


dataFrame2 <- VirusNorm(number_bees = 1, dataFrame2)






