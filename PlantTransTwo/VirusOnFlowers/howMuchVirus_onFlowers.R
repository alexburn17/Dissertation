###########################################################################################
# Figuring out how much virus to put on flowers for plant trans experiment part II
# P. Alexander Burnham
# May 16, 2018
##########################################################################################


# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
# setwd("~/Dissertation/PlantTransTwo/VirusOnFlowers") Sam Working Dirrectory 
setwd("~/Documents/GitHub/Dissertation/PlantTransTwo/VirusOnFlowers")

###########################################################################################
# Read in Virus Data:
PlantVir <- read.table("PlantqPCRResults.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

PlantDil <- read.table("RNAdilutions_NanodropResults.csv", header=TRUE, sep = ",", stringsAsFactors = FALSE) 

PlantEco <- read.csv("plantsEcoDat.csv", header=TRUE, stringsAsFactors = FALSE) 


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

# clean virus data set:
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
  GITCperbee <- 500
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



###########################################################################
# function name: CT_Threash
# description: creates binary data and makes genome copy 0 if below Ct threash
# parameters: dataframe
###########################################################################

CT_Threash <- function(data=data){
  
  splitDF <- split(data, data$target_name)
  
  # make DWV norm_genome_copbee 0 if Ct value is > 32.918
  splitDF$DWV$genomeCopy[which(splitDF$DWV$Ct_mean > 32.918)] <- 0
  splitDF$BQCV$genomeCopy[which(splitDF$BQCV$Ct_mean > 32.525)] <- 0
  
  splitDF$DWV$virusBINY  <- ifelse(splitDF$DWV$Ct_mean > 32.918, 0, 1)
  splitDF$BQCV$virusBINY  <- ifelse(splitDF$BQCV$Ct_mean > 32.525, 0, 1)

  # merge split dataframe back into "BombSurv" dataframe:
  data <- rbind(splitDF$DWV, splitDF$BQCV)
  
  # norm_genomeCopy is 0 if NA
  data$virusBINY[is.na(data$virusBINY)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################


# removed loads below limit of detection:
# dataFrame2 <- CT_Threash(dataFrame2)

dat <- dataFrame2

# Read in Virus Data:
#dat <- read.csv("CleanPlantVirus.csv", header=TRUE, stringsAsFactors = FALSE) 

# only positive values:
datPos <- dat[!dat$genomeCopy==0,]

# create data frames for each exp:
datSurv <- dat[dat$Exp=="Survey",]
datTrans <- dat[dat$Exp=="PlantTrans",]

# create data frames for each exp that are psotive:
datSurvPos <- datPos[datPos$Exp=="Survey",]
datTransPos <- datPos[datPos$Exp=="PlantTrans",]


library("plyr")
library("ggplot2")

#ddply summarize:
expPlot <- ddply(datPos, c("Exp"), summarise, 
              n = length(genomeCopy),
              mean = mean(genomeCopy, na.rm=TRUE),
              sd = sd(genomeCopy, na.rm=TRUE),
              se = sd / sqrt(n))

#choosing color pallet
colors <- c("blue", "grey")

#creating the figure of max gens:
plot1 <- ggplot(expPlot, aes(x=Exp, y=mean, fill=colors)) +
  geom_bar(stat="identity",
           position=position_dodge()) + labs(y="Viral Load (genome copies)", x="Experiment") + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.2),position=position_dodge(.9))

plot1 + theme_minimal(base_size = 17) + theme(legend.position=c(2, 2)) + scale_fill_manual(values=colors)



# print plant trans viral load data:
print(datTransPos)
#write.csv(datTransPos, "datTransPos.csv")
