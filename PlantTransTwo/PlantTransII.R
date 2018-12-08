###################################################################
# Analysis of Plant Trans II
# P. Alexander Burnham
# August 23, 2018
###################################################################

#Preliminaries:
# Clear memoxry of characters
ls()
rm(list=ls())


# set working Directory:
setwd("~/Documents/GitHub/Dissertation/PlantTransTwo")

# Read in Data:
dat <- read.csv("PlantTransData.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE)



###########################################################################
# function name: PrelimClean
# description: removes unneed PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=MigVirus){
  
  # take only columns that we want:
  library(dplyr)
  data <- dplyr::select(data, ID, Cq_mean, target_name, quantity_mean, Plate)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="No Sample"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################

dat <- PrelimClean(data=dat)
table(dat$ID) 


# Read in Data:
dil <- read.csv("DilutionsPlantTrans.csv", 
                header=TRUE, 
                sep = ",", 
                stringsAsFactors = FALSE)


dat <- merge(x=dat, y=dil, by="ID")


###########################################################################
# function name: VirusNorm
# description: normalizes virus data with dilutions and constants 
# parameters: number of bees and a data frame 
# returns a dataframe with normalized virus values 
###########################################################################

VirusNorm <- function(number_bees = 1, data=data){
  
  # set constant values for genome copies per bee calculations:
  crude_extr <- 400
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





###########################################################################
# function name: actinNormal
# description: normalizes virus data with actin values 
# parameters: a data frame with actin values
# returns a dataframe with normalized virus values 
###########################################################################

actinNormal <- function(data=MigVirus){
  
  # pull only actin values out of dataframe
  ActinOnly <- data[which(data$target_name=="ACTIN"),]
  
  # create DF of ACTIN genome copies and lab ID:
  ActinDF <- data.frame(ActinOnly$ID, ActinOnly$genomeCopy)
  colnames(ActinDF) <- c("ID", "ACT_genomeCopy")
  
  # merge ACTIN dataframe with main dataframe:
  #Need rownames and all.x=TRUE because data frames are different sizes.
  data <- merge(data, ActinDF, by=c("ID"), all.x=TRUE)
  
  # find mean of all ACTIN values:
  ActinMean <- mean(ActinOnly$genomeCopy, na.rm = TRUE)
  
  # create column for normalized genome copies per bee:
  data$NormGenomeCopy <- (data$genomeCopy/data$ACT_genomeCopy)*ActinMean
  
  return(data)
}


###########################################################################
# END OF FUNCITON
###########################################################################





dat <- VirusNorm(number_bees = 1, data=dat)

dat <- actinNormal(data=dat)

dat <- dat[dat$target_name=="DWV",]



# Read in Data:
EcoDat <- read.csv("PlantTransII_Dats.csv", 
                header=TRUE, 
                sep = ",", 
                stringsAsFactors = FALSE)


dat <- merge(x=dat, y=EcoDat, by="ID")


write.csv(dat, "dat.csv")




















dat <- dat[dat$target_name=="DWV",]

dat1 <- dat[dat$exp==1,]

dat2 <- dat[dat$exp==2,]

plot(dat1$time, log10(dat1$NormGenomeCopy))


library(ggplot2)
library(plyr)



ggplot(dat1, aes(x=(time), y=log10(NormGenomeCopy+1))) + 
  geom_point(size=3) + theme_minimal(base_size = 17) + labs(x="Foraging Time (seconds)", y = "DWV Load log(genome copies)") + coord_cartesian(xlim = c(0, 100))
  

plot(dat2$Treatment, log10(dat2$NormGenomeCopy))


dat2$logCop <- log10(dat2$NormGenomeCopy + 1)
# summary stats for plotting purposes:
dose <- ddply(dat2, c("Treatment"), summarise, 
                        n = length(logCop),
                        mean = mean(logCop, na.rm = TRUE),
                        sd = sd(logCop, na.rm = TRUE),
                        se = sd / sqrt(n))

# remove NAs an reorder facotrs
dose <- dose[complete.cases(dose),]
dose$Treatment <- factor(dose$Treatment, levels = c(1,3,5,10))

library("RColorBrewer")

colors <- brewer.pal(n = 5, name = "PuBuGn")

#Create plot in ggplot 
ggplot(dose, aes(x=Treatment, y=mean, fill=Treatment)) + 
  geom_bar(stat="identity",  
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), 
                width=.4,
                position=position_dodge(.9)) + labs(x= "Dose (copies in milions)", y = "DWV Load log(genome copies)") + theme_classic(base_size = 17) + coord_cartesian(ylim = c(0, 8)) + theme(legend.position="none") + scale_fill_manual(values=colors[2:5])


















######################################################################################################
############################################# - FLOWERS - ############################################
######################################################################################################



# Read in Data:
dat3 <- read.csv("FlowersDat.csv", 
                header=TRUE, 
                sep = ",", 
                stringsAsFactors = FALSE)

dat3 <- dat3[dat3$Plate==40,]
dat3 <- dat3[dat3$target_name=="DWV",]

###########################################################################
# function name: PrelimClean
# description: removes unneed PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=MigVirus){
  
  # take only columns that we want:
  library(dplyr)
  data <- dplyr::select(data, ID, Cq_mean, target_name, quantity_mean, Plate)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="No Sample"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################

dat3 <- PrelimClean(data=dat3)
table(dat3$ID) 


# Read in Data:
dil <- read.csv("DilutionsPlantTrans.csv", 
                header=TRUE, 
                sep = ",", 
                stringsAsFactors = FALSE)


dat <- merge(x=dat, y=dil, by="ID")


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
  dil.factor <- 1
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  

  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################


dat3 <- VirusNorm(number_bees = 1, data=dat3)

dat3$PA <- c(1,1,1,1,1,1,1,0,0,0,0,1,0,0,1,1,1)
dat3$Experiment <- c(rep("HIT", 7), rep("Random", 7), rep("HBI", 3)) 

dat3$load <- ifelse(dat3$PA==1, dat3$genomeCopy, 0)
dat3$logLoad <- log10(dat3$load+1)











######################################################################################################
############################################# - FAKE FLOWERS - #######################################
######################################################################################################



# Read in Data:
dat4 <- read.csv("CottonFlower.csv", 
                 header=TRUE, 
                 sep = ",", 
                 stringsAsFactors = FALSE)


###########################################################################
# function name: PrelimClean
# description: removes unneed PCR file columns, gblock, NTC and duplicates 
# parameters: 
# data = data frame, 
# returns a DF that is partially cleaned
###########################################################################

PrelimClean <- function(data=MigVirus){
  
  # take only columns that we want:
  library(dplyr)
  data <- dplyr::select(data, ID, Cq_mean, target_name, quantity_mean, Match)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="No Sample"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################

dat4 <- PrelimClean(data=dat4)
table(dat4$ID) 





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
  dil.factor <- 1
  
  #create column for total_extr_vol
  total_extr_vol <- (GITCperbee * number_bees)
  
  
  # create column for genome copies per bee:
  data$genomeCopy <- ((((((data$quantity_mean / cDNA_eff) / rxn_vol) * dil.factor) * eluteRNA) / crude_extr) * total_extr_vol) / number_bees
  
  # norm_genomeCopy is 0 if NA
  data$genomeCopy[is.na(data$genomeCopy)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################


dat4 <- VirusNorm(number_bees = 1, data=dat4)


dat4$logLoad <- log10(dat4$genomeCopy+1)
Reverse$Match <- Reverse$ID

t <- merge(x=Reverse, y=dat4, by="Match")

write.csv(t, "Reverse.csv")


plot(t$logLoad.y, t$logLoad.x, ylim=c(2,6))
