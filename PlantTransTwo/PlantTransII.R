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
dat <- read.csv("PlantTransII.csv", 
                    header=TRUE, 
                    sep = ",", 
                    stringsAsFactors = FALSE)

# Read in Data:
dat2 <- read.csv("PlantTransII_Exp2.csv", 
                header=TRUE, 
                sep = ",", 
                stringsAsFactors = FALSE)

dat <- merge(dat, dat2, by.x = c("ID"), all.x=TRUE)





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
  data <- dplyr::select(data, ID, Cq_mean, target_name, quantity_mean, dil.factor, Colony, Treatment, time, exp)
  
  # remove duplicate rows
  data<-data[!duplicated(data), ]
  
  # remove NTC rows from dataframe:
  data<-data[!(data$ID=="No Sample"),]
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################



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





###########################################################################
# function name: CT_Threash
# description: creates binary data and makes genome copy 0 if below Ct threash
# parameters: dataframe
###########################################################################

CT_Threash <- function(data=data){
  
  splitDF <- split(data, data$target_name)
  
  # make DWV norm_genome_copbee 0 if Ct value is > 32.918
  splitDF$DWV$NormGenomeCopy[which(splitDF$DWV$Cq_mean > 32.918)] <- 0
  splitDF$BQCV$NormGenomeCopy[which(splitDF$BQCV$Cq_mean > 32.525)] <- 0
  
  splitDF$DWV$virusBINY  <- ifelse(splitDF$DWV$Cq_mean > 32.918, 0, 1)
  splitDF$BQCV$virusBINY  <- ifelse(splitDF$BQCV$Cq_mean > 32.525, 0, 1)
  
  # merge split dataframe back into "BombSurv" dataframe:
  data <- rbind(splitDF$DWV, splitDF$BQCV, splitDF$IAPV)
  
  # norm_genomeCopy is 0 if NA
  data$virusBINY[is.na(data$virusBINY)] <- 0
  
  return(data)
  
}

###########################################################################
# END OF FUNCITON
###########################################################################


dat <- PrelimClean(data=dat)

dat <- VirusNorm(number_bees = 1, data=dat)

dat <- dat[which(dat$ID %in% c(3:200)),]

dat <- actinNormal(data=dat)




dat <- dat[dat$target_name=="DWV",]

dat1 <- dat[dat$exp==1,]

dat2 <- dat[dat$exp==2,]

plot(dat1$time, log10(dat1$NormGenomeCopy))


library(ggplot2)
library(plyr)



ggplot(dat1, aes(x=log10(time), y=log10(NormGenomeCopy+1))) + 
  geom_point(size=3) + theme_minimal(base_size = 17) + labs(x="Foraging Time log(seconds)", y = "DWV Load log(genome copies)") 
  

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

