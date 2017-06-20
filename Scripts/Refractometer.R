###########################################################################################
# P. Alexander Burnham
# May 25, 2017
# Test between two Refractometers 
###########################################################################################

# Preliminaries:
# Clear memory of characters:
ls()
rm(list=ls())

###########################################################################################
# write in data:

Trial <- rep(c(1:5),4)
Concentration <- c(rep("Ten", 10), rep("Twenty", 10))
Refractometer <- rep(c(rep("New", 5), rep("Old", 5)),2)
RefReading <- c(12.5,12.5,12.5,12.5,12.5,10,10,10,10,10,21,21,22,21,21,18,19,19,19,19)

data <- data.frame(Trial, Refractometer, Concentration, RefReading)

###########################################################################################
# data analysis 

splitDat <- split(data, data$Concentration)

mod20 <- aov(splitDat$Twenty$RefReading~splitDat$Twenty$Refractometer)
summary(mod20)

mod10 <- aov(splitDat$Ten$RefReading~splitDat$Ten$Refractometer)
summary(mod10)


splitDatConc <- split(data, data$Refractometer)

reading <- splitDatConc$New$RefReading - splitDatConc$Old$RefReading
conc <- c(rep("Ten",5), rep("Twenty",5))

df <- data.frame(conc, reading)

m <- aov(data=df, reading~conc)
summary(m)

mean(reading)
###########################################################################################
# data visualization 




