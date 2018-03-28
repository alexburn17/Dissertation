###########################################################################################
# Data Analysis for Temporal Variation Coinfection Study (Chapter 1)
# P. Alexander Burnham
# July 27, 2017
###########################################################################################

# Preliminaries:
# Clear memory of characters:
rm(list=ls())

# Set Working Directory: 
setwd("~/BurnhamAlexPrivate/TempVar_Coinfection_Data")

# read in data:
TempVar <- read.csv("TempVarData.csv",header=TRUE,sep=",",stringsAsFactors=FALSE, comment.char = '#')

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

# make sampling event a factor:
TempVar$SamplingEventFact <- as.character(TempVar$SamplingEvent)

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


###########################################################################
# function name: TempVarFunc (requires plyr and ggplot2)
# description: creates summary stat DFs and plots rough figures plotting a
# continuous or binary varable by varroa presence of absence 
# parameters: 
# name = character name of variable within data
# data = data frame, 
# returns a plot and a dataframe
###########################################################################

TempVarFunc <- function(name="DWVbinary", data = TempVar){
  
  # prepare character for use in dplyr environment
  data$x <- data[,name]
  
  # creating a data frame with means, sd and se for each variable
  Temp <- ddply(data, c("VarroaBinary", "SamplingEvent"), summarise, 
                n = length(x),
                mean = mean(x, na.rm=TRUE),
                sd = sd(x, na.rm = TRUE),
                se = sd / sqrt(n))
  
  Temp<-Temp[!(Temp$SamplingEvent==3),]
  
  # create a character string of infected with varroa or not
  Temp$Var <- ifelse(Temp$VarroaBinary==1, "Infected", "Uninfected")
  
# plot data frame
p <- ggplot( data = Temp, 
        aes(x = SamplingEvent, 
            y = mean, 
            color = Var)
) + geom_point(size=4.5) + labs(x = "Sampling Event", y = name) + coord_cartesian(xlim = c(1,2)) + geom_line(size=1.5) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + scale_color_manual(values=c("black", "blue")) + theme_classic(base_size = 19) + theme(legend.position=c(.22, .87)) + labs(color="Varroa Status:") + scale_x_continuous(breaks=c(1,2))

z <- list(Temp, p)

return(z)

}
###########################################################################
# END OF FUNCITON
###########################################################################


###############################################################################################
######## PROGRAM BODY #########################################################################
###############################################################################################

###############################################################################################
# Plotting Temporal Variation PREVALENCE data:

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


# creating seperate plots for prevalence and population metrics:
p1sep <- ggplot(data = Path, 
                aes(x = SamplingEvent, 
                    y = mean, 
                    color = variableName)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Prevalance") + coord_cartesian(ylim = c(0, 1), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.05)) + geom_line(size=1) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.85, .24)) + labs(color="Pathogen:") + scale_x_continuous(breaks=c(1,2,3))

print(p1sep)



# plotting brood measures through time
p2sep <- ggplot( data = Brood, 
                 aes(x = SamplingEvent, 
                     y = mean, 
                     group = variableName)
) + geom_point(size=4) + labs(x = "Sampling Event", y = "Rel. Intensity") + coord_cartesian(ylim = c(0, 1), xlim = c(1,3)) + geom_line(aes(linetype=variableName), size=1) + scale_fill_brewer(palette = "Paired") + theme_classic(base_size = 17) + theme(legend.position=c(.22, .87), panel.border = element_blank(), axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'), axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(linetype="Population:") + scale_x_continuous(breaks=c(1,2,3)) 

print(p2sep)

###############################################################################################


# create list of all variable names
namesList <- list("NosemaBinary", "VarroaBinary", "DWVbinary", "BQCVbinary", "VarroaLoad", "NosemaLoad", "logBQCV", "logDWV", "FOB", "BroodPattern")

# apply function to all variable names
resultList <- lapply(namesList, TempVarFunc, data=TempVar)


###############################################################################################
library("MuMIn")

# model of loBQCV by varroa load using a gamma distribution:
Full <- lmer(data = TempVar, formula = logBQCV ~ VarroaLoad * SamplingEvent + (1|labID) + (SamplingEvent|Yard))

# running an anova on the model for determining significance:
Anova(Full)
r.squaredGLMM(Full)




# main plot parameters
p <- ggplot(TempVar,aes(VarroaLoad, logBQCV)) + geom_point(aes(color=SamplingEventFact), size = 3.5) + coord_cartesian(xlim = c(0,25), ylim=c(10,25))  

# aestetic parameters 
p + theme_classic(base_size = 17) + theme(legend.position=c(.85, .8), legend.background = element_rect(color = "black", size = .5)) + labs(color="Time Point:", y="BQCV log(genome copies/bee)", x="Varroa Load (mites/100 bees)") + scale_color_manual(values = c("slategrey", "blue", "black")) + geom_smooth(aes(color=SamplingEventFact), method  = lm, se = FALSE)


###############################################################################################
# Same plot as above but without time point groups:

# main plot parameters
p1 <- ggplot(TempVar,aes(VarroaLoad, logBQCV)) + geom_point(size = 3.5) + coord_cartesian(xlim = c(0,25), ylim=c(10,25))  

# aestetic parameters 
p1 + theme_classic(base_size = 17) + theme(legend.position=c(.17, .8), legend.background = element_rect(color = "black", size = .5)) + labs(y="BQCV log(genome copies/bee)", x="Varroa Load (mites/100 bees)") + geom_smooth(method  = lm, se = FALSE)

# model of loBQCV by varroa load using a gamma distribution:
Full1 <- lmer(data = TempVar, formula = logBQCV ~ VarroaLoad + (1|labID) + (SamplingEvent|Yard))

# running an anova on the model for determining significance:
Anova(Full1)

r.squaredGLMM(Full1)


###############################################################################################


Full1 <- lmer(data = TempVar, formula = FOB ~  (DWVbinary + VarroaBinary) * SamplingEvent + (1|labID) + (1|Yard))

Anova(Full1)


###############################################################################################

# create pathogen richness varaible 
TempVar$PathRich <- TempVar$NosemaBinary + TempVar$VarroaBinary + TempVar$DWVbinary + TempVar$BQCVbinary + TempVar$ChalkBrood + TempVar$AFB + TempVar$EFB + TempVar$PMS + TempVar$SBV + TempVar$Snot + TempVar$BSB

# run anova and post hoc test on Pathogen Richness variable
x <- aov(TempVar$PathRich~as.factor(TempVar$SamplingEvent))
summary(x)
TukeyHSD(x)

# pathoengen richness by sampling Event
kruskal.test(TempVar$PathRich~as.factor(TempVar$SamplingEvent))


#ddply summarize:
Temp <- ddply(TempVar, c("SamplingEvent"), summarise, 
                        n = length(PathRich),
                        mean = mean(PathRich, na.rm=TRUE),
                        sd = sd(PathRich, na.rm=TRUE),
                        se = sd / sqrt(n))


#creating the figure
#choosing color pallet
colors <- c("steelblue", "blue", "grey")

plot1 <- ggplot(Temp, aes(x=SamplingEvent, y=mean, fill=colors)) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(y="Pathogen Richness", x="Sampling Event") + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.2),position=position_dodge(.9))

plot1 + theme_bw(base_size = 17) + scale_fill_manual(values=colors) + theme(legend.position="none") + coord_cartesian(ylim = c(0, 6)) + annotate(geom = "text", x = 1, y = 3.5, label = "A",cex = 6) + annotate(geom = "text", x = 2, y = 3.6, label = "A",cex = 6) + annotate(geom = "text", x = 3, y = 4.8, label = "B",cex = 6)






###############################################################################################

Temp1 <- ddply(TempVar, c("PathRich"), summarise, 
               n = length(VarroaLoad),
               mean = mean(VarroaLoad, na.rm=TRUE),
               sd = sd(VarroaLoad, na.rm=TRUE),
               se = sd / sqrt(n))

# edit temporary data frame to remove low sample size points
Temp1$PathRich <- as.factor(Temp1$PathRich)
Temp1 <- Temp1[-1,]
Temp1 <- Temp1[-4,]

# change levels to low medium and high
levels(Temp1$PathRich) <- c("A","Low","Medium","High","F")

#creating the figure
#choosing color pallet
colors <- c("darkolivegreen4", "darkgreen", "grey")

plot2 <- ggplot(Temp1, aes(x=PathRich, y=mean, fill=rev(colors))) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(y="Varroa Load", x="Pathogen Richness") + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.2),position=position_dodge(.9))

plot2 + theme_bw(base_size = 17) + scale_fill_manual(values=colors) + theme(legend.position="none") + coord_cartesian(ylim = c(0, 7)) + annotate(geom = "text", x = 1, y = 1.8, label = "A",cex = 6) + annotate(geom = "text", x = 2, y = 4.1, label = "B",cex = 6) + annotate(geom = "text", x = 3, y = 6.7, label = "B",cex = 6)  + coord_flip()



# data analysis for var laod vs path rich
threeVal <- TempVar[!TempVar$PathRich > 4,] 
threeVal <- threeVal[!TempVar$PathRich < 2,] 
z <- aov(threeVal$VarroaLoad~as.factor(threeVal$PathRich))
summary(z)
TukeyHSD(z)

kruskal.test(threeVal$VarroaLoad~as.factor(threeVal$PathRich))



fisher.test(TempVar$VarroaBinary, as.factor(TempVar$PathRich))

Temp2 <- ddply(TempVar, c("PathRich"), summarise, 
              n = length(VarroaBinary),
              mean = mean(VarroaBinary, na.rm=TRUE),
              sd = sd(VarroaBinary, na.rm=TRUE),
              se = sd / sqrt(n))

# change levels to low medium and high
Temp2$PathRich <- as.factor(Temp2$PathRich)
levels(Temp2$PathRich) <- c("A","Low","Medium","High","F")
Temp2 <- Temp2[-1,]
Temp2 <- Temp2[-4,]



#creating the figure
#choosing color pallet
colors <- c("slategrey", "darkorange4", "grey")

plot3 <- ggplot(Temp2, aes(x=PathRich, y=mean, fill=rev(colors))) +
  geom_bar(stat="identity", color="black",
           position=position_dodge()) + labs(y="Varroa Percent Prevalence ", x="Pathogen Richness")

plot3 + theme_bw(base_size = 17) + scale_fill_manual(values=colors) + theme(legend.position="none") + coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent) + coord_flip()




# creat ggplot stock colors
gg_color_hue <- function(n) {
  hues <- seq(15,375,length=n+1)
  hcl(h=hues,l=65,c=100)[1:n]
}

  n <- 4
cols <- gg_color_hue(n)









VarLoad <- ddply(TempVar, c("SamplingEvent"), summarise, 
               n = length(VarroaLoad),
               mean = mean(VarroaLoad, na.rm=TRUE),
               sd = sd(VarroaLoad, na.rm=TRUE),
               se = sd / sqrt(n))

NosLoad <- ddply(TempVar, c("SamplingEvent"), summarise, 
            n = length(NosemaLoad),
            mean = mean(NosemaLoad, na.rm=TRUE),
            sd = sd(NosemaLoad, na.rm=TRUE),
            se = sd / sqrt(n))

# create cirus column in data frame
ParaDat <- rbind(VarLoad,NosLoad)
Parasite <- c(rep("Varroa",3), rep("Nosema",3))
ParaDat <- cbind(ParaDat, Parasite)

# colors for parasites
cols[3]
cols[4]


ggplot(data = ParaDat, 
       aes(x = SamplingEvent, 
           y = mean, 
           color = Parasite)
) + geom_point(size=4.5) + labs(x = "Sampling Event", y = "Parasite Load") + coord_cartesian(ylim = c(0, 17), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.1)) + geom_line(size=1.7) + scale_fill_brewer(palette = "Paired") + theme_bw(base_size = 17) + theme(legend.position=c(.85, .85)) + labs(linetype="Parasite") + scale_color_manual(values = c(cols[3], cols[4])) + scale_x_continuous(breaks=c(1,2,3))















BQCVLoad <- ddply(TempVar, c("SamplingEvent"), summarise, 
                 n = length(logBQCV),
                 mean = mean(logBQCV, na.rm=TRUE),
                 sd = sd(logBQCV, na.rm=TRUE),
                 se = sd / sqrt(n))


DWVLoad <- ddply(TempVar, c("SamplingEvent"), summarise, 
                  n = length(logDWV),
                  mean = mean(logDWV, na.rm=TRUE),
                  sd = sd(logDWV, na.rm=TRUE),
                  se = sd / sqrt(n))


# create cirus column in data frame
VirusDat <- rbind(DWVLoad, BQCVLoad)
Virus <- c(rep("DWV",3), rep("BQCV",3))
VirusDat <- cbind(VirusDat, Virus)

# colors for viruses
cols[1]
cols[2]


ggplot(data = VirusDat, 
       aes(x = SamplingEvent, 
           y = mean, 
           color = Virus)
) + geom_point(size=4.5) + labs(x = "Sampling Event", y = "Viral Load log(copies/bee)") + coord_cartesian(ylim = c(0, 25), xlim = c(1,3)) + geom_errorbar(aes(ymin = mean - se, ymax = mean + se, width = 0.1)) + geom_line(size=1.7) + scale_fill_brewer(palette = "Paired") + theme_bw(base_size = 17) + theme(legend.position=c(.2, .85)) + labs(linetype="Virus") + scale_color_manual(values = c(cols[2], cols[1])) + scale_x_continuous(breaks=c(1,2,3))

