###################################################################
# Cellualar Automata Model of Spillover HB -> Bombus
# P. Alexander Burnham
# June 22, 2018
###################################################################

#Preliminaries:
# Clear memoxry of characters
ls()
rm(list=ls())

# set working Directory:
setwd("~/Documents/GitHub/Dissertation/PlantTransTwo")

# Paramters:
#____________________________________________________________________
TimeSteps <- 50 # number of time steps

xDim <- 50 # x dimension of matrix
yDim <- 50 # y dimension of matrix

probHBbirth <- 80 # susceptable honey bee birth rate 
probBBbirth <- 80 # susceptable bumble bee birth rate
probHBdeath <- 1 # susceptable honey bee death rate
probBBdeath <- 1 # susceptable bumble bee death rate

probHBbirthINF <- 80 # infected honey bee birth rate 
probBBbirthINF <- 80 # infected bumble bee birth rate
probHBdeathINF <- 1 # infected susceptable honey bee death rate
probBBdeathINF <- 1 # infected susceptable bumble bee death rate

#____________________________________________________________________

require(RColorBrewer)

cols <- colorRampPalette(brewer.pal(5, "Blues"))(6)
colsINF <- colorRampPalette(brewer.pal(5, "Reds"))(6)
# vector of colors for each state:
cols <- c("white", cols[2:6], colsINF[2:6])


# initialize random matrix with starting proportions of HB, BB colonies, and flowers:
beeVec <-sample(0:10, xDim*yDim, replace=T, prob = c(0.8, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02))
beeMat <- matrix(data = beeVec, nrow = yDim, ncol = xDim)

# create a static copy of above matrix:
beeMatStatic <- beeMat

# function for creating file name with leading zeros
# makes it easier to process them sequentially
rename <- function(x){
  if (x < 10) {
    return(name <- paste('000',t,'plot.png',sep=''))
  }
  if (x < 100 && i >= 10) {
    return(name <- paste('00',t,'plot.png', sep=''))
  }
  if (x >= 100) {
    return(name <- paste('0', t,'plot.png', sep=''))
  }
}

# initialize a matrix to store number of each state in matrix at each time step
counter <- matrix(nrow=TimeSteps, ncol=11) 

# loop through time steps
for (t in 1:TimeSteps){
  
  # store table of counts in each row of counter matrix
  counter[t,] <- as.vector(table(beeMat))
  
  # loop through each cell in the matrix
  for (i in 1:xDim){
    for (j in 1:yDim){

        # create random unit movement steps for i and j
        temp1 <- sample(c(1, -1), 1 ,replace=T, prob = c(0.5, 0.5))
        temp2 <- sample(c(1, -1), 1 ,replace=T, prob = c(0.5, 0.5))

        # honeybee rules:  
        if(beeMat[i,j]==1){
          beeMat[i,j] <- 0
          if(runif(1, 1,100)>=probHBdeath){
          beeMat[ifelse(sum(i,temp1) %in% 1:xDim, 
                        sum(i,temp1), 
                        round(runif(1, 1, xDim))), 
                 ifelse(sum(j,temp2) %in% 1:yDim, 
                        sum(j,temp2), 
                        round(runif(1, 1, yDim)))] <- 1
           }
         }
 
        # bumble bee rules: 
        if(beeMat[i,j]==2){
          beeMat[i,j] <- 0 
          if(runif(1, 1,100)>=probBBdeath){
          beeMat[ifelse(sum(i,temp1) %in% 1:xDim, 
                       sum(i,temp1), 
                        round(runif(1, 1, xDim))), 
                ifelse(sum(j,temp2) %in% 1:yDim, 
                        sum(j,temp2), 
                        round(runif(1, 1, yDim)))] <- 2
           } 
        }

        # INFECTED honeybee rules:  
        if(beeMat[i,j]==6){
          beeMat[i,j] <- 0
          if(runif(1, 1,100)>=probHBdeathINF){
            beeMat[ifelse(sum(i,temp1) %in% 1:xDim, 
                          sum(i,temp1), 
                          round(runif(1, 1, xDim))), 
                   ifelse(sum(j,temp2) %in% 1:yDim, 
                          sum(j,temp2), 
                          round(runif(1, 1, yDim)))] <- 6
          }
        }

        # INFECTED bumble bee rules:  
        if(beeMat[i,j]==7){
          beeMat[i,j] <- 0
          if(runif(1, 1,100)>=probHBdeath){
            beeMat[ifelse(sum(i,temp1) %in% 1:xDim, 
                          sum(i,temp1), 
                          round(runif(1, 1, xDim))), 
                   ifelse(sum(j,temp2) %in% 1:yDim, 
                          sum(j,temp2), 
                          round(runif(1, 1, yDim)))] <- 7
          }
        }

        # ensure flowers and colonies remain static using static matrix
        if(beeMatStatic[i,j] == 3 || beeMatStatic[i,j] == 4 || 
           beeMatStatic[i,j] == 5 || beeMatStatic[i,j] == 8 ||
           beeMatStatic[i,j] == 9 || beeMatStatic[i,j] == 10){
          beeMat[i,j] <- beeMatStatic[i,j]
        }

        # honeybee colony rules
        if(beeMatStatic[i,j] == 3){
          if(runif(1, 1,100)<=probHBbirth){
          beeMat[ifelse(sum(i,temp1) %in% 1:xDim, 
                        sum(i,temp1), 
                        round(runif(1, 1, xDim))), 
                 ifelse(sum(j,temp2) %in% 1:yDim, 
                        sum(j,temp2), 
                        round(runif(1, 1, yDim)))] <- 1
          }
        }

        # bumble bee colony rules
        if(beeMatStatic[i,j] == 4){
          if(runif(1, 1,100)<=probBBbirth){
            beeMat[ifelse(sum(i,temp1) %in% 1:xDim, 
                          sum(i,temp1), 
                          round(runif(1, 1, xDim))), 
                   ifelse(sum(j,temp2) %in% 1:yDim, 
                          sum(j,temp2), 
                          round(runif(1, 1, yDim)))] <- 2
          }
        }
        

        # INFECTED honeybee colony rules
        if(beeMatStatic[i,j] == 9){
          if(runif(1, 1,100)<=probHBbirthINF){
            beeMat[ifelse(sum(i,temp1) %in% 1:xDim, 
                          sum(i,temp1), 
                          round(runif(1, 1, xDim))), 
                   ifelse(sum(j,temp2) %in% 1:yDim, 
                          sum(j,temp2), 
                          round(runif(1, 1, yDim)))] <- 6
          }
        }

        # INFECTED bumble bee colony rules
        if(beeMatStatic[i,j] == 10){
          if(runif(1, 1,100)<=probBBbirthINF){
            beeMat[ifelse(sum(i,temp1) %in% 1:xDim, 
                          sum(i,temp1), 
                          round(runif(1, 1, xDim))), 
                   ifelse(sum(j,temp2) %in% 1:yDim, 
                          sum(j,temp2), 
                          round(runif(1, 1, yDim)))] <- 7
          }
        }

    } # end j loop
} # end i loop

# create names for each png:
name <- rename(t)
png(name)
 
# create each image to visulaize as a matrix 
image(1:nrow(beeMat), 1:ncol(beeMat), as.matrix(beeMat), col=cols, asp=1, xaxt='n', yaxt='n', ann=FALSE, bty='n')

dev.off()

} # end of time step for loop

#run ImageMagick: creates a gif of all images
my_command <- 'convert *.png -delay 3 -loop 0 animation.gif'
system(my_command)































