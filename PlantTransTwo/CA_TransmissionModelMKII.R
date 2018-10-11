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
#----------------------------------------------------------------------------
TimeSteps <- 300 # number of time steps

xDim <- 50 # x dimension of matrix
yDim <- 50 # y dimension of matrix

probBirth <- 5 # bee birth rate 
probDeath <- 3 # bee death rate
probDep <- 30 # probability of depositing virus on flower

probAquireInfected <- 3
probScen <- 2
probFlow <- 3

colsBees <- c("white", "yellow", "orange", "pink", "red")
#----------------------------------------------------------------------------




# initialize random matrix with starting proportions of HB, BB individs:
beeVec <-sample(0:4, xDim*yDim, replace=T, prob = c(0.74, 0.05, 0.15, 0.1, 0.01))
beeMat <- matrix(data = beeVec, nrow = yDim, ncol = xDim)


# initialize random matrix with starting proportions of HB, BB colonies, and flowers:
colVec <-sample(0:6, xDim*yDim, replace=T, prob = c(0.85, 0.01, 0.01, 0.01, 0.01, 0.1, 0.01))
colMat <- matrix(data = colVec, nrow = yDim, ncol = xDim)


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
counter <- matrix(nrow=TimeSteps, ncol=5) 
counterCols <- matrix(nrow=TimeSteps, ncol=7)

# loop through time steps
for (t in 1:TimeSteps){
  
  # store table of counts in each row of counter matrix
  counter[t,] <- as.vector(table(beeMat))
  counterCols[t,] <- as.vector(table(colMat))
  
  # loop through each cell in the matrix
  for (i in 1:xDim){
    for (j in 1:yDim){

        # create random unit movement steps for i and j
        temp1 <- sample(c(1, -1), 1 ,replace=T, prob = c(0.5, 0.5))
        temp2 <- sample(c(1, -1), 1 ,replace=T, prob = c(0.5, 0.5))

        # Move the Bees:  
        if(beeMat[i,j]>=1){
          
          a <- ifelse(sum(i,temp1) %in% 1:xDim, sum(i,temp1), round(runif(1, 1, xDim)))
          b <- ifelse(sum(j,temp2) %in% 1:yDim, sum(j,temp2), round(runif(1, 1, yDim)))
          
          if(beeMat[a, b]==0){ 
            beeMat[a, b] <- beeMat[i,j]
            beeMat[i,j] <- 0
          }
        }
        
        # Kill the Bees:  
        if(beeMat[i,j]>=1){
          if(runif(1, 1,100)<=probDeath){
            beeMat[i,j] <- 0
          }
        }
        
        # Make the Bees:  
        if(colMat[i,j]==1 || colMat[i,j]==2 || colMat[i,j]==3 || colMat[i,j]==4){
          if(runif(1, 1,100)<=probBirth){
            beeMat[i,j] <- colMat[i,j]
          }
        }
        
        # Honey Bees Deposite Virus on flowers
        if(beeMat[i,j]==3 & colMat[i,j]==5){
          if(runif(1, 1,100)<=probDep){
            colMat[i,j] <- 6
          }
        }
        
        # Bumble Bees Pick up Virus
        if(beeMat[i,j]==2 & colMat[i,j]==6){
          if(runif(1, 1,100)<=probAquireInfected){
            beeMat[i,j] <- 4
          }
        }
        
        # Kill the Flowers:  
        if(colMat[i,j]>=5){
          if(runif(1, 1,100)<=probScen){
            colMat[i,j] <- 0
          }
        }
        
        # Make the Flowers:  
        if(colMat[i,j]==0){
          if(runif(1, 1,100)<=probFlow){
            colMat[i,j] <- 5
          }
        }
      
        
    } # end j loop
} # end i loop

  
# create names for each png:
name <- rename(t)
png(name)
 
# create each image to visulaize as a matrix 
image(1:nrow(beeMat), 1:ncol(beeMat), as.matrix(beeMat), col=colsBees, asp=1, xaxt='n', yaxt='n', ann=FALSE, bty='n')

dev.off()

} # end of time step for loop

#run ImageMagick: creates a gif of all images
my_command <- 'convert *.png -delay 3 -loop 0 animation.gif'
system(my_command)

dev.off()






ResultMat <- cbind(counter[,2:5], counterCols[,6:7])
x <- as.data.frame(ResultMat)
names(x) <- c("SHB", "SBB", "IHB", "IBB", "SF", "IF")

x$SHBprev <- x$SHB/(x$SHB+x$IHB)
x$SBBprev <- x$SBB/(x$SBB+x$IBB)
x$IHBprev <- x$IHB/(x$SHB+x$IHB)
x$IBBprev <- x$IBB/(x$SBB+x$IBB)
x$SFprev <- x$SF/(x$IF+x$SF)
x$IFprev <- x$IF/(x$IF+x$SF)

ResultMat <- as.matrix(x[,c(7:12)])

colors <- c("purple", "orange", "green", "red", "blue", "black")

matplot(y=ResultMat, type = "l", xlab = "Time",
        ylab = "Prevelence", lwd=3, 
        col=colors, ylim = c(0, 1), 
        lty=1, cex.lab = 1.3)

grid()

legend(197, 1.04, 
       legend=c("S HB", "S BB", "I HB"), 
       col=colors, lty=1, cex=.8)

legend(250, 1.04, 
       legend=c("I BB", "S FL", "I FL"), 
       col=colors[4:6], lty=1, cex=.8)





















