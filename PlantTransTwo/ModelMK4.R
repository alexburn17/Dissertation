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



ModelFunc <- function(xDim=xDim1, 
                      yDim=yDim1,
                      TimeSteps=TimeSteps1,
                      probBirth=probBirth1,
                      probDeath=probDeath1,
                      probDep=probDep1,
                      probAquireInfected=probAquireInfected1,
                      probScen=probScen1,
                      probFlow=probFlow1, 
                      beeMat=beeMat1,
                      colMat=colMat1, 
                      counter=counter1, 
                      counterCols=counterCols1){
  
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
        if(colMat[i,j]==1 || colMat[i,j]==2 || colMat[i,j]==3){
          if(runif(1, 1,100)<=probBirth){
            beeMat[i,j] <- colMat[i,j]
          }
        }
        
        # Honey Bees Deposite Virus on flowers
        if(beeMat[i,j]==3 & colMat[i,j]==5){
          if(runif(1, 1,100)<=sample(probDep, 1, replace=T, prob = c(rep(1/length(probDep), 
                                                                         length(probDep))))){
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
    
    
    
    
  } # end of time step for loop
  
  mu <- mean(counter[-c(1:50),5]/(counter[-c(1:50),5]+counter[-c(1:50),3]))  
  
  
  #return(list(counter, counterCols))
  return(mu)
  
}




# Paramters:
#----------------------------------------------------------------------------
TimeSteps1 <- 300 # number of time steps

xDim1 <- 50 # x dimension of matrix
yDim1 <- 50 # y dimension of matrix

probBirth1 <- 6 # bee birth rate 
probDeath1 <- 5 # bee death rate
probDep1 <- c(rnorm(10, mean=22, sd=10)) # probability of depositing virus on flower

probAquireInfected1 <- 30
probScen1 <- 3
probFlow1 <- 4

#----------------------------------------------------------------------------


# initialize random matrix with starting proportions of HB, BB individs:
beeVec <-sample(0:4, xDim1*yDim1, replace=T, prob = c(0.74, 0.05, 0.15, 0.1, 0.01))
beeMat1 <- matrix(data = beeVec, nrow = yDim1, ncol = xDim1)


# initialize random matrix with starting proportions of HB, BB colonies, and flowers:
colVec <-sample(0:6, xDim1*yDim1, replace=T, prob = c(0.85, 0.01, 0.01, 0.01, 0.01, 0.1, 0.01))
colMat1 <- matrix(data = colVec, nrow = yDim1, ncol = xDim1)




# initialize a matrix to store number of each state in matrix at each time step
counter1 <- matrix(nrow=TimeSteps1, ncol=5) 
counterCols1 <- matrix(nrow=TimeSteps1, ncol=7)


vecDiv <- vector()

for(i in 1:30){

probDep1 <- c(rnorm(10, mean=22, sd=20))

vecDiv[i] <- ModelFunc(probAquireInfected = 30, probDep = probDep1)

}





vecMono <- vector()

probDep1 <- c(22, 22)

for(i in 1:30){
  
  vecMono[i] <- ModelFunc(probAquireInfected = 30, probDep = probDep1)
  
}













counter <- x[[1]]
counterCols <- x[[2]]

ResultMat <- cbind(counter[,2:5], counterCols[,6:7])
x <- as.data.frame(ResultMat)
names(x) <- c("SHB", "SBB", "IHB", "IBB", "SF", "IF")

#x$SHBprev <- x$SHB/(x$SHB+x$IHB)
#x$SBBprev <- x$SBB/(x$SBB+x$IBB)
x$IHBprev <- x$IHB/(x$SHB+x$IHB)
x$IBBprev <- x$IBB/(x$SBB+x$IBB)
#x$SFprev <- x$SF/(x$IF+x$SF)
x$IFprev <- x$IF/(x$IF+x$SF)

ResultMat <- as.matrix(x[,c(7,8,9)])

colors <- c("green", "red", "blue")

matplot(y=ResultMat, type = "l", xlab = "Time",
        ylab = "Prevalence", lwd=3, 
        col=colors, ylim = c(0, 1), 
        lty=1, cex.lab = 1.3)

grid()

legend(230, 1, 
       legend=c("Infected HB", "Infected BB", "Infected FL"), 
       col=colors, lty=1, cex=.8)






















