


ls()
rm(list=ls())
ls()

library(deSolve)

NosemaModel2 <- function(t, state, parameters){
  with(as.list(c(state, parameters)), {
    
    dSdt <- -(S * P * beta) - (S * muA) 
    dI1dt <- (S * P * beta) - (I1 * muA) - (I1 * gamma) 
    dI2dt <- (I1 * P * gamma) - (I2 * muB)
    dPdt <- (I1 * alpha1) + (I2 * alpha2) - (P * theta) 
    
    return(list(c(dSdt,dI1dt,dI2dt,dPdt)))
  })
}

# set parameters
state<-c(S=1, I1=0.02, I2=0.00, P=0.0)
parameters <- c(#r=0.05,
  beta=0.4,
  alpha1=0.10,
  alpha2=0.15,
  gamma=0.05,
  muA=0.01,
  muB=0.02,
  theta=0.05
)
times <- seq(0,150,by=1)
out <- ode(y=state,times=times, func=NosemaModel2, parms=parameters)
out<-as.data.frame(out)
out$time <- NULL
out$S <- NULL
NosemaLS <- out
NosemaLS <- (NosemaLS[,3])
out$P <- NULL
Surv <- 1 - (out[,1]+out[,2])
out <- cbind(out,Surv)
head(out,10)

plot(x=times,y=(2.5*(1-Surv)),
     type="l",
     #xlab="Time", 
     #ylab="Abundance", 
     #main="Abundance and Caste Composition of Bombus spp.",
     ylim=c(0,1),
     lwd=3,
     lty=1,
     font.lab=2,
     bty="l", 
     col="black",
     xaxt = "n", 
     yaxt = "n")



