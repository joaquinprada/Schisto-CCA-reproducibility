## Generate a truncated distribution - from Keith Goldfeld in R-bloggers
rnormt <- function(n, range, mu, s = 1) {
  
  # range is a vector of two values
  
  F.a <- pnorm(min(range), mean = mu, sd = s)
  F.b <- pnorm(max(range), mean = mu, sd = s)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qnorm(u, mean = mu, sd = s)
}

## Function to calculate the individual level infection intensity 
modelData <- function(nInd, ID){
  nInf <- rbinom(1,nInd,prob=P[ID])
  
  Infectlvl <- rgamma(nInf,shape = sh[ID],rate = rt[ID])
  
  Intensityday <- sapply(1:3,function(x){c(rnormt((nInd-nInf), range = c(-3,Inf), mu = 0, s = interCCAsd[ID]),
                                           rnormt(nInf, range = c(0,Inf), mu = 9 / (1 + exp(-multim[1]*(Infectlvl-multim[2]))), 
                                                  s = interCCAsd[ID]))})
  
  mIntensity1 <- Intensityday[,1]
  mIntensity2 <- rowMeans(Intensityday[,1:2])
  mIntensity3 <- rowMeans(Intensityday)
  
  return(c(length(which(mIntensity1>1))/nInd,
           length(which(mIntensity2>1))/nInd,
           length(which(mIntensity3>1))/nInd,
           length(which(mIntensity1>1.5))/nInd,
           length(which(mIntensity2>1.5))/nInd,
           length(which(mIntensity3>1.5))/nInd,
           length(which(mIntensity1>2))/nInd,
           length(which(mIntensity2>2))/nInd,
           length(which(mIntensity3>2))/nInd,
           length(which(mIntensity1>3))/nInd,
           length(which(mIntensity2>3))/nInd,
           length(which(mIntensity3>3))/nInd))
}

## Function to simulate sens and spec
modelSensSpec <- function(nInd, ID){
  nInf <- rbinom(1,nInd,prob=P[ID])
  
  Infectlvl <- rgamma(nInf,shape = sh[ID],rate = rt[ID])
  
  Intensityday <- sapply(1:3,function(x){c(rnormt((nInd-nInf), range = c(-3,Inf), mu = 0, s = interCCAsd[ID]),
                                           rnormt(nInf, range = c(0,Inf), mu = 9 / (1 + exp(-multim[1]*(Infectlvl-multim[2]))), 
                                                  s = interCCAsd[ID]))})
  
  mIntensity1 <- Intensityday[,1]
  mIntensity2 <- rowMeans(Intensityday[,1:2])
  mIntensity3 <- rowMeans(Intensityday)
  
  return(c(length(which(mIntensity1[(nInd-nInf+1):nInd]>1,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<=1,1,0))/(nInd-nInf),
           length(which(mIntensity1[(nInd-nInf+1):nInd]>1.5,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<=1.5,1,0))/(nInd-nInf),
           length(which(mIntensity1[(nInd-nInf+1):nInd]>2,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<=2,1,0))/(nInd-nInf),
           length(which(mIntensity1[(nInd-nInf+1):nInd]>3,1,0))/nInf,
           length(which(mIntensity1[1:(nInd-nInf)]<=3,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>1,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<=1,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>1.5,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<=1.5,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>2,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<=2,1,0))/(nInd-nInf),
           length(which(mIntensity2[(nInd-nInf+1):nInd]>3,1,0))/nInf,
           length(which(mIntensity2[1:(nInd-nInf)]<=3,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>1,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<=1,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>1.5,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<=1.5,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>2,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<=2,1,0))/(nInd-nInf),
           length(which(mIntensity3[(nInd-nInf+1):nInd]>3,1,0))/nInf,
           length(which(mIntensity3[1:(nInd-nInf)]<=3,1,0))/(nInd-nInf)))
}
