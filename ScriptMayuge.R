###### Mayuge data analysis #####
### By JM Prada

## Load library and set workspace
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

## Load data
dt <- read.csv("DataMayuge.csv")

## Save variables of interest
## Note warning will appear when transforming CCA data to numeric - can be ignored
KK1 <- as.matrix(dt[,2:3])
KK2 <- as.matrix(dt[,6:7])
KK3 <- as.matrix(dt[,10:11])

CCA1 <- as.matrix(cbind(as.numeric(dt[,14]),dt[,15:16]))-1
CCA2 <- as.matrix(cbind(as.numeric(dt[,17]),dt[,18:19]))-1
CCA3 <- as.matrix(cbind(as.numeric(dt[,20]),dt[,21:22]))-1

## Load data for priors of logistic function for CCA
dtprior <- read.csv("kints.csv")

## Fit a multivariate normal distribution for the two parameters
## of the logistic function for CCA
library(MGMM)
fitParams <- FitGMM(data = as.matrix(dtprior[,c(4,5)]))

multim <- fitParams@Mean
multicov <- fitParams@Covariance

## Set seed JAGS model ##
.RNG.seed <- function(chain)
  return( switch(chain, "1"= 1, "2"= 2) )
.RNG.name <- function(chain)
  return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )

## Initialize Status ##
N <- nrow(KK1)
KKR <- 2 # number of repeated KK
CCAR <- 3 # number of repeated CCA
status <- rep(1,N)
rt_inter <- .2
rt_intra <- 1

## Model definition ##
m <- "model {
  # Prior random walk #
  
  for (n in 1:N){
  
    ## status
    status[n] ~ dbern(P)
    
    ## KK component
     for (r in 1:KKR){
      KK1[n,r] ~ dnegbin(rt_intra/(egg_excretion[n,status[n]+1,1]+rt_intra),rt_intra)
      KK2[n,r] ~ dnegbin(rt_intra/(egg_excretion[n,status[n]+1,2]+rt_intra),rt_intra)
      KK3[n,r] ~ dnegbin(rt_intra/(egg_excretion[n,status[n]+1,3]+rt_intra),rt_intra)
    }
    
    for(d in 1:3){
      egg_excretion[n,1,d] <- 0
      egg_excretion[n,2,d] ~ dnorm(infect_intensity[n],rt_inter)T(0,)
      
      CCAintensity[n,1,d] ~ dnorm(0,tau_inter)T(-3,)
      CCAintensity[n,2,d] ~ dnorm(9 / (1 + exp(-multiParam[1]*(infect_intensity[n]-multiParam[2]))),tau_inter)T(0,)

    }
    
    infect_intensity[n] ~ dgamma(sh,rt)
    
    ## CCA component
    for (r in 1:CCAR){
      CCA1[n,r] ~ dnorm(CCAintensity[n,status[n]+1,1],rt_intra_CCA)T(0,9) 
      CCA2[n,r] ~ dnorm(CCAintensity[n,status[n]+1,2],rt_intra_CCA)T(0,9) 
      CCA3[n,r] ~ dnorm(CCAintensity[n,status[n]+1,3],rt_intra_CCA)T(0,9) 
    }
  }
  
  ## Priors ##
  P ~ dbeta(1, 1)
  rt_inter ~ dnorm(.3,10)T(0,)
  rt_intra ~ dunif(0,10^3)
  tau_inter ~ dnorm(.3,10)T(0,)
  rt_intra_CCA ~ dgamma(0.001,0.001)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  multiParam[1:2] ~ dmnorm.vcov(multim, multicov)
  
  #inits# .RNG.seed, .RNG.name, status, rt_intra, rt_inter
  #data# N, KK1, KK2, KK3, CCA1, CCA2, CCA3, KKR, CCAR, multim, multicov
  #monitor# P, sh, rt, rt_intra, rt_intra_CCA, rt_inter, tau_inter
}"

Results <- run.jags(m, burnin=5000, sample=5000, thin=1, n.chains=2, jags.refresh = 1, method = 'parallel',
                    plots = F, silent.jags = F)

plot(Results)

## Save image after model run, commented to avoid accidentally overwriting results
#save.image("Results.RData")

#####################################################################
### Load RData file for analyses of model runs 
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
load("Results.RData")

## Draw samples from the posteriors
set.seed(1)
ndraws <- 1000
IDs <- sample.int(20000,ndraws)

## See range Prevalence estimate
quantile(c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs],c(0.025,0.5,0.975))

## Calculate vars for KK
## variance negbin = mu x (mu + rt)/rt
## variance norm = 1/rt (as JAGS uses precision)
var_inter_KK <- 1/c(Results$mcmc[[1]][,"rt_inter"],Results$mcmc[[2]][,"rt_inter"])[IDs]

## For intra KK we just take the mean of a given day for the average individual
## (so mean of gamma at pop level)
iinten <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]/
  c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]

var_intra_KK <- iinten * (iinten + c(Results$mcmc[[1]][,"rt_intra"],
                                     Results$mcmc[[2]][,"rt_intra"])[IDs])/
  c(Results$mcmc[[1]][,"rt_intra"],
    Results$mcmc[[2]][,"rt_intra"])[IDs]

mean(var_intra_KK/var_inter_KK) #15 times the variance of intra is higher than inter
sd(var_intra_KK/var_inter_KK) #~8.5 sd in the mean

## calculate vars for CCA
var_inter_CCA <- 1/c(Results$mcmc[[1]][,"tau_inter"],Results$mcmc[[2]][,"tau_inter"])[IDs]
var_intra_CCA <- 1/c(Results$mcmc[[1]][,"rt_intra_CCA"],Results$mcmc[[2]][,"rt_intra_CCA"])[IDs]

mean(var_intra_CCA/var_inter_CCA) #0.26 times the variance of intra is higher than inter
sd(var_intra_CCA/var_inter_CCA) #0.022 sd in the mean

## Calculate 95% quantiles
## KK * 24 to transform from eggs to Eggs per gram (EPGs)
quantile(var_inter_KK,c(0.025,0.5,0.975))*24
quantile(var_intra_KK,c(0.025,0.5,0.975))*24
quantile(var_inter_CCA,c(0.025,0.5,0.975))
quantile(var_intra_CCA,c(0.025,0.5,0.975))

####################################################################
## Do plot
library(scales)

pdf("SDKKandCCAMayuge.pdf",width=7*2,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfcol=c(1,2))

plot(NA,xlim = c(0,100),ylim = c(0,1),axes = F, xlab = "", ylab = "",pch=19)
polygon(density(var_intra_KK)$x,density(var_intra_KK)$y/max(density(var_intra_KK)$y),
        col = alpha("grey",alpha=0.5), border = NA)
polygon(density(var_inter_KK)$x,density(var_inter_KK)$y/max(density(var_inter_KK)$y),
        col = alpha("red",alpha=0.5), border = NA)

axis(1, at = seq(0,100,by=20), labels = seq(0,100,by=20)*24)
axis(2)

mtext("Scaled Density",side=2,cex=1,line=1.2)
mtext("Variance of KK (epg)",side=1,cex=1,line=1.2)

legend("topright",c("Inter-day","Intra-day"),col=alpha(c("red","grey"),alpha=0.5),
       bty='n',cex=0.75,lty=1)

### G-Score
plot(NA,xlim = c(0,8),ylim = c(0,1),axes = F, xlab = "", ylab = "",pch=19)
polygon(density(var_intra_CCA)$x,density(var_intra_CCA)$y/max(density(var_intra_CCA)$y),
        col = alpha("grey",alpha=0.5), border = NA)
polygon(density(var_inter_CCA)$x,density(var_inter_CCA)$y/max(density(var_inter_CCA)$y),
        col = alpha("red",alpha=0.5), border = NA)

axis(1)
axis(2)

mtext("Scaled Density",side=2,cex=1,line=1.2)
mtext("Variance of POC-CCA",side=1,cex=1,line=1.2)

legend("topright",c("Inter-day","Intra-day"),col=alpha(c("red","grey"),alpha=0.5),
       bty='n',cex=0.75,lty=1)

dev.off()
