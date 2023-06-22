###### Tororo data analysis #####
### By JM Prada

## Load library and set workspace
## Load results from Mayuge as some posteriors will be used as priors here
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
load("Results.RData")

## Load data
dt <- read.csv("DataTororo.csv")

## Save variables of interest
## Note warning will appear when transforming CCA data to numeric - can be ignored
KK1 <- as.matrix(dt[,2:3])
KK2 <- as.matrix(dt[,4:5])
KK3 <- as.matrix(dt[,6:7])

CCA1 <- as.matrix(cbind(as.numeric(dt[,8]),dt[,9:10]))-1
CCA2 <- as.matrix(cbind(as.numeric(dt[,11]),dt[,12:13]))-1
CCA3 <- as.matrix(cbind(as.numeric(dt[,14]),dt[,15:16]))-1

## Prior for tau
library(fitdistrplus)

## Sample IDs from posterior
set.seed(1)
ndraws <- 1000
IDs <- sample.int(20000,ndraws)

## Prior for rt_intra_CCA
fit <- fitdist(c(Results$mcmc[[1]][,"rt_intra_CCA"],Results$mcmc[[2]][,"rt_intra_CCA"])[IDs],
               distr = "gamma", method = "mle")
plot(fit)
fit
## rt_intra_CCA fixed to 
## mean(rgamma(1000,fit$estimate[1],fit$estimate[2]))
## which is 1.41

## Prior rt_inter_CCA
fit <- fitdist(c(Results$mcmc[[1]][,"tau_inter"],
                 Results$mcmc[[2]][,"tau_inter"])[IDs],
               distr = "gamma", method = "mle")

plot(fit)
fit

## Prior for rt_inter
fit <- fitdist(c(Results$mcmc[[1]][,"rt_inter"],Results$mcmc[[2]][,"rt_inter"])[IDs],
               distr = "gamma", method = "mle")
plot(fit)
fit

## Prior rt_intra
fit <- fitdist(c(Results$mcmc[[1]][,"rt_intra"],
                 Results$mcmc[[2]][,"rt_intra"])[IDs],
               distr = "gamma", method = "mle")

plot(fit)
fit

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
      egg_excretion[n,2,d] ~ dnorm(infect_intensity[n],rt_inter)T(0,) #rt_inter, 0.2428324
      
      CCAintensity[n,1,d] ~ dnorm(0,tau_inter)T(-3,)
      CCAintensity[n,2,d] ~ dnorm(9 / (1 + exp(-multiParam[1]*(infect_intensity[n]-multiParam[2]))),tau_inter)T(0,)

    }
    
    infect_intensity[n] ~ dgamma(sh,rt)
    
    ## CCA component
    for (r in 1:CCAR){
      CCA1[n,r] ~ dnorm(CCAintensity[n,status[n]+1,1],1.41)T(0,9) #rt_intra_CCA
      CCA2[n,r] ~ dnorm(CCAintensity[n,status[n]+1,2],1.41)T(0,9) 
      CCA3[n,r] ~ dnorm(CCAintensity[n,status[n]+1,3],1.41)T(0,9) 
    }
  }
  
  ## Priors ##
  P ~ dbeta(1, 1)
  rt_inter ~ dgamma(7.289622,24.702766)
  rt_intra ~ dgamma(97.79942,102.28184)
  tau_inter ~ dgamma(204.6194,956.8694)
  rt_intra_CCA ~ dgamma(300.1179,212.0805)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  multiParam[1:2] ~ dmnorm.vcov(multim, multicov)
  
  #inits# .RNG.seed, .RNG.name, status
  #data# N, KK1, KK2, KK3, CCA1, CCA2, CCA3, KKR, CCAR, multim, multicov
  #monitor# P, sh, rt, rt_intra, tau_inter, rt_inter 
  #, rt_intra_CCA
}"

Results2 <- run.jags(m, burnin=10000, sample=10000, thin=10, n.chains=2, jags.refresh = 1, method = 'parallel',
                    plots = F, silent.jags = F)

plot(Results2)

## Save image after model run, commented to avoid accidentally overwriting results
#save.image("Results2.RData")

############################################################################
### Load RData file for analyses of model runs 
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
load("Results2.RData")

## Draw samples from the posteriors
set.seed(1)
ndraws <- 1000
IDs <- sample.int(20000,ndraws)

## See range Prevalence estimate
quantile(c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs],c(0.025,0.5,0.975))

## Calculate vars for KK
## variance negbin = mu x (mu + rt)/rt
## variance norm = 1/rt (as JAGS uses precision)
var_inter_KK <- 1/c(Results2$mcmc[[1]][,"rt_inter"],Results2$mcmc[[2]][,"rt_inter"])[IDs]

## For intra KK we just take the mean of a given day for the average individual
## (so mean of gamma at pop level)
iinten <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]/
  c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]

var_intra_KK <- iinten * (iinten + c(Results2$mcmc[[1]][,"rt_intra"],
                                     Results2$mcmc[[2]][,"rt_intra"])[IDs])/
  c(Results2$mcmc[[1]][,"rt_intra"],
    Results2$mcmc[[2]][,"rt_intra"])[IDs]

mean(var_intra_KK/var_inter_KK) #12 times the variance of intra is higher than inter
sd(var_intra_KK/var_inter_KK) #~3.8 sd in the mean

## calculate vars for CCA
var_inter_CCA <- 1/c(Results2$mcmc[[1]][,"tau_inter"],Results2$mcmc[[2]][,"tau_inter"])[IDs]
var_intra_CCA <- 1/1.41

mean(var_intra_CCA/var_inter_CCA) #0.11 times the variance of intra is higher than inter
sd(var_intra_CCA/var_inter_CCA) #0.0056 sd in the mean

## Calculate 95% quantiles
## KK * 24 to transform from eggs to Eggs per gram (EPGs)
quantile(var_inter_KK,c(0.025,0.5,0.975))*24
quantile(var_intra_KK,c(0.025,0.5,0.975))*24
quantile(var_inter_CCA,c(0.025,0.5,0.975))



####################################################################
## Do plots
library(scales)

pdf("SDKKandCCATororo.pdf",width=7*2,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfcol=c(1,2))

plot(NA,xlim = c(0,100),ylim = c(0,1),axes = F, xlab = "", ylab = "",pch=19)
polygon(density(var_intra_KK)$x,density(var_intra_KK)$y/max(density(var_intra_KK)$y),
        col = alpha("grey",alpha=0.5), border = NA)
#segments(1/0.2428324,0,1/0.2428324,1,col = alpha("red",alpha=0.5))
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
#polygon(density(var_intra_CCA)$x,density(var_intra_CCA)$y/max(density(var_intra_CCA)$y),
#        col = alpha("grey",alpha=0.5), border = NA)
segments(1/1.41,0,1/1.41,1,col = alpha("grey",alpha=0.5))
polygon(density(var_inter_CCA)$x,density(var_inter_CCA)$y/max(density(var_inter_CCA)$y),
        col = alpha("red",alpha=0.5), border = NA)

axis(1)
axis(2)

mtext("Scaled Density",side=2,cex=1,line=1.2)
mtext("Variance of POC-CCA",side=1,cex=1,line=1.2)

legend("topright",c("Inter-day","Intra-day"),col=alpha(c("red","grey"),alpha=0.5),
       bty='n',cex=0.75,lty=1)

dev.off()

##############################################
### Prevalence Figure

KKprev <- mean(sapply(1:nrow(KK1),function(x) 
{ifelse(sum(c(KK1[x,],
              KK2[x,],
              KK3[x,]))>0,1,0)}
)
, na.rm = T)

CCAprevH <- mean(sapply(1:nrow(KK1),function(x)
{ifelse(mean(c(CCA1[x,],
               CCA2[x,],
               CCA3[x,]), na.rm = T)>3,1,0)})
, na.rm = T)

CCAprevL <- mean(sapply(1:nrow(KK1),function(x)
{ifelse(mean(c(CCA1[x,],
               CCA2[x,],
               CCA3[x,]), na.rm = T)>2,1,0)})
, na.rm = T)

ModelP <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
ModelPrev <- quantile(ModelP,c(0.025,0.5,0.975))

pdf("PrevTororo.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim = c(0,2),ylim = c(0,1),axes = F, xlab = "", ylab = "",pch=19)

points(0,KKprev, pch=19,col="black")
points(1,CCAprevH, pch = 19, col = "red")
points(1,CCAprevL, pch = 19, col = "green")
points(2,ModelPrev[2], pch = 19)
arrows(2,ModelPrev[1],2,ModelPrev[3],angle=90,code = 3, length = 0.1)

axis(1, at= 0:2, labels = c("Kato-Katz", "POC-CCA", "Model"))
axis(2, at=seq(0,1,by=.2), labels = seq(0,1,by=.2)*100)

mtext("Prevalence (%)",side=2,cex=1,line=1.2)
#mtext("Variance of G-score",side=1,cex=1,line=1.2)

dev.off()