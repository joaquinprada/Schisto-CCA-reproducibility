###### Simulations to evaluate the error #####
## By JM Prada

## Load library and set workspace
library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

## Load data
load("Results2.RData")
source("SimulationErrorFunctions.R")

## Set draws
set.seed(100)
ndraws <- 1000
IDs <- sample.int(20000,ndraws)

#### Mayuge
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tau_inter"],Results$mcmc[[2]][,"tau_inter"])[IDs])
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]

## Simulate prevalence estimates for Mayuge
PrevEstimatesMayuge <- sapply(1:length(IDs),function(x){modelData(1000,x)})

## Calculate quantiles for squared error
sderH <- apply((P-PrevEstimatesMayuge)^2,1,quantile,c(.025,.5,.975))

#### Tororo
interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tau_inter"],Results2$mcmc[[2]][,"tau_inter"])[IDs])
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]

## Simulate prevalence estimates for Tororo
PrevEstimatesTororo <- sapply(1:length(IDs),function(x){modelData(1000,x)})

## Calculate squared error
sderL <- apply((P-PrevEstimatesTororo)^2,1,quantile,c(.025,.5,.975))


######################################################
### Plot
pdf("ErrorPlot.pdf",width=7*2,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfcol=c(1,2))

plot(NA,xlim = c(1,12),ylim = c(0,0.1),axes = F, xlab = "", ylab = "")
points(sderH[2,], pch=19)
arrows(x0=1:12,y0=sderH[1,],y1=sderH[3,],angle = 90, code = 3, length = .1)
abline(v=3.5, lty = "dashed")
abline(v=6.5, lty = "dashed")
abline(v=9.5, lty = "dashed")

axis(1, at = 1:12, labels = rep(1:3,4))
axis(2)

mtext("Squared error",side=2,cex=1,line=1.2)
mtext("Number of days of sampling",side=1,cex=1,line=1.2)

text(2,0.095,labels = "Threshold = G2")
text(5,0.095,labels = "Threshold = G2.5")
text(8,0.095,labels = "Threshold = G3")
text(11,0.095,labels = "Threshold = G4")

#### Tororo
plot(NA,xlim = c(1,12),ylim = c(0,0.1),axes = F, xlab = "", ylab = "")
points(sderL[2,], pch=19)
arrows(x0=1:12,y0=sderL[1,],y1=sderL[3,],angle = 90, code = 3, length = .1)
abline(v=3.5, lty = "dashed")
abline(v=6.5, lty = "dashed")
abline(v=9.5, lty = "dashed")

axis(1, at = 1:12, labels = rep(1:3,4))
axis(2)

mtext("Squared error",side=2,cex=1,line=1.2)
mtext("Number of days of sampling",side=1,cex=1,line=1.2)

text(2,0.095,labels = "Threshold = G2")
text(5,0.095,labels = "Threshold = G2.5")
text(8,0.095,labels = "Threshold = G3")
text(11,0.095,labels = "Threshold = G4")

dev.off()


#############################################
#### Mayuge Sens/Spec
interCCAsd <- sqrt(1/c(Results$mcmc[[1]][,"tau_inter"],Results$mcmc[[2]][,"tau_inter"])[IDs])
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]

## Simulate Sensitivity and Specificity for Mayuge 
SensSpecEstimatesMayuge <- sapply(1:length(IDs),function(x){modelSensSpec(1000,x)})

## Calculate quantiles
QsenspcMayuge <- apply(SensSpecEstimatesMayuge,1,quantile,c(.025,.5,.975))

#### Tororo
interCCAsd <- sqrt(1/c(Results2$mcmc[[1]][,"tau_inter"],Results2$mcmc[[2]][,"tau_inter"])[IDs])
P <- c(Results2$mcmc[[1]][,"P"],Results2$mcmc[[2]][,"P"])[IDs]
sh <- c(Results2$mcmc[[1]][,"sh"],Results2$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results2$mcmc[[1]][,"rt"],Results2$mcmc[[2]][,"rt"])[IDs]

## Simulate Sensitivity and Specificity for Tororo
SensSpecEstimatesTororo <- sapply(1:length(IDs),function(x){modelSensSpec(1000,x)})

## Calculate quantiles
QsenspcTororo <- apply(SensSpecEstimatesTororo,1,quantile,c(.025,.5,.975))

### Plot ROC curve ###
pdf("Figure 2 - ROCcurve.pdf",width=7*2,height=4.64)
par(font=2, cex.axis=1, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))
par(mfcol=c(1,2))

## Mayuge
plot(NA,xlim = c(0,.6),ylim = c(.4,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcMayuge[2,c(2,4,6,8)],QsenspcMayuge[2,c(1,3,5,7)], pch=19)
lines(1-QsenspcMayuge[2,c(2,4,6,8)],QsenspcMayuge[2,c(1,3,5,7)])
arrows(x0=1-QsenspcMayuge[2,c(2,4,6,8)],y0=QsenspcMayuge[1,c(1,3,5,7)],
        y1=QsenspcMayuge[3,c(1,3,5,7)], angle=90, code = 3,length = .05)
arrows(x0=1-QsenspcMayuge[1,c(2,4,6,8)],y0=QsenspcMayuge[2,c(1,3,5,7)],
       x1=1-QsenspcMayuge[3,c(2,4,6,8)], angle=90, code = 3,length = .05)

## 2 Days
points(1-QsenspcMayuge[2,c(10,12,14,16)],QsenspcMayuge[2,c(9,11,13,15)], pch=19, col="red")
lines(1-QsenspcMayuge[2,c(10,12,14,16)],QsenspcMayuge[2,c(9,11,13,15)], col="red")
arrows(x0=1-QsenspcMayuge[2,c(10,12,14,16)],y0=QsenspcMayuge[1,c(9,11,13,15)],
       y1=QsenspcMayuge[3,c(9,11,13,15)], angle=90, code = 3,length = .05, col="red")
arrows(x0=1-QsenspcMayuge[1,c(10,12,14,16)],y0=QsenspcMayuge[2,c(9,11,13,15)],
       x1=1-QsenspcMayuge[3,c(10,12,14,16)], angle=90, code = 3,length = .05, col="red")

## 3 Days
points(1-QsenspcMayuge[2,c(18,20,22,24)],QsenspcMayuge[2,c(17,19,21,23)], pch=19, col="darkgreen")
lines(1-QsenspcMayuge[2,c(18,20,22,24)],QsenspcMayuge[2,c(17,19,21,23)], col="darkgreen")
arrows(x0=1-QsenspcMayuge[2,c(18,20,22,24)],y0=QsenspcMayuge[1,c(17,19,21,23)],
       y1=QsenspcMayuge[3,c(17,19,21,23)], angle=90, code = 3,length = .05, col="darkgreen")
arrows(x0=1-QsenspcMayuge[1,c(18,20,22,24)],y0=QsenspcMayuge[2,c(17,19,21,23)],
       x1=1-QsenspcMayuge[3,c(18,20,22,24)], angle=90, code = 3,length = .05, col="darkgreen")

axis(1)
axis(2)

mtext("Sensitivity",side=2,cex=1.25,line=1.2)
mtext("1 - Specificity",side=1,cex=1.25,line=1.2)

legend("bottomright",c("Sampling 1 Day","Sampling 2 Days","Sampling 3 Days"),
       col=c("black","red","darkgreen"),
       bty='n',cex=1.25,lty=1)

## Tororo
plot(NA,xlim = c(0,.6),ylim = c(.4,1),axes = F, xlab = "", ylab = "")
## 1 Day
points(1-QsenspcTororo[2,c(2,4,6,8)],QsenspcTororo[2,c(1,3,5,7)], pch=19)
lines(1-QsenspcTororo[2,c(2,4,6,8)],QsenspcTororo[2,c(1,3,5,7)])
arrows(x0=1-QsenspcTororo[2,c(2,4,6,8)],y0=QsenspcTororo[1,c(1,3,5,7)],
       y1=QsenspcTororo[3,c(1,3,5,7)], angle=90, code = 3,length = .05)
arrows(x0=1-QsenspcTororo[1,c(2,4,6,8)],y0=QsenspcTororo[2,c(1,3,5,7)],
       x1=1-QsenspcTororo[3,c(2,4,6,8)], angle=90, code = 3,length = .05)

## 2 Days
points(1-QsenspcTororo[2,c(10,12,14,16)],QsenspcTororo[2,c(9,11,13,15)], pch=19, col="red")
lines(1-QsenspcTororo[2,c(10,12,14,16)],QsenspcTororo[2,c(9,11,13,15)], col="red")
arrows(x0=1-QsenspcTororo[2,c(10,12,14,16)],y0=QsenspcTororo[1,c(9,11,13,15)],
       y1=QsenspcTororo[3,c(9,11,13,15)], angle=90, code = 3,length = .05, col="red")
arrows(x0=1-QsenspcTororo[1,c(10,12,14,16)],y0=QsenspcTororo[2,c(9,11,13,15)],
       x1=1-QsenspcTororo[3,c(10,12,14,16)], angle=90, code = 3,length = .05, col="red")

## 3 Days
points(1-QsenspcTororo[2,c(18,20,22,24)],QsenspcTororo[2,c(17,19,21,23)], pch=19, col="darkgreen")
lines(1-QsenspcTororo[2,c(18,20,22,24)],QsenspcTororo[2,c(17,19,21,23)], col="darkgreen")
arrows(x0=1-QsenspcTororo[2,c(18,20,22,24)],y0=QsenspcTororo[1,c(17,19,21,23)],
       y1=QsenspcTororo[3,c(17,19,21,23)], angle=90, code = 3,length = .05, col="darkgreen")
arrows(x0=1-QsenspcTororo[1,c(18,20,22,24)],y0=QsenspcTororo[2,c(17,19,21,23)],
       x1=1-QsenspcTororo[3,c(18,20,22,24)], angle=90, code = 3,length = .05, col="darkgreen")

axis(1)
axis(2)

mtext("Sensitivity",side=2,cex=1.25,line=1.2)
mtext("1 - Specificity",side=1,cex=1.25,line=1.2)

legend("bottomright",c("Sampling 1 Day","Sampling 2 Days","Sampling 3 Days"),
       col=c("black","red","darkgreen"),
       bty='n',cex=1.25,lty=1)

dev.off()


