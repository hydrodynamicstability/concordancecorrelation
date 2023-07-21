##########################################################################################
#
# Generate plots for Correlation manuscript
# 
# C. Tong -- 22 Jan 2011
#
# R version 2.12.1, packages listed below
#
# Slightly revised 30 June - 20 July 2023 using R version 3.6.1 and updated packages.
#
# Copyright by C. Tong 2023, all rights reserved.
# 
##########################################################################################

# An earlier version of this code was used to generate all the plots and simulation results shown in this JSM conference paper:
# C. Tong (2011), "Concordance correlation decomposed into the product of precision and accuracy", JSM 2011 Proceedings, Biometrics Section, pp. 649-662.

# For context, you may download a free copy of this paper here:  
# https://www.academia.edu/8840828/Concordance_correlation_coefficient_decomposed_into_the_product_of_precision_and_accuracy
#

options(help_type ="html")

# setwd("...") 

library(mvtnorm) # ver 0.4.2
library(ellipse) # version 1.1-1 

sessionInfo()


##########################################################################
###
### R function to generate simulated data from multivariate normal
### Downloaded from http://maven.smith.edu/~nhorton/R/
### Site of Prof. Nicholas J. Horton, Smith College
### August 17, 2007
### See Horton et al., Amer. Stat., 48:  343-357 (2004).
###
### 2023 update:  the above server does not seem to be available anymore,
### as Prof. Horton has moved to another institution.  I have reproduced
### his function here, with his gracious permission.
###
##########################################################################



rmultnorm <- function(n, mu, vmat, tol = 1e-07)  
    # a function to generate random multivariate Gaussians 
    {
       p <- ncol(vmat)  
       if (length(mu)!=p)
           stop("mu vector is the wrong length")  
       if (max(abs(vmat - t(vmat))) > tol) 
           stop("vmat not symmetric") 
       vs <- svd(vmat)  
       vsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d))) 
       ans <- matrix(rnorm(n * p), nrow = n) %*% vsqrt
       ans <- sweep(ans, 2, mu, "+")
       dimnames(ans) <- list(NULL, dimnames(vmat)[[2]])  
       return(ans)
   }
   

##########################################################################
###
### R function to calculate first principal component (ODR regression line)
### C. Tong 18 June 2010
### Updated 26 dec 2010 to return eval's and Coleman correlation
###
### Note: this is an early version of the odr.r function posted at my GitHub archive.
### It is included here to document how the graphs in the manuscript were generated,
### but users are advised to use the more rugged version found at the GitHub archive
### https://github.com/hydrodynamicstability/SLR.when.both.variables.random
###
##########################################################################

odr <- function(x,y)
{
    # Variances and covariances
    syy <- var(y,na.rm=TRUE)
    sxx <- var(x,na.rm=TRUE)
    sxy <- cov(x,y,use="pairwise.complete.obs")
    
    # Slope and intercept:
    m2 <- (syy - sxx + sqrt((syy - sxx)^2 + 4*sxy^2))/(2*sxy)
    int <- mean(y,na.rm=TRUE) - m2*mean(x,na.rm=TRUE)
    
    # Eigenvalues
    eval1 <- (syy + sxx + sqrt((syy - sxx)^2 + 4*sxy^2))/2
    eval2 <- (syy + sxx - sqrt((syy - sxx)^2 + 4*sxy^2))/2
    
    # Ratio of square roots of eigenvalues; eccentricity
    gof <- sqrt(eval2/eval1)
    eccentricity <- sqrt(1 - (eval2/eval1))
    
    # Coleman correlation on raw data
    rc <- 1 - sqrt(eval2/eval1)
    
    # return
    list(int=int,slope=m2,gof=gof,eccentricity=eccentricity,eval1=eval1,eval2=eval2,rc=rc)
}



##########################################################################################
#
# Simulate data set from multivariate normal
# 
##########################################################################################

# Properties of the population distribution
center <- c(5,5)
vmat <- matrix(c(1.5,1.8,1.8,1.5),2,2,byrow=TRUE)

# run simulation
set.seed(33)
sim.data <- rmultnorm(50,center,vmat)

# clean up results
colnames(sim.data) <- c("x","y")
sim.data <- as.data.frame(sim.data)


##########################################################################################
#
# Plot simulated data set (building Fig. 1a of my paper)
# 
##########################################################################################

# Settings for eps file:
cexpt <- 1.5
xpos <- 4.5

# for JPG it was xpos=4, cexpt=3

postscript(file="rotation.eps",width=12,height=4)
# jpeg(file="rotation.jpg",width=1440,height=480)

par(font.lab=cexpt,font.axis=2)
par(cex.lab=cexpt,cex.axis=cexpt)
bmax <- 11 # to define plotting region

par(mfrow=c(1,3))
plot(sim.data$x,sim.data$y,xlab="x",ylab="y",main="",xlim=c(0,bmax),ylim=c(0,bmax))
abline(0,1,lty=1)

text(xpos,10.5,"(a)\nPearson correlation = 0.78",cex=cexpt,font=2)

##########################################################################################
#
# First rotation (Fig. 1b of my paper)
# 
##########################################################################################

# setup rotation matrix
angle <- 15*pi/180
rot.mtx <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),2,2,byrow=TRUE)

# create rotated version of simulated data
sim.data2 <- rot.mtx %*% t(as.matrix(sim.data))
sim.data2 <- as.data.frame(t(sim.data2))
colnames(sim.data2) <- c("x","y")

# plot results
plot(sim.data2$x,sim.data2$y,xlab="x",ylab="y",main="", xlim=c(0,bmax),ylim=c(0,bmax))
abline(0,1,lty=1)

abline(0,tan((45*pi/180 - angle)),lty=2)

text(xpos,10.5,"(b)\nPearson correlation = 0.67",cex=cexpt,font=2)


##########################################################################################
#
# Second rotation (Fig. 1c of my paper)
# 
##########################################################################################

# setup rotation matrix
angle <- 30*pi/180
rot.mtx <- matrix(c(cos(angle),sin(angle),-sin(angle),cos(angle)),2,2,byrow=TRUE)

# created rotated version of simulated data
sim.data3 <- rot.mtx %*% t(as.matrix(sim.data))
sim.data3 <- as.data.frame(t(sim.data3))
colnames(sim.data3) <- c("x","y")

# plot results
plot(sim.data3$x,sim.data3$y,xlab="x",ylab="y",main="", xlim=c(0,bmax),ylim=c(0,bmax))
abline(0,1,lty=1)

abline(0,tan((45*pi/180 - angle)),lty=2)

sim.odr3 <- odr(sim.data3$x,sim.data3$y)
# abline(sim.odr3$int,sim.odr3$slope,lty=2)
# I will not add the ODR regression line to this plot.

text(xpos,10.5,"(c)\nPearson correlation = 0.33",cex=cexpt,font=2)

dev.off()




##########################################################################################
#
# Plot figure 2:  ellipse with theta; transformation for standardized Coleman correlation
# 
##########################################################################################


# Settings for eps file:
cexpt <- 1.5

# This is equation 11 in my paper
rcs.trans <- function(r)
{
    1 - sqrt((1 - abs(r))/(1+abs(r)))
}

# parameters for ellipse in Fig. 2a.
cvh <- (-121 - 49*sqrt(3) + 11*sqrt(2*(134+49*sqrt(3))))/108
cvhy <- (255+121*sqrt(3) - 11*sqrt(804+294*sqrt(3)))/108
cmatrix2 <- matrix(c(1,cvh,cvh,cvhy),nrow=2)
m <- (cvhy - 1 + sqrt((cvhy - 1)^2 + 4*cvh^2))/(2*cvh)
mm <- (cvhy - 1 - sqrt((cvhy - 1)^2 + 4*cvh^2))/(2*cvh)

elps.s1 <- ellipse(cmatrix2,centre=c(0,0))

# create fig. 2
postscript(file="ellipseRotation.eps",width=12,height=6,onefile=TRUE)

par(mfrow=c(1,2),cex.lab=cexpt,font=2)
# fig 2a
plot(elps.s1,type="l",xlim=c(-3,3),ylim=c(-3,3))
abline(h=0)
abline(v=0)
abline(0,m,lty=2)
abline(0,mm,lty=2)
text(1,0.2,substitute(theta) )
text(-2.5,2.5,"(a)",cex=cexpt,font=2)

# create fig 2b
curve(rcs.trans(x),from=0,to=1,xlab="r",ylab="Standardized Coleman correlation")
abline(0,1,lty=2)
text(0.05,0.95,"(b)",cex=cexpt,font=2)

dev.off()



##########################################################################
###
### Simulation of performance of Coleman correlation (Sec. 5 of my paper)
###
##########################################################################



##########################################################################################
#
# Simulate data set
# 
##########################################################################################



# parameters of the simulation

set.seed(19)
center <- c(5,5)
bmax <- 11 # for plotting region
theta <- 45
# theta <- 0
k <- 0.5
target.rc <- 1 - sqrt(k)
num.pts <- 1000

angle <- theta*pi/180
rot.mtx <- matrix(c(cos(angle),-sin(angle),sin(angle),cos(angle)),2,2,byrow=TRUE)

vmat <- matrix(c(1,0,0,k),2,2,byrow=TRUE)
vvmat <- rot.mtx %*% vmat %*% solve(rot.mtx)


# run trial simulation, clean up results
sim.data <- rmultnorm(num.pts,center,vvmat)
colnames(sim.data) <- c("x","y")
sim.data <- as.data.frame(sim.data)

# graph results
plot(sim.data$x,sim.data$y,xlab="x",ylab="y",main="",xlim=c(0,bmax),ylim=c(0,bmax))
sodr <- odr(sim.data$x,sim.data$y)



# Simulation for real
n.runs <- 10000 # you may want to reduce this if you attempt to run this code, otherwise it may take a long time to run!
kvec <- seq(0.1,0.9,by=0.1)
angle.vec <- seq(0,90,by=10)
nvec <- seq(10,100,by=10)

# to store the results:
datablock <- numeric(length(kvec)*length(angle.vec)*length(nvec)*n.runs)
dim(datablock) <- c(length(kvec),length(angle.vec),length(nvec),n.runs)

center <- c(5,5)

set.seed(17)
system.time(
for (i in 1:length(kvec))
{
    target.rc <- 1-sqrt(kvec[i])
    for (j in 1:length(angle.vec))
    {
        angle <- angle.vec[j]*pi/180
        rot.mtx <- matrix(c(cos(angle),-sin(angle),sin(angle),cos(angle)),2,2,byrow=TRUE)
        vmat <- matrix(c(1,0,0,kvec[i]),2,2,byrow=TRUE)
        vvmat <- rot.mtx %*% vmat %*% solve(rot.mtx)
        for (k in 1:length(nvec))
        {
            for (p in 1:n.runs)
            {
            sim.data <- rmultnorm(nvec[k],center,vvmat)
            colnames(sim.data) <- c("x","y")
            sim.data <- as.data.frame(sim.data)
            sodr <- odr(sim.data$x,sim.data$y)
            datablock[i,j,k,p] <- sodr$rc - target.rc
            }
        }
    }
}
)

# 50 runs:  45 sec.
# 1000 runs:  15 min.
# 10000 runs:  3 hours

# save(datablock,file="db50.rda")
# save(datablock,file="db1000.rda")
# save(datablock,file="db10000.rda") # not used
# load("db50.rda")
# load("db1000.rda")
# load("db10000.rda")


# Visualize results


# Pick 40 degrees angle:

# ymin <- min(datablock[,5,,])
# ymax <- max(datablock[,5,,])
ymin <- -1
ymax <- 1

median.block <- fq.block <- tq.block <- numeric(length(kvec)*length(nvec))
dim(median.block) <- dim(fq.block) <- dim(tq.block) <- c(length(kvec),length(nvec))

par(mfrow=c(3,3))
for (i in 1:length(kvec))
{
bp <- boxplot(datablock[i,5,1,],
    datablock[i,5,2,],
    datablock[i,5,3,],
    datablock[i,5,4,],
    datablock[i,5,5,],
    datablock[i,5,6,],
    datablock[i,5,7,],
    datablock[i,5,8,],
    datablock[i,5,9,],
    datablock[i,5,10,],xaxt="n",xlab="Sample size",ylab="Estimated - Target Rc",ylim=c(ymin,ymax),
    main=paste("Distribution of estimation error for\nColeman correlation, theta=40 degrees, eigenvalue ratio=",kvec[i]))
axis(1,at=c(1:10),nvec)
abline(h=0,lty=2)

median.block[i,] <- bp$stats[3,]
fq.block[i,] <- bp$stats[2,]
tq.block[i,] <- bp$stats[4,]
}


# Is it similar for other angles?

# Pick 0 degrees
a <- 1
#ymin <- min(datablock[,a,,])
#ymax <- max(datablock[,a,,])
ymin <- -1
ymax <- 1
x11()
par(mfrow=c(3,3))
for (i in 1:length(kvec))
{
boxplot(datablock[i,a,1,],
    datablock[i,a,2,],
    datablock[i,a,3,],
    datablock[i,a,4,],
    datablock[i,a,5,],
    datablock[i,a,6,],
    datablock[i,a,7,],
    datablock[i,a,8,],
    datablock[i,a,9,],
    datablock[i,a,10,],xaxt="n",xlab="Sample size",ylab="Estimated - Target Rc",ylim=c(ymin,ymax),
    main=paste("Distribution of estimation error for\nColeman correlation, theta=0 degrees, eigenvalue ratio=",kvec[i]))
axis(1,at=c(1:10),nvec)
abline(h=0,lty=2)

}


# Get quantiles for all data sets.
median.block <- fq.block <- tq.block <- numeric(length(angle.vec)*length(kvec)*length(nvec))
dim(median.block) <- dim(fq.block) <- dim(tq.block) <- c(length(angle.vec),length(kvec),length(nvec))

for (a in 1:length(angle.vec))
{
    for (i in 1:length(kvec))
    {
        bp <- boxplot(datablock[i,a,1,],
    datablock[i,a,2,],
    datablock[i,a,3,],
    datablock[i,a,4,],
    datablock[i,a,5,],
    datablock[i,a,6,],
    datablock[i,a,7,],
    datablock[i,a,8,],
    datablock[i,a,9,],
    datablock[i,a,10,],plot=FALSE)
    median.block[a,i,] <- bp$stats[3,]
    fq.block[a,i,] <- bp$stats[2,]
    tq.block[a,i,] <- bp$stats[4,]

    }
}


# Visualize quantile data.  Medians

par(mfrow=c(3,3))
for (i in 1:length(kvec))
{
plot(median.block[1,i,],xlab="Sample size",ylab="Median error",ylim=c(-1,1),xaxt="n",type="b",main=kvec[i])
axis(1,at=c(1:10),nvec)
colvec <- rainbow(length(angle.vec))
for (a in 1:length(angle.vec))
{
    points(median.block[a,i,],type="b",col=colvec[a])
}
}

# First quartiles
par(mfrow=c(3,3))
for (i in 1:length(kvec))
{
plot(fq.block[1,i,],xlab="Sample size",ylab="First quartile error",ylim=c(-1,1),xaxt="n",type="b",main=kvec[i])
axis(1,at=c(1:10),nvec)
colvec <- rainbow(length(angle.vec))
for (a in 1:length(angle.vec))
{
    points(fq.block[a,i,],type="b",col=colvec[a])
}
}

# Third quartiles
par(mfrow=c(3,3))
for (i in 1:length(kvec))
{
plot(tq.block[1,i,],xlab="Sample size",ylab="Third quartile error",ylim=c(-1,1),xaxt="n",type="b",main=kvec[i])
axis(1,at=c(1:10),nvec)
colvec <- rainbow(length(angle.vec))
for (a in 1:length(angle.vec))
{
    points(tq.block[a,i,],type="b",col=colvec[a])
}
}


# Plots for publication

# Settings for eps file:
cexpt <- 1.5
xpos <- 4.5

ymin <- -1.2
ymax <- 1.2
mainvec <- c("(a) k = 0.1","(b) k = 0.5","(c) k = 0.9")



postscript(file="rcsimulation.eps",width=6,height=90,onefile=TRUE)

par(mfrow=c(3,1),cex.lab=cexpt,cex.axis=cexpt,font=2,font.lab=cexpt,font.axis=2)
k <- 1
for (i in c(1,5,9))
{
bp <- boxplot(datablock[i,5,1,],
    datablock[i,5,2,],
    datablock[i,5,3,],
    datablock[i,5,4,],
    datablock[i,5,5,],
    datablock[i,5,6,],
    datablock[i,5,7,],
    datablock[i,5,8,],
    datablock[i,5,9,],
    datablock[i,5,10,],xaxt="n",xlab="Sample size",ylab="Estimated - Population",ylim=c(ymin,ymax))
axis(1,at=c(1:10),nvec,font.axis=2,cex.axis=cexpt)
abline(h=0,lty=2)

text(1,1,mainvec[k],cex=cexpt,font=2)
k <- k+1
}

dev.off()

# Copyright (C) 2023 by Christopher Tong

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
