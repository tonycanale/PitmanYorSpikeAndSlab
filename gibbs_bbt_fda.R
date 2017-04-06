#============================================================
# On the Pitman–Yor process with spike and slab base measure
# by A. CANALE, A. LIJOI, B. NIPOTI AND I. PRUENSTER
# bbt curve application
# (c) A. Canale, May 2016
#============================================================
rm(list=ls())

# load auxiliary functions
source("aux_fda.R")
source("splines_fda.R")
source("clustering_fda.R")
source("gen_fact.R")

# load libraries
library(mnormt)
library(splines)

# load data (Verona center)
# The original data bank may be accessed upon written request to the Department of Statistical Sciences
# of the University of Padua. See www.stat.unipd.it/en/fare-ricerca/basi-di-dati
#load("datiCICLOVerona.RData")
# a fake 2-row dataset looks like the following
temp <- matrix(NA, 2, 49)
temp[1,] <- c(rep(NA,6), 36, 35.8, 36.4, 36.45, 36.4,36.45, 36.55, 36.65, 36.75, 36.85, 36.8, 36.75, rep(NA, 31))
temp[2,] <- c(rep(NA,7), 35.2, 35.2, 35.3, 36, 36.1,36.2, 36.3, 36.45, 36.5, 36.65, 36.6, 36.6, 36.55, 36.6, rep(NA, 28))

# other variables 
T <- 49
t <- 1:T
tcont <- c(t,t[-T]+0.5)
# for the splines expansion put fixed knots over a reasonable domain
#knots <- seq(-24,24, length=7)
knots <- c(-18,-12,-6,0,6,12,18) #come bruno!

# mcmc interations
nrep <- 8000
nb <- 2000

#prior paramters
#gamma parameters
a_sigma2 <- 1/2
b_sigma2 <- 1/2
# PY parameters
sigma <- 0.25 
theta <- 1
#spike-n-slab paramters
a <- 0.10
c <- 0.90
#covariance matrix for (tau1, tau2)
OmegaInv0 <-  diag(1,2)
# next par are for Inv-Wishar (currently not uesed)
rho1 <- 1
rho2 <- 1
R <- diag(c(rho1,rho2))
#spline basis coefficients prior
d <- length(knots)+3
beta0 = rep(0, d)
SigmaInv0 = diag(1, d)
#delta prior parameter (beta)
a_delta <- 1
b_delta <- 1
  
#mcmc quantities  
tau1 <- tau2 <- omega <- lambda <- matrix(NA, N, nrep)
alpha_i <- array(NA, dim=c(n,2,nrep))
sigma2 <- rep(NA, nrep)
OmegaInv <- array(NA, dim=c(2,2,nrep))
clusters <- matrix(NA, N, nrep)
atoms <- list()
storefij <- array(NA, dim=c(N,49,nrep))
parclus <- rep(N, nrep)
delta <- rep(0.95, nrep)
  
#mcmc initializations
tau1[,1] <- apply(temp,1,min, na.rm=TRUE)
tau2[,1] <- apply(temp,1,max, na.rm=TRUE) - tau1[,1]
alpha <- c(mean(tau1[,1]), mean(tau2[,1]))
omega[,1] <- 1
ovulazione <- function(x) diff(range(which(!is.na(x))))
lambda[,1] <- apply(temp, 1, ovulazione)
alpha_i[,1,1] <- tapply(tau1[,1],FUN = mean,INDEX = ID)
alpha_i[,2,1] <- tapply(tau2[,1],FUN = mean,INDEX = ID)
sigma2[1] <- 2
OmegaInv <- 2*diag(2)
clusters[,1] <- 1
fbase <- function(x) plogis(x)
atoms[[1]] <- fbase
applyfbase <- function(cyclepar, domain=t) cyclepar[3] + cyclepar[4]*fbase((domain-cyclepar[1])/cyclepar[2])
fij <- t(apply(cbind(lambda[,1],omega[,1],0,1), 1, applyfbase))
tildefij <- t(apply(cbind(lambda[,1],omega[,1],tau1[,1],tau2[,1]), 1, applyfbase))
storefij[,,1] <- tildefij

# generalized factorial coefficients (just once for signa fixed)
C <- gen_fact_coeff(N, sigma)

# start the MCMC sampler
for(ite in 2:nrep){
#ite = 2

# cycle specific follicular and lutheal phase temperatures
taus <- sapply(1:N, sampletaus, fij, temp, OmegaInv, sigma2[ite-1], alpha_i[,,ite-1], ID)
tau1[,ite] <- taus[1,]
tau2[,ite] <- taus[2,]

# cycle specific first day of high temperature (ovulation)
lambda[,ite] <- sapply(1:N, samplelambdas, temp, t, tau1[,ite], tau2[,ite], 
                       omega[,ite-1], sigma2[ite-1], cluster=clusters[,ite-1], atoms)

# cycle specifc increase in temperature
omega[,ite] <- sapply(1:N, sampleomegas, temp, t, 1/2, 1, tau1=tau1[,ite], tau2=tau2[,ite], 
                      omegaprec=omega[,ite-1], lambda=lambda[,ite], sigma2=sigma2[ite-1], 
                      cluster=clusters[,ite-1], atoms=atoms)

# subject specific means for tau1 and tau2
alpha_i[,,ite] <- t(sapply(1:n, samplealphas, ni, tau1[,ite], tau2[,ite], 
                           OmegaInv, R, alpha, ID))

# cluster allocation 
clusters[,ite] <- clusterallocation(temp, atoms, clusters[,ite-1], a, c, theta, sigma, C, 
           tau1[,ite], tau2[,ite], lambda[,ite], omega[,ite], sigma2[ite-1], 
           beta0, SigmaInv0)

#update a and c 
parclus[ite] <- sum(clusters[,ite]==1)
zeta[ite] <- rbeta(1, a_delta + parclus[ite], b_delta + N - parclus[ite])
c <- delta[ite]
a <- 1-delta[ite]

# reshuffling atoms 
k <- length(table(clusters[,ite]))
if(k>1)
{
reshuffledatoms <- sapply(2:k, reshuffle, clusters=clusters[,ite], yij=temp, 
                          tau1=tau1[,ite], tau2=tau2[,ite], lambda=lambda[,ite], omega=omega[,ite], 
                          sigma2=sigma2[ite-1], SigmaInv0=SigmaInv0, beta0=beta0, 
                          domain=t, knots=knots)
atoms <- append(atoms[1], reshuffledatoms)
}
if(k==1)  atoms <- atoms[1]
#new cycle specific trajectories are
fij <- t(sapply(1:N, applyfij, tau1=rep(0,N), tau2=rep(1,N), lambda=lambda[,ite], 
                   omega=omega[,ite], clusters=clusters[,ite], atoms=atoms, domain=t))
tildefij <- t(sapply(1:N, applyfij, tau1=tau1[,ite], tau2=tau2[,ite], lambda=lambda[,ite], 
                   omega=omega[,ite], clusters=clusters[,ite], atoms=atoms, domain=t))
if(ite>nb) storefij[,,ite] <- tildefij

# measurement error common variance
sigma2[ite] <- 1/rgamma(1, a_sigma2 + 0.5*sum(!is.na(temp)), b_sigma2 + 0.5*sum((temp[!is.na(temp)] - tildefij[!is.na(temp)])^2))

# just print "i" if i-th iteration is OK!
cat(ite)
}


