#============================================================
# On the Pitman–Yor process with spike and slab base measure
# by A. CANALE, A. LIJOI, B. NIPOTI AND I. PRUENSTER
# simulation experiment on density estimation
# (c) A. Canale, Nov 2016
#============================================================

n <- 50
atom.tau <- (1/0.2)^2
atom.mu <- 0 

set.seed(1)
ind1 <- sample(1:5, n, prob=c(0.1,0.1,0.2,0.2,0.4), replace=TRUE)
sample1 <- function(ind, mu, sigma) rnorm(1, mu[ind], sigma[ind])
y <- sapply(ind1, sample1, mu = c(-3.5,3.5,1,-1,0), sigma = c(1,1,0.8,0.8,0.2))

library(bayesm)
source('gibbs_density.R')
source("gen_fact.R")

# inner mixture with random zeta DP case
res1 <- mix.fit(y,  prior="DP", theta=1, sigma=0, a.tau=.5, b.tau=2, 
                  mu0=0, k=20, atom.mu=atom.mu, atom.tau=atom.tau, nrep=12000, nb=2000, randomzeta=TRUE)

# inner mixture  with random zeta PY case
res2 <- mix.fit(y, prior="PY", theta=1, sigma=.5, a.tau=.5, b.tau=2, 
                  mu0=0, k=20, atom.mu=atom.mu, atom.tau=atom.tau, nrep=12000, nb=2000, randomzeta=TRUE)

# inner mixture with fixed (wrong) zeta 
res3 <- mix.fit(y, prior="PY", theta=1, sigma=0.5, a.tau=.5, b.tau=2, 
                   mu0=0, k=20, atom.mu=atom.mu, atom.tau=atom.tau, nrep=12000, nb=2000, zeta0=0.8, randomzeta=FALSE)

# outer mixture with fixed (wrong) zeta
res4 <- mix.fit(y, prior="outer", theta=1, sigma=0.5, a.tau=.5, b.tau=2, 
                  mu0=0, k=20, atom.mu=atom.mu, atom.tau=atom.tau, nrep=12000, nb=2000, zeta0=0.8, randomzeta=FALSE)