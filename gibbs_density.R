#============================================================
# On the Pitman–Yor process with spike and slab base measure
# by A. CANALE, A. LIJOI, B. NIPOTI AND I. PRUENSTER
# Gibbs samplers for density estimation 
# (c) A. Canale, Nov 2016
#============================================================


mix.fit <- function(y, prior="DP", theta, sigma=0, a.tau, b.tau, kappa=var(y), mu0, 
                    k=20, atom.mu=0, atom.tau=1, nrep, nb, zeta0=0.5, randomzeta=TRUE)
{
  # y     vector with the data
  # prior             string of text for the prior type: "DP", "PY" for inner model "outer" for outer model
  # theta             parameter of the PY
  # sigma             prameter of the PY
  # a.tau, b.tau      parameters for the diffuse gamma density
  # kappa, mu0        parameters for the diffuse normal density
  # k                 upper bound on the number of clusters (for memory allocation purpose only)
  # atom.mu, atom.tau values for the spike
  # zeta0             mass for the spike
  # nrep, nb          MCMC iterations and burn-in
  # randomzeta        logical. if TRUE, uniform prior on zeta.
  
  N <- length(y)
  clusters <- matrix(NA, N, nrep)
  clusters[,1] <- 1
  mu <- tau <- matrix(NA, nrep, k)
  mu[1, ] <- c(mean(y), rep(NA, k-1))
  tau[1, ] <- c(1/var(y), rep(NA, k-1))
  C <- gen_fact_coeff(N, sigma)
  parclus <- rep(N, length=nrep)
  zeta <- rep(zeta0, length=nrep)
  a_zeta <- b_zeta <- 1
  zeta[1] <- 0.5
  for(t in 2:nrep)
  {
    #where we are?
    cat(paste(ifelse(round(t/10,0)==t/10,t,"")),"")
    
    #Update cluster allocation
    if(prior=="PY") clusters[,t] <- clusterallocation.py(y, mu[t-1,], tau[t-1,], clusters[,t-1], 
                                                         zeta[t-1], theta, sigma, C, a.tau, b.tau, mu0, kappa)
    if(prior=="DP") clusters[,t] <- clusterallocation.dp(y, mu[t-1,], tau[t-1,], clusters[,t-1], 
                                                         zeta[t-1], theta, a.tau, b.tau, mu0, kappa)
    if(prior=="outer") clusters[,t] <- clusterallocation.outer.py(y, mu[t-1,], tau[t-1,], clusters[,t-1], 
                                                                  zeta[t-1], theta, sigma, a.tau, b.tau, mu0, kappa)
    
    #update spike mass
    parclus[t] <- sum(clusters[,t]==1)
    if(randomzeta)
    {
      zeta[t] <- rbeta(1, a_zeta + parclus[t], b_zeta + N - parclus[t])
    }
    
    #Update mu and tau from conditional posterior
    pars <- reshuffle.norm.mix(y, clusters[,t], k, a.tau, b.tau, mu0, kappa, atom.mu, atom.tau)
    tau[t,] <- pars[(k+1):(2*k)]
    mu[t, ] <- pars[1:k]
  }
  list(mu=mu, tau=tau,zeta.m=mean(zeta[nb:nrep]), prop0.m = mean(parclus[nb:nrep])/n, 
       kn.m=mean(apply(clusters[,nb:nrep],2,lengthtable)))
}
#########################################################
dmixnorm<-function(x,pi,mu,tau){	
			sum(pi*dnorm(x,mu,sqrt(tau^(-1))))
			}
#########################################################
rprecision<-function(ind,ahat,bhat)
{
#random generate the precision tau from each cluster from a gamma distribution
rgamma(1,ahat[ind],bhat[ind])
}

rmean<-function(ind,muhat,kappahat,tau)
{
#random generate the mean mu for each cluster from a normal distribution
rnorm(1,muhat[ind], sqrt( kappahat[ind]*tau[ind]^(-1) ) )
}

reshuffle.norm.mix <- function(y, S, k, a.tau, b.tau, mu0, kappa, atom.mu, atom.tau)
{
  #Update mu and tau from conditional posterior
  nh  <- c(tapply(y,factor(S,levels=1:k),length))
  nh   <- replace(nh,is.na(nh),0)	 
  ybar <- c(tapply(y,factor(S,levels=1:k),mean)  )
  ybar <- replace(ybar,is.na(ybar),0)	 
  ydev <- (nh-1) * c(tapply(y,factor(S,levels=1:k),var)  )
  ydev <- replace(ydev,is.na(ydev),0)	 
  ahat.tau <- a.tau + nh/2
  bhat.tau <- b.tau + 1/2 * ( ydev + (nh/(1+kappa*nh))*(ybar-mu0)^2 )
  tau <- c(atom.tau, as.numeric(lapply(2:k,rprecision,ahat=ahat.tau,bhat=bhat.tau)))
  kappahat <- (kappa^(-1) + nh)^(-1)
  muhat <- kappahat * ( (kappa^(-1)) * mu0 +nh*ybar)
  mu <- c(atom.mu ,as.numeric(lapply(2:k,rmean,muhat=muhat, kappahat=kappahat, tau=tau ) ))
  tau[nh==0] <- NA
  mu[nh==0] <- NA
  c(mu, tau)
}

# custer allocation 
clusterallocation.py <- 
  function(y, mu, tau, clusters, zeta, theta, sigma, C, a.tau, b.tau, mu0, kappa)
  {  
    N <- length(y)
    K <- length(mu)
    for(i in 1:N)
    {
      hi <- clusters[i]
      nh <- table(clusters)
      singleton <- nh[hi]==1
      nh[hi] <- nh[hi]-1 
      k <- length(nh)-singleton# length(table(clusters[-i]))
      n0 <- nh[1]
      #newatom
      if(!singleton) 
      {
        tau[k+1] <- rgamma(1, a.tau, b.tau)
        mu[k+1] <- rnorm(1, mu0, sqrt(kappa/tau[k+1]))
      }
      #likelihoods
      likelihoods <- dnorm(y[i], mu, 1/sqrt(tau), log=FALSE)
      likelihoods <- likelihoods[!is.na(likelihoods)]
      
      #predictive weigths
      if(n0==0)
      {
        pnew <-  (1-zeta) * (theta + (k-1)*sigma)/(theta + N-1) 
        pzero <-  zeta * (theta + (k-1)*sigma)/(theta + N-1) 
      }
      else
      {
      if((theta/sigma + k - 1)>0)
        {
        ldenom <- (1:n0)*log(zeta) + C[n0,1:n0] + 
        lgamma(theta/sigma + k - 1 + 1:n0) - lgamma(theta/sigma + k - 1)
      lnum1 <-  (1:n0)*log(zeta) + C[n0,1:n0] + 
        lgamma(theta/sigma + k     + 1:n0) - lgamma(theta/sigma + k)
      lnum2 <-  (1:(n0+1))*log(zeta) + C[(n0+1),1:(n0+1)] + 
        lgamma(theta/sigma + k - 1 + 1:(n0+1)) - lgamma(theta/sigma + k - 1)
      M <- max(max(ldenom), max(lnum1),max(lnum2))
      pnew <- (1-zeta) * (theta + (k-1)*sigma)/(theta + N-1) * sum(exp(lnum1-M))/sum(exp(ldenom-M))
      pzero <- 1/(theta + N-1) *  sum(exp(lnum2-M))/sum(exp(ldenom-M))
        }
      else{
        ldenom <- (1:n0)*log(zeta) + C[n0,1:n0] + 
          lgamma(theta/sigma + k - 1 + 1:n0) #- lgamma(theta/sigma + k - 1)
        lnum1 <-  (1:n0)*log(zeta) + C[n0,1:n0] + 
          lgamma(theta/sigma + k     + 1:n0) - lgamma(theta/sigma + k)
        M <- max(max(ldenom), max(lnum1))
        pnew <- (1-zeta) * (theta + (k-1)*sigma)/(theta + N-1) * sum(exp(lnum1-M))/sum(exp(ldenom-M)/gamma(theta/sigma + k - 1))
        pzero <- 1-pnew#1/(theta + N-1) *  sum(exp(lnum2-M))/sum(exp(ldenom-M)/gamma(theta/sigma + k - 1))
      }
      }
      pocc <- (nh[-1]-sigma)/(theta+N-1)
      predictive.weights <- c(pzero,pocc)
      if(singleton)  predictive.weights[hi] <- pnew
      if(!singleton) predictive.weights <- c(predictive.weights, pnew)
      clusters[i] <- sample(1:(k+1), 1, prob=predictive.weights*likelihoods)
      
      #fix the list of atoms if some cluster are now empty
      nh.new <- table(factor(clusters, levels=1:(k+1)))
      #nh.new <- nh
      #if(singleton) nh.new[clusters[i]] <- nh.new[clusters[i]] + 1 
      #if((!singleton) & (clusters[i] == (k+1)) ) nh.new[k+1] <- 1
      #if((!singleton) & (clusters[i] < (k+1)) ) nh.new[clusters[i]] <- nh.new[clusters[i]] + 1
      is.empty <- sum(nh.new==0)
      new.empty <- which(nh.new==0)
      #   cat("(nh=",nh.new, "new.em=",new.empty,"),")
      if(is.empty)
      {
        populated <- nh.new>0
        mu <- mu[populated]
        tau <- tau[populated]
        clusters[clusters>which(nh.new==0)] = clusters[clusters>which(nh.new==0)] - 1
      }
      else
      {
        populated <- nh.new>0
        mu <- mu[populated]
        tau <- tau[populated]
      }
      #reallocation for subject i done, restart for subject i+1
    }
    clusters
  }

# custer allocation 
clusterallocation.dp <- 
  function(y, mu, tau, clusters, zeta, theta, a.tau, b.tau, mu0, kappa)
  {  
    N <- length(y)
    K <- length(mu)
    for(i in 1:N)
    {
      hi <- clusters[i]
      nh <- table(clusters)
      singleton <- nh[hi]==1
      nh[hi] <- nh[hi]-1 
      k <- length(nh)-singleton
      n0 <- nh[1]
      #newatom
      if(!singleton) 
      {
        tau[k+1] <- rgamma(1, a.tau, b.tau)
        mu[k+1] <- rnorm(1, mu0, sqrt(kappa/tau[k+1]))
      }
      #likelihoods
      likelihoods <- dnorm(y[i], mu, 1/sqrt(tau), log=FALSE)
      likelihoods <- likelihoods[!is.na(likelihoods)]
      
      #predictive weigths
      pnew <- (1-zeta) * theta/(theta + N - 1) 
      pzero <- (zeta * theta + nh[1])/(theta + N - 1) 
      pocc <- nh[-1]/(theta+N-1)
      predictive.weights <- c(pzero,pocc)
      if(singleton)  predictive.weights[hi] <- pnew
      if(!singleton) predictive.weights <- c(predictive.weights, pnew)
      clusters[i] <- sample(1:(k+1), 1, prob=predictive.weights*likelihoods)
      
      #fix the list of atoms if some cluster are now empty
      nh.new <- table(factor(clusters, levels=1:(k+1)))
      is.empty <- sum(nh.new==0)
      new.empty <- which(nh.new==0)
      if(is.empty)
      {
        populated <- nh.new>0
        mu <- mu[populated]
        tau <- tau[populated]
        clusters[clusters>which(nh.new==0)] = clusters[clusters>which(nh.new==0)] - 1
      }
      else
      {
        populated <- nh.new>0
        mu <- mu[populated]
        tau <- tau[populated]
      }
      #reallocation for subject i done, restart for subject i+1
    }
    clusters
  }

# custer allocation for outer spike-&-slab mixture
clusterallocation.outer.py <- 
  function(y, mu, tau, clusters, zeta, theta, sigma, a.tau, b.tau, mu0, kappa)
  {  
    N <- length(y)
    K <- length(mu)
    for(i in 1:N)
    {
      hi <- clusters[i]
      nh <- table(clusters)
      singleton <- nh[hi]==1
      not.spike <- hi!=1
      nh[hi] <- nh[hi]-1 
      k <- length(nh)-singleton
      n0 <- nh[1]
      #newatom
      if(!singleton) 
      {
        tau[k+1] <- rgamma(1, a.tau, b.tau)
        mu[k+1] <- rnorm(1, mu0, sqrt(kappa/tau[k+1]))
      }
      #likelihoods
      likelihoods <- dnorm(y[i], mu, 1/sqrt(tau), log=FALSE)
      likelihoods <- likelihoods[!is.na(likelihoods)]
      
      #predictive weigths
      pnew <- (1-zeta) * (theta + k*sigma)/(theta + N-nh[1]-not.spike)
      pzero <- zeta
      pocc <- (1-zeta) * (nh[-1]-sigma)/(theta + N-nh[1]-not.spike)
      predictive.weights <- c(pzero,pocc)
      if(singleton)  predictive.weights[hi] <- pnew
      if(!singleton) predictive.weights <- c(predictive.weights, pnew)
      clusters[i] <- sample(1:(k+1), 1, prob=predictive.weights*likelihoods)
      
      #fix the list of atoms if some cluster are now empty
      nh.new <- table(factor(clusters, levels=1:(k+1)))
      is.empty <- sum(nh.new==0)
      new.empty <- which(nh.new==0)
      if(is.empty)
      {
        populated <- nh.new>0
        mu <- mu[populated]
        tau <- tau[populated]
        clusters[clusters>which(nh.new==0)] = clusters[clusters>which(nh.new==0)] - 1
      }
      else
      {
        populated <- nh.new>0
        mu <- mu[populated]
        tau <- tau[populated]
      }
      #reallocation for subject i done, restart for subject i+1
    }
    clusters
  }

######
lengthtable <- function(x) length(table(x))
######