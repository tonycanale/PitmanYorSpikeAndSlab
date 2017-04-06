# sample a new atom from the prior
samplenewatom <- function(theta0, SigmaInv0, knots)
{
  V <- solve(SigmaInv0)
  newbeta <- rmnorm(1, theta0, V)
  newbeta
  f <- function(x) bs(x, knots=knots) %*% newbeta
  f
}

# compute the value of the cycle specific likelihood
cluters_logliks <- function(h, yij, atoms, tau1, tau2, lambda, omega, sigma2, t)
{
  keep <- !is.na(yij)
  f <- atoms[[h]]
  mean <- tau1 + tau2*f( (t[keep] - lambda)/omega )
  sum(dnorm(yij[keep], mean, sqrt(sigma2), log=TRUE))
}

compute.likelihoods <- function(yij, atoms, tau1, tau2, lambda, omega, sigma2, t)
{
  kp1 <- dim(summary(atoms))[1]
  logliks <- sapply(1:kp1, cluters_logliks, yij=yij, atoms=atoms, 
                    tau1=tau1, tau2=tau2, lambda=lambda, omega=omega, sigma2=sigma2, t=t)
  exp(logliks)
}

# custer allocation 
clusterallocation <- 
  function(yij, atoms, clusters, a, c, theta, sigma, C, 
           tau1, tau2, lambda, omega, sigma2, theta0, SigmaInv0)
{  
  N <- nrow(yij)  
  #  cat("subject: ")
  for(i in 1:N)
  {
    #preprocess
    hi <- clusters[i]
    nh <- table(clusters)
    singleton <- nh[hi]==1
    nh[hi] <- nh[hi]-1 
    k <- length(nh)-singleton# length(table(clusters[-i]))
    n0 <- nh[1]
#    cat("\xn\nk=", k, "with length(table(clusters[-i]))=", length(table(clusters[-i])))
#    cat(", cluster = ", hi, "of (new) size ", nh[hi], "singleton=",singleton)
    #newatom
    if(!singleton) atoms[[k+1]] <- samplenewatom(theta0, SigmaInv0, knots)
    
    #likelihoods
    likelihoods <- compute.likelihoods(yij[i,], atoms, tau1[i], tau2[i], 
                                       lambda[i], omega[i], sigma2, t)
    
    #predictive weigths
    ldenom <- (1:n0)*log(c/(a+c)) + C[n0,1:n0] + 
      lgamma(theta/sigma + k - 1 + 1:n0) - lgamma(theta/sigma + k - 1)
    lnum1 <-  (1:n0)*log(c/(a+c)) + C[n0,1:n0] + 
      lgamma(theta/sigma + k     + 1:n0) - lgamma(theta/sigma + k)
    lnum2 <-  (1:(n0+1))*log(c/(a+c)) + C[(n0+1),1:(n0+1)] + 
      lgamma(theta/sigma + k - 1 + 1:(n0+1)) - lgamma(theta/sigma + k - 1)
    M <- max(max(ldenom), max(lnum1),max(lnum2))
    pnew <- a/(a+c) * (theta + (k-1)*sigma)/(theta + N-1) * sum(exp(lnum1-M))/sum(exp(ldenom-M))
    pzero <- 1/(theta + N-1) *  sum(exp(lnum2-M))/sum(exp(ldenom-M))
    pocc <- (nh[-1]-sigma)/(theta+N-1)
    predictive.weights <- c(pzero,pocc)
    if(singleton)  predictive.weights[hi] <- pnew
    if(!singleton) predictive.weights <- c(predictive.weights, pnew)
    #cat(predictive.weights)
    #cat("\n")
    #cat(likelihoods)
    #then update the cluser membership
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
      atoms <- atoms[populated]
      clusters[clusters>which(nh.new==0)] = clusters[clusters>which(nh.new==0)] - 1
    }
    else
    {
      populated <- nh.new>0
      atoms <- atoms[populated]
    }
    #reallocation for subject i done, restart for subject i+1
}
 clusters
}


clusterallocationDP <- 
  function(yij, atoms, clusters, a, c, theta, sigma, C, 
           tau1, tau2, lambda, omega, sigma2, theta0, SigmaInv0)
  {  
    N <- nrow(yij)  
    #  cat("subject: ")
    for(i in 1:N)
    {
      #preprocess
      hi <- clusters[i]
      nh <- table(clusters)
      singleton <- nh[hi]==1
      nh[hi] <- nh[hi]-1 
      k <- length(nh)-singleton# length(table(clusters[-i]))
      n0 <- nh[1]
      #    cat("\xn\nk=", k, "with length(table(clusters[-i]))=", length(table(clusters[-i])))
      #    cat(", cluster = ", hi, "of (new) size ", nh[hi], "singleton=",singleton)
      #newatom
      if(!singleton) atoms[[k+1]] <- samplenewatom(theta0, SigmaInv0, knots)
      
      #likelihoods
      likelihoods <- compute.likelihoods(yij[i,], atoms, tau1[i], tau2[i], 
                                         lambda[i], omega[i], sigma2, t)
      
      #predictive weigths
      
      pnew <- a/(a+c) * theta/(theta + N - 1) 
      pzero <- (c/(a+c) * theta + nh[1])/(theta + N - 1) 
      pocc <- nh[-1]/(theta+N-1)
      predictive.weights <- c(pzero,pocc)
      if(singleton)  predictive.weights[hi] <- pnew
      if(!singleton) predictive.weights <- c(predictive.weights, pnew)

      #cat(predictive.weights)
      #cat("\n")
      #cat(likelihoods)
      #then update the cluser membership
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
        atoms <- atoms[populated]
        clusters[clusters>which(nh.new==0)] = clusters[clusters>which(nh.new==0)] - 1
      }
      else
      {
        populated <- nh.new>0
        atoms <- atoms[populated]
      }
      #reallocation for subject i done, restart for subject i+1
    }
    clusters
  }


# cluster allocation following Scarpa & Dunson

clusterallocation <- 
  function(yij, atoms, clusters, a, c, theta, sigma, C, 
           tau1, tau2, lambda, omega, sigma2, theta0, SigmaInv0)
  {  
    N <- nrow(yij)  
    #  cat("subject: ")
    for(i in 1:N)
    {
      #preprocess
      hi <- clusters[i]
      nh <- table(clusters)
      singleton <- nh[hi]==1
      nh[hi] <- nh[hi]-1 
      k <- length(nh)-singleton# length(table(clusters[-i]))
      n0 <- nh[1]
      #    cat("\xn\nk=", k, "with length(table(clusters[-i]))=", length(table(clusters[-i])))
      #    cat(", cluster = ", hi, "of (new) size ", nh[hi], "singleton=",singleton)
      #newatom
      if(!singleton) atoms[[k+1]] <- samplenewatom(theta0, SigmaInv0, knots)
      
      #likelihoods
      likelihoods <- compute.likelihoods(yij[i,], atoms, tau1[i], tau2[i], 
                                         lambda[i], omega[i], sigma2, t)
      
      #predictive weigths
      ldenom <- (1:n0)*log(c/(a+c)) + C[n0,1:n0] + 
        lgamma(theta/sigma + k - 1 + 1:n0) - lgamma(theta/sigma + k - 1)
      lnum1 <-  (1:n0)*log(c/(a+c)) + C[n0,1:n0] + 
        lgamma(theta/sigma + k     + 1:n0) - lgamma(theta/sigma + k)
      lnum2 <-  (1:(n0+1))*log(c/(a+c)) + C[(n0+1),1:(n0+1)] + 
        lgamma(theta/sigma + k - 1 + 1:(n0+1)) - lgamma(theta/sigma + k - 1)
      M <- max(max(ldenom), max(lnum1),max(lnum2))
      pnew <- a/(a+c) * (theta + (k-1)*sigma)/(theta + N-1) * sum(exp(lnum1-M))/sum(exp(ldenom-M))
      pzero <- 1/(theta + N-1) *  sum(exp(lnum2-M))/sum(exp(ldenom-M))
      pocc <- (nh[-1]-sigma)/(theta+N-1)
      predictive.weights <- c(pzero,pocc)
      if(singleton)  predictive.weights[hi] <- pnew
      if(!singleton) predictive.weights <- c(predictive.weights, pnew)
      #cat(predictive.weights)
      #cat("\n")
      #cat(likelihoods)
      #then update the cluser membership
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
        atoms <- atoms[populated]
        clusters[clusters>which(nh.new==0)] = clusters[clusters>which(nh.new==0)] - 1
      }
      else
      {
        populated <- nh.new>0
        atoms <- atoms[populated]
      }
      #reallocation for subject i done, restart for subject i+1
    }
    clusters
  }