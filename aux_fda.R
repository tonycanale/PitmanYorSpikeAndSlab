sampletaus <- function(ind, fij, temp, OmegaInv, sigma2, alpha_i, ID)
{
  woman <- ID[ind]
  yij <- temp[ind,]
  keep <- !is.na(yij)
  yij <- matrix(yij[keep])
  xij <- fij[ind,]
  xij <- matrix(xij[keep])
  xij <- cbind(1,xij)
  alpha <- matrix(alpha_i[woman,])
  XtX <- t(xij) %*% xij
  Xty <- t(xij) %*% yij
  V = solve(OmegaInv + XtX/sigma2)
  a = V %*% (OmegaInv %*% alpha + (Xty/sigma2))
  taus <- rmnorm(n = 1, mean = a, varcov=V)
  taus
}

samplelambdas <- function(ind, temp, t, tau1, tau2, omega, sigma2, cluster, atoms)
{
  #conditional posterior for lambda just the likelihook at discrete values
  f <- atoms[[ cluster[ind] ]]
  yij <- temp[ind,]
  w <- omega[ind]
  t1 <- tau1[ind]
  t2 <- tau2[ind]
  keep <- !is.na(yij)
  yij <- yij[keep]
  tij <- t[keep]
  values <- min(which(keep)):max(which(keep))
  tildefij <- matrix(NA, length(values), sum(keep))
  logplambda <- rep(1, length(values))
  for(l in 1:length(values))
  {
    tildefij[l, ] <- t1 + t2*f((tij-values[l])/w)
    logplambda[l] <- sum(dnorm(yij, tildefij[l, ], sqrt(sigma2), log=TRUE))
  }
  plambda <- exp(logplambda)/sum(exp(logplambda))
  sample(values, 1, prob=plambda, replace=TRUE)
}

sampleomegas <- function(ind, temp, t, a=1/2, b=1, tau1, tau2, omegaprec, lambda, sigma2, cluster, atoms)
{
  #conditional posterior for omega
  f <- atoms[[ cluster[ind] ]]
  yij <- temp[ind,]
  wprec <- omegaprec[ind]
  lamb <- lambda[ind]
  t1 <- tau1[ind]
  t2 <- tau2[ind]
  keep <- !is.na(yij)
  yij <- yij[keep]
  tij <- t[keep]
  wstar <- rgamma(1, a, b) #prior=proposal --> likratio
  tildefij_old  <- t1 + t2*f((tij-lamb)/wprec)
  tildefij_star <- t1 + t2*f((tij-lamb)/wstar)
  loglikold <-  sum(dnorm(yij, tildefij_old,  sqrt(sigma2), log=TRUE)) 
  loglikstar <- sum(dnorm(yij, tildefij_star, sqrt(sigma2), log=TRUE)) 
  loglikratio <- loglikstar - loglikold
  if(loglikratio>=0) newval <- wstar
  if(loglikratio<0)  newval <- sample(c(wstar, wprec), 1, prob=c(exp(loglikratio), 1-exp(loglikratio)))
  newval
}
#
samplealphas <- function(ind, ni, tau1, tau2, OmegaInv, R, alpha, ID)
{
  t1sum <- sum(tau1[ID==ind])
  t2sum <- sum(tau2[ID==ind])
  taus <- matrix(c(t1sum,t2sum),2,1)
  V2 = solve(R + ni[ind]*OmegaInv)
  a2 = V2 %*% (R %*% alpha + OmegaInv %*% taus)
  keep.on.rocking <- 1
  while(keep.on.rocking) #in the free world!
  {
  resval <- rmnorm(n = 1, mean = a2, varcov=V2)
  keep.on.rocking <- sum(resval>0)<2
  }
  resval
}
  
# calculate the new fij
applyfij <- function(i, tau1, tau2, lambda, omega, clusters, atoms, domain)
{
  hi <- clusters[i]
  f <- atoms[[hi]]
  tau1[i] + tau2[i]*f((domain-lambda[i])/omega[i])
}
