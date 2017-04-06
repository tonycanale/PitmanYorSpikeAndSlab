#============================================================
# Spike & slab project
# Splines auxiliary functions
# May 2016
#============================================================
XtransposeX <- function(x) t(x) %*% x
Mscalar <- function(i,M,scalar) M[,,i]*scalar[i]
applyBytau2 <- function(i,B,y,tau1,tau2) 
{
  keep = !is.na(y[i,])
  tau2[i] * t(B[keep,,i]) %*% matrix(y[i,keep] - tau1[i])
}
#
betasampling <- function(yij, nh, tau1, tau2, lambda, omega, sigma2, SigmaInv0, beta0=0, domain=t, knots)
{
  d <- length(knots)+3
  if(nh>1)
  {
  tij <- t(matrix(rep(domain, nh), length(domain), nh))
  tijnorm <- (tij - lambda)/omega
  Bij <- array(apply(tijnorm, 1, bs, knots=knots),dim=c(length(domain),d,nh))
  Bijtau2 <- array(sapply(1:nh, FUN=Mscalar, M=Bij, scalar=tau2), dim=c(length(domain),d,nh))
  BtB <- array(apply(Bijtau2, 3, XtransposeX),dim=c(d,d,nh))
  sumBtB <- apply(BtB,c(1,2),sum)
  sumBytau2 <- rowSums(sapply(1:nh, FUN=applyBytau2, B=Bij, y=yij, tau1=tau1, tau2=tau2))
  
  Vinv <- SigmaInv0 + sumBtB/sigma2
  V <- solve(Vinv)
  m <- V %*% (SigmaInv0 %*% beta0 + sumBytau2/sigma2)
  newbeta <- rmnorm(1, m, V)
  }
  if(nh==1)
  {
    #sistemare la cosa quando un cluster e' con un singolo utente
    tij <- domain
    tijnorm <- (tij - lambda)/omega
    Bij <- bs(tijnorm, knots=knots)
    Bijtau2 <- Bij*tau2
    BtB <- t(Bijtau2) %*% Bijtau2 
    keep <- !is.na(yij)
    Bytau2 <- tau2 * t(Bij[keep,]) %*% matrix(yij[keep] - tau1)
        
    Vinv <- SigmaInv0 + BtB/sigma2
    V <- solve(Vinv)
    m <- V %*% (SigmaInv0 %*% beta0 + Bytau2/sigma2)
    newbeta <- rmnorm(1, m, V)
  }
  newbeta
}

#sample new atoms from the posterior conditionally on the cluster allocations 
reshuffle <- function(ind, clusters, yij, tau1, tau2, lambda, omega, 
                      sigma2, SigmaInv0, beta0=0, domain=t, knots)
{
  keep <- which(clusters==ind)
  resampledbeta <- betasampling(yij[keep,], length(keep), tau1[keep], tau2[keep], 
                                  lambda[keep], omega[keep], sigma2, 
                                  SigmaInv0, beta0, domain, knots)
  f <- function(x) bs(x, knots=knots) %*% resampledbeta
  f
}


