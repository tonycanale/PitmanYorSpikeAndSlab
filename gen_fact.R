# computes a nxn matrix the log of the generalized factorial coefficients,
# namely the element in position (i,j) of the matrix gen_fact_coeff
# coincides with log(C(i,j,sigma)). 

gen_fact_coeff <- function(n, sigma)
{
	M <- matrix(NA, n, n)
	M[1,1] <- log(sigma)
	for(i in 2:n)
	{
		M[i,1] <- log(-(sigma-(i-1))) + M[i-1,1]
		for(j in 2:i)
		{ 
			if (i==j) {M[i,j] <- log(sigma) + M[i-1,j-1]}
            else
			{
		        a <- log(sigma) + M[i-1,j-1];
		        b <- log(-( (sigma*j) - (i-1))) + M[i-1,j];
		        M[i,j] <- max(a,b) + log( 1 + exp(min(a,b ) - max(a,b)));
			}

		}
	}
	M[upper.tri(M)] <- t(M)[upper.tri(M)]
	return(M)
}
