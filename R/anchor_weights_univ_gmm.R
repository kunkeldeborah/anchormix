weights_univ_gmm <- function(y, Am, permM, theta_samp, sigsq_samp)
{
  pvec <- sapply(1:nrow(permM), function(q){ prod(dnorm(y[Am[,1]], theta_samp[permM[q,]][Am[,2]],
             sqrt(sigsq_samp)[permM[q,]][Am[,2]])) } ) 
  pvec/sum(pvec)
}

