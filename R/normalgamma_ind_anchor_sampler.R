# normalgamma_ind_anchor_sampler.R
# function for Gibbs sampling of anchored finite mixture
# univariate with independent normal gamma prior


sampleA_indprior = function(Ydata, Am=NULL,Hyp.list, K, n_params=2, Nsamples=1e4, nburn = 1e3, nthin = 25, nchains=10,
                            perm.step=TRUE)
{
  nsamples = ceiling(Nsamples/nchains)*nchains
  nsamples.per.chain = ceiling(Nsamples/nchains)
  n = length(Ydata)
  gamma_samples = array(NA, dim=c(nsamples,K,n_params) )
  eta_samples = matrix(NA, nrow=nsamples, ncol = K)
  S_samples = matrix(NA, nrow=nsamples, ncol = n)
  sampled_perms <- NULL
  
  if( perm.step)
  {
    perm_matrix <- permutations(K)
  }
  
  if( !is.null(Am))
  {
  S_samples[,Am[,1]] = S_samples[,Am[,2]]
  }
  b0_samples = rep(NA, nsamples)
  accept_samples = rep(NA, nsamples);   prob_samples = rep(NA, nsamples)
  
  for( ch in 1:nchains)
  {
  it = 1; samp = 1;
  
  # initialize model parameters
  current_b0 = rgamma(1, Hyp.list$g0, Hyp.list$h0)
  current_eta = rdirichlet(1,rep(Hyp.list$alpha,K))
  current_S = sample(1:K, n, prob=current_eta, replace=T)
  if( !is.null(Am))
  {
  current_S = enforce_anchors(Am, current_S)
  }
  current_gamma = sample_gamma_indprior_mat(K, Hyp.list$kappa, Hyp.list$mu0, Hyp.list$a0, current_b0)
  
  while( samp <= nsamples.per.chain)
  {  
  current_gamma = update_gamma_indprior_SFM(K, Ydata, current_S, current_gamma, Hyp.list$kappa,Hyp.list$mu0,  Hyp.list$a0, current_b0)
  current_eta = rdirichlet(1, Hyp.list$alpha+ nvec_count(current_S,K))
  
  if( perm.step & is.null(Am))
  {
    perm <- sample(1:K, K, replace=F)
    current_gamma <- current_gamma[perm,]
    current_eta <- current_eta[perm]
  }
  
  if(perm.step & !is.null(Am) )
  {
    p.vec <- weights_univ_gmm(Ydata, Am, perm_matrix, current_gamma[,1], 1/current_gamma[,2])
    qq <- sample(1:nrow(perm_matrix), 1, prob=p.vec, replace=F)
    perm <- perm_matrix[qq,]
    current_gamma <- current_gamma[perm,]
    current_eta <- current_eta[perm]
  }
  current_S = update_allocations(Ydata, current_gamma[,1], 1/sqrt(current_gamma[,2]), current_eta)
  if( !is.null(Am))
  {
  current_S = enforce_anchors(Am, current_S)
  }
  current_b0 = update_b0_indprior_SFM(current_gamma, K, Hyp.list$a0, Hyp.list$g0,Hyp.list$h0)
  

  # save current sample ever nthin^th time.
  if(  it > (nburn - 1) & it%%nthin == 0 )
  {
    gamma_samples[nsamples.per.chain*(ch-1)+samp,,] = current_gamma
    eta_samples[nsamples.per.chain*(ch-1)+samp,] = current_eta
    S_samples[nsamples.per.chain*(ch-1)+samp, ] = current_S
    b0_samples[nsamples.per.chain*(ch-1)+samp] = current_b0
    if(perm.step & !is.null(Am))
    {
    sampled_perms <- c(sampled_perms,qq)
    }
    samp = samp+1;
  }
  it = it + 1
  
  }
  }

  
  return(list(gamma = gamma_samples, eta = eta_samples, S = S_samples, b0 = b0_samples))
}



