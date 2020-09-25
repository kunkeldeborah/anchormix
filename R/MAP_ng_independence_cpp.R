# MAP_normal_wishart

# MAP estimation for two models:
# independence: theta ~ N(mu,kappa.inv*tausq); sigma.inv^2 ~ gamma(a,b)
# conjugate: theta ~ N(mu, kappa.inv*sigma.inv^2); sigma.inv^2 ~ gamma(a,b)

update_allocation_probs = function(Y, temp.eta, temp.theta, temp.sigma )
{
  n = length(Y)
  C_out = matrix(NA, nrow=n, ncol=length(temp.eta))
  for( i in 1:n)
  {
    C_out[i,] = sapply(1:length(temp.eta), function(j){ temp.eta[j]*dnorm(Y[i], temp.theta[j], temp.sigma[j])  })
  }
  C_out
}



MAP_indng_EM_cpp = function(Y.data, my.k, my.priorlist, n.iter = 5, diff.threshold=1e-4, update_b = FALSE)
{
  # check that b0 hyperparameters are given, if necessary
  if(update_b)
  {
    if( is.null(my.priorlist$g0) | is.null(my.priorlist$h0)) {stop("provide g0 and h0 if update_b=TRUE.")}
  }

  # check that b0 hyperparameters are given, if necessary
  if(!update_b)
  {
    if( is.null(my.priorlist$b0) ) {stop("provide b0 if update_b=FALSE.")}
  }
  # repeat EM for niter trials
  niter = n.iter
  tvec = rep(NA, niter)
  n = length(Y.data)
  k = my.k
  
  means.list= list() 
  precisions.list= list() 
  etas.list = list() 
  LowerB = list()
  LowerBmax = list()
  b.list <- list()
  
  for( iter in 1:niter)
  {
    LowerBound = -99999

    diff=100
    t=1
    thetat = NULL
    precisiont =NULL
    etat = rep(1/my.k,my.k)
    C.full = matrix(NA, nrow=n, ncol=my.k)
    
    # Initialize allocation vector and component means
    stemp = rmultinom(n, 1, etat); stemp = matrix(1:my.k, nrow=1)%*%(stemp)
    thetat = runif(my.k,min(Y.data),max(Y.data))

 
    if( !update_b)
    {
    b0_t = my.priorlist$b0
    }
    else
    {
      b0_t <- 2
    }

    precisiont = rgamma(my.k,my.priorlist$a0, b0_t)
 
    # uncomment these line to use deterministic initial values
    stemp = rep(1:my.k, ceiling(n/my.k))[1:n]; #thetat = tapply(Y.data,as.factor(stemp),mean) 
    precisiont = rep(1.3, my.k)
   
    # current parameter values
    thetat_temp = thetat
    precisiont_temp = precisiont
    etat_temp = etat
    
    while( abs(diff) > diff.threshold)
    #while(t < 2e3)
    {
      C.full = allocationprobs_univ_normal(Y.data, thetat_temp, 1/sqrt(precisiont_temp),etat_temp )


      # Update theta and sigma
      for( j in 1:k)
      {
        thetat_temp[j] = update_theta_indprior( C.full[,j], Y.data,  precisiont_temp[j],my.priorlist$kappa,
my.priorlist$mu0)
        precisiont_temp[j] = update_precision_indprior(C.full[,j], Y.data, thetat_temp[j], my.priorlist$a0, b0 = b0_t)

      }


      thetat = rbind(thetat , thetat_temp)
      precisiont = rbind(precisiont , precisiont_temp)
      
      
      # Update eta
      etat_temp = apply(C.full,2,sum) 
      etat_temp = ( etat_temp + my.priorlist$alpha - 1 )/( sum(etat_temp + my.priorlist$alpha - 1) )
      if ( sum(etat_temp)!=1 )
      {
        etat_temp = etat_temp / (sum(etat_temp))
      }
      etat = rbind(etat, etat_temp)
      
      if( update_b)
      {
        b0_t <- update_b0_gammaprior(my.priorlist$a0,my.priorlist$g0, my.priorlist$h0, precisiont_temp)
        LowerBound = c(LowerBound, lb_ind_MAP_b(Y.data, thetat_temp, 1/sqrt(precisiont_temp), etat_temp, 
		my.priorlist$kappa, my.priorlist$mu0, my.priorlist$a0, b0_t,  my.priorlist$g0,
		my.priorlist$h0)+log(ddirichlet(etat_temp, rep(my.priorlist$alpha,k)) ) )

      }

      else
      {
      LowerBound = c(LowerBound, lb_ind_MAP(Y.data, thetat_temp, 1/sqrt(precisiont_temp), etat_temp, 
		my.priorlist$kappa, my.priorlist$mu0, my.priorlist$a0, b0_t) )
      }
        t = t+1
        #diff = 1
        #if( t > 1e3) { diff = 1e-5 }
      diff = LowerBound[t] - LowerBound[t-1]
      
    }
    tvec[iter] = t
    means.list[[iter]] = thetat_temp
    precisions.list[[iter]] = precisiont_temp
    etas.list[[iter]] = etat_temp
    LowerB[[iter]] = LowerBound
    LowerBmax[[iter]] = LowerBound[t]
    b.list[[iter]] <- b0_t
    
  }
  
  return( list( theta.list = means.list, precisions.list = precisions.list, etas.list=etas.list, 
                LowerB = LowerB, LowerBmax = LowerBmax, b.list=b.list) )
}

# anchored EM for the same model
AEM_indng_cpp = function(Y.data, my.k, m.all, my.priorlist, n.iter = 5, diff.threshold=1e-4, update_b=FALSE)
{
  # check that b0 hyperparameters are given, if necessary
  if(update_b)
  {
    if( is.null(my.priorlist$g0) | is.null(my.priorlist$h0)) {stop("provide g0 and h0 if update_b=TRUE.")}
  }

  # check that b0 hyperparameters are given, if necessary
  if(!update_b)
  {
    if( is.null(my.priorlist$b0) ) {stop("provide b0 if update_b=FALSE.")}
  }
  # repeat EM for niter trials
  niter = n.iter
  tvec = rep(NA, niter)
  n = length(Y.data)
  k = my.k
  
  means.list= list() 
  precisions.list= list() 
  etas.list = list() 
  Anchors.list = list()
  LowerB = list()
  LowerBmax = list()
  b.list <- list()
  
  for( iter in 1:niter)
  {
    LowerBound = -99999
    
    diff=100
    t=1
    thetat = NULL
    precisiont =NULL
    Amt = NULL
    etat = rep(1/k,k)
    C.full = matrix(NA, nrow=n, ncol=k)
    
    if( !update_b)
    {
      b0_t = my.priorlist$b0
    }
    else
    {
      b0_t <- 2
    }
    
    # Initialize allocation vector and component means
    stemp = rmultinom(n, 1, etat); stemp = matrix(1:k, nrow=1)%*%(stemp)
    thetat = runif(k,min(Y.data),max(Y.data))
    precisiont = rgamma(k,my.priorlist$a0, b0_t)
    
    # current parameter values
    thetat_temp = thetat
    precisiont_temp = precisiont
    etat_temp = etat
    
    while( abs(diff) > diff.threshold)
    {
     C.full = allocationprobs_univ_normal(Y.data, thetat_temp, 1/sqrt(precisiont_temp),etat_temp )
      
      Atemp <- sort_probmatrix(melt_probmatrix(C.full),0)
      Anchors <- find_anchored(Atemp, rep(m.all,k),k)
      #Anchors2 = matrix(NA, ncol=m.all, nrow=k)
      # select mj anchor points based on C.full
      #for( r in 1:m.all)
      #{
        #sort indices of mth biggest C value for each component
      #  order.temp = order(sapply(1:k, function(j){C.full[order(C.full[,j], decreasing=T)[r],j]}), decreasing=T)
      #  for( j in order.temp)
      #  {
          # for components ordered by order.temp
          # vector of matches (logical) between ordered probs for jth component
          # and vector of points already anchored
      #    match.vec = match( order(C.full[,j], decreasing=T),c( Anchors[order.temp,1:r]))
      #    Anchors[j,r]= order(C.full[,j], decreasing=T)[which( is.na(match.vec ))[1]]
      #  }
      #}
      
      # Anchored matrix of allocation probabilities
      C.anch = C.full
      for( j in 1:k)
      {
        C.anch[Anchors[j,],c(1:k)[-j]]=0
        C.anch[Anchors[j,],j]=1
      }
      # store anchors in array
      Amt = abind(Amt, Anchors, along=3)
      
      # Update theta and sigma
      for( j in 1:k)
      {
        thetat_temp[j] = update_theta_indprior( C.anch[,j], Y.data,  precisiont_temp[j],my.priorlist$kappa,
my.priorlist$mu0)
        precisiont_temp[j] = update_precision_indprior(C.anch[,j], Y.data, thetat_temp[j], my.priorlist$a0, b0 = b0_t)

      }
      thetat = rbind(thetat , thetat_temp)
      precisiont = rbind(precisiont , precisiont_temp)
      
      
      # Update eta
      etat_temp = apply(C.anch,2,sum) 
      etat_temp = ( etat_temp + my.priorlist$alpha - 1 )/( sum(etat_temp + my.priorlist$alpha - 1) )
      if ( sum(etat_temp)!=1 )
      {
        etat_temp = etat_temp / (sum(etat_temp))
      }
      etat = rbind(etat, etat_temp)
      
      if( update_b)
      {
        b0_t <- update_b0_gammaprior(my.priorlist$a0,my.priorlist$g0, my.priorlist$h0, precisiont_temp)
        LowerBound = c(LowerBound, lb_ind_MAP_b_A(Y.data, thetat_temp, 1/sqrt(precisiont_temp), etat_temp, 
		my.priorlist$kappa, my.priorlist$mu0, my.priorlist$a0, b0_t,  my.priorlist$g0,
		my.priorlist$h0, C.full, Anchors)+log(ddirichlet(etat_temp, rep(my.priorlist$alpha,k)) ) )

      }

      else
      {
      LowerBound = c(LowerBound, lb_ind_MAP_A(Y.data, thetat_temp, 1/sqrt(precisiont_temp), etat_temp, 
		my.priorlist$kappa, my.priorlist$mu0, my.priorlist$a0, b0_t, C.full, Anchors) )
      }
      t = t+1
      diff = LowerBound[t] - LowerBound[t-1]
      
    }
    tvec[iter] = t
    means.list[[iter]] = thetat_temp
    precisions.list[[iter]] = precisiont_temp
    etas.list[[iter]] = etat_temp
    Anchors.list[[iter]] = Anchors
    LowerB[[iter]] = LowerBound
    LowerBmax[[iter]] = LowerBound[t]
    b.list[[iter]] <- b0_t
    
  }
  
  return( list( theta.list = means.list, precisions.list = precisions.list, etas.list=etas.list, 
                Anchors.list = Anchors.list,
                LowerB = LowerB, LowerBmax = LowerBmax, b.list=b.list) )
}
