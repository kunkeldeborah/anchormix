// conjugate_normal.cpp: functions needed for updates of a conjugate GMM. 
//The component means thetaj are N(mu,sigma2j/kappa) and sigma2j is gamma(a,b).  
//Gamma is parameterized to have mean a/b. */ 

# include <RcppArmadillo.h>
# include <omp.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include "helpers.h"
# include "distributions.h"


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec elppd(arma::vec Y, arma::mat theta,
                     arma::mat sigma, arma::mat eta)
{
  int n_samples = theta.n_rows;
  int k = eta.n_cols;
  int n = Y.n_elem;
  arma::vec out_sum(1); out_sum(0) = 0;
  arma::vec out(n); out.fill(0);
  arma::vec prec(k);
  for( int t=0; t<n_samples; t++)
  {

  out = out + mixture_likelihoods_univ_normal(Y,theta.row(t).t(),sigma.row(t).t(),eta.row(t).t(), false, false) ;
    
  }
  
  for( int i=0; i<n; i++)
  {
    out(i) = log(out(i)/n_samples);
    out_sum += out(i);
  }
  
  return(out_sum);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec elppd_M(arma::mat Y, arma::mat theta,
                arma::mat sigma, arma::mat eta)
{

  int N = Y.n_rows;
  arma::vec out(1); out(0)=0;
  for( int h=0; h<N; h++)
  {
    
    out = out + elppd(Y.row(h).t(),theta,sigma,eta)/(N*1.0) ;
    
  }
  

  return(out);
}

