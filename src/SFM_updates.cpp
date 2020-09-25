# include <RcppArmadillo.h>
# include <omp.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <iostream>
# include <RcppArmadilloExtensions/sample.h>
# include "helpers.h"
# include "distributions.h"

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_sticks(const arma::vec Svec, const arma::vec avec, const arma::vec bvec)
{
  int K = avec.n_elem; int sumNkplus;
  arma::vec Nvec(K); Nvec= nvec_count(Svec, K);
  if( K != Nvec.n_elem ) { Rcout<< "problem:length Nvec not equal K." << std::endl; }
  arma::vec vout(K);
  vout(K-1) =1;
  sumNkplus = sum(Nvec)-Nvec(0);
  for( int k = 0; k < (K-1); k++)
  {
    vout(k) = R::rbeta(avec(k)+Nvec(k), bvec(k) + sumNkplus );
    sumNkplus += -Nvec(k+1);
  }
  
  return(vout);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec calculate_weights(const arma::vec sticks)
{
  int K = sticks.n_elem; int sumNkplus;
  arma::vec wout(K);
  wout(0) =sticks(0);
  double prod = (1-sticks(0));
  wout(1) = sticks(1)*prod;
  for( int k = 2; k < (K); k++)
  {
    prod = prod*(1-sticks(k-1));
    wout(k) = sticks(k)*prod;
  }
  
  return(wout);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_allocations(const arma::vec y, const arma::vec theta, const arma::vec sigma,
                             const arma::vec eta)
{
  arma::mat updateP(y.n_elem, eta.n_elem);
  arma::vec sout(y.n_elem);
  updateP = allocationprobs_univ_normal(y, theta, sigma, eta);
  for(int i =0 ; i<y.n_elem; i++ )
  {
    arma::vec probi(eta.n_elem);
    probi = updateP.row(i).t();
    sout(i) = RcppArmadillo::sample(int_seq(eta.n_elem), 1, false, probi)(0);
  }
  return(sout);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec enforce_anchors(const arma::mat A, arma::vec s)
{
  for(int i =0 ; i<A.n_rows; i++ )
  {
    s(A(i,0)-1) = A(i,1);
  }
  return(s);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_post_e0(double e0, const int K, const arma::vec Nvec, const double e0_ha0, const double e0_hb0)
{
  int Kplus = 0;
  double post = 0;
  for( int j=0; j < K; j ++)
  {
    if( Nvec(j)>0 )
    {
    Kplus += 1;
    post += lgamma(Nvec(j)+ e0) - (lgamma(e0 + 1) - log(e0));
    }
  }
  post += lgamma(K+1) - lgamma(K - Kplus +1) + lgamma(K*e0 +1) - log(K*e0) -lgamma(sum(Nvec)+ K*e0);
  post += log_gamma_pdf( e0, e0_ha0, e0_hb0);
  return(post);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double e0_accept_prob(double proposed_e0, double e0, const int K, const arma::vec Nvec, const double e0_ha0, const double e0_hb0)
{
  double log_prob = 0;
  log_prob = log_post_e0(proposed_e0, K, Nvec, e0_ha0, e0_hb0) - log_post_e0( e0, K, Nvec, e0_ha0, e0_hb0);
  return(exp(log_prob));
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List update_e0(const double current_e0, const int K, const arma::vec Svec, const double e0_ha0, const double e0_hb0,
                 const double prop_std = 5e-2)
{
  double prob=-9e3;
  double new_e0 = -9e3;
  double accept = -9e3;
  arma::vec Nvec(K); Nvec = nvec_count(Svec,K);
  double proposed_e0 = -9e3;
  proposed_e0 = current_e0 + arma::randn()*prop_std; 
  if( proposed_e0 <= 0 ) { new_e0 = current_e0; accept = 0;  }
  else
    {
    prob = e0_accept_prob(proposed_e0, current_e0, K, Nvec, e0_ha0, e0_hb0);
    double u = arma::randu();
    if( u > prob ) { new_e0 = current_e0; accept = 0; }
    if( u < prob ) { new_e0 = proposed_e0; accept = 1; }
    }
  return List::create(_["e0"] = new_e0, _["accept"] = accept,
                      _["prob"] = prob ); 
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int Kplus(const int K, const arma::vec Svec)
{
  int K_plus = 0;
  arma::vec n(K);
  n = nvec_count(Svec,K);
  for( int k = 0; k < K; k ++)
  {
    if(n(k) > 0) K_plus +=1;
  }
 return(K_plus);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double update_b0_indprior_SFM(const arma::mat gamma, const double k, const double a0, const double ghyp,
                              const double hhyp)
{
  double out = -INFINITY;
  double precisionsum = 0;
  for( int j=0; j < gamma.n_rows; j ++)
  {
    precisionsum += gamma(j,1);
  }
  out = R::rgamma(ghyp + k*a0,1.0/(hhyp + precisionsum));
  return(out);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double update_lambda_NG_SFM(const arma::mat gamma, const double k, const double mu0, const double R0, const double nu1,
                              const double nu2)
{
  double out = -INFINITY;
  double sumsq = 0;
  for( int j=0; j < gamma.n_rows; j ++)
  {
    sumsq += (gamma(j,0)-mu0)*(gamma(j,0) -mu0 ) ;
  }
  out = R::rgamma(nu1 + k/2.0,1.0/(nu2 + .5*sumsq/pow(R0,.5)));
  return(out);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double update_mu0_NG_SFM(const arma::mat gamma, const double k, const double R0, const double lambda,
                            const double m0,const double M0)
{
  double out = -INFINITY; double mun = 0; double vn = 0;
  double muj_sum = 0;
  for( int j=0; j < gamma.n_rows; j ++)
  {
    muj_sum += gamma(j,0);
  }
  mun = (lambda*R0*m0 + M0*muj_sum)/(lambda*R0 + k*M0);
  vn = (M0*lambda*R0)/( k*M0 + lambda*R0 );
  out = R::rnorm(mun,pow(vn, .5));
  return(out);
}
