# include <RcppArmadillo.h>
# include <omp.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include "helpers.h"

using namespace Rcpp;

const double log2pi = std::log(2.0 * M_PI);
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double std_normal_pdf(const double y)
{
  return(exp(-.5*y*y - .5*log2pi));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double normal_pdf(const double y, const double theta,const double sigma)
{
  return(std_normal_pdf((y-theta)/sigma)/sigma);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_normal_pdf(const double y, const double theta, const double sigma)
{
  double z = (y-theta)/sigma;
  return( -.5*z*z - .5*log2pi - log(sigma));
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double gamma_pdf(const double y, const double a, const double b )
{
  return( pow(y, a-1)*exp(-b*y)*pow(b,a)/tgamma(a) );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double log_gamma_pdf(const double y, const double a, const double b )
{
  return( (a-1)*log(y) -b*y + a*log(b) - lgamma(a) );
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec normal_pdf_arma(const arma::vec y, const arma::vec theta, const arma::vec sigma,
                          const bool sum=true, const bool log=true)
{
  int n = y.n_elem;
  arma::vec logout(n); logout.fill(0);   arma::vec out(n); out.fill(0);
  arma::vec sum_out(1); sum_out.fill(0);
  for( int i=0; i <n; i++)
  {
    logout(i) += log_normal_pdf(y(i),theta(0),sigma(0));
    out(i) += normal_pdf(y(i),theta(0),sigma(0));
    sum_out += logout(i);
  }
  if( sum & log ) { return(sum_out); }
  if( !sum & log ) { return(logout); }
  if( sum & !log ) { return(exp(sum_out)); }
  if( !sum & !log ) { return(out); }
  return(0);
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec gamma_pdf_arma(const arma::vec y, const double a, const double b,
                          const bool sum=true, const bool log=true)
{
  int n = y.n_elem;
  arma::vec logout(n); logout.fill(0);   arma::vec out(n); out.fill(0);
  arma::vec sum_out(1); sum_out.fill(0);
  for( int i=0; i <n; i++)
  {
    logout(i) += log_gamma_pdf(y(i),a,b);
    out(i) += gamma_pdf(y(i),a,b);
    sum_out += logout(i);
  }
  if( sum & log ) { return(sum_out); }
  if( !sum & log ) { return(logout); }
  if( sum & !log ) { return(exp(sum_out)); }
  if( !sum & !log ) { return(out); }
  return(0);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec rnormArma(int n, arma::vec mu, arma::vec sigma) {
  arma::vec Y = arma::randn(n, 1);
  arma::vec error(1); error(0) = -INFINITY;
  if( (1 == mu.n_elem) &  (1 ==sigma.n_elem) )
  {
    return arma::repelem(mu,1,n).t() + Y%arma::repelem(sigma,1,n).t();
    
  }
  if ( (n == mu.n_elem) &  (n ==sigma.n_elem) )
  {
    return mu + Y%sigma;
  }
  else { return(error); }
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat allocationprobs_univ_normal(const arma::vec y, const arma::vec theta, const arma::vec sigma,
                                   const arma::vec eta)
{
  int k = theta.n_elem;
  arma::mat out_Prob(y.n_elem,k); out_Prob.fill(-9999);
  if( theta.n_elem != eta.n_elem) {Rcout << "length mistmatch." << std::endl; }
  for( int i=0; i < y.n_elem; i++)
  {
    double normj = -INFINITY;
    arma::rowvec out_prob(k); 
    
    for( int j=0; j < k; j++)
    {
      out_prob(j) = log(eta(j)) + log_normal_pdf(y(i), theta(j),sigma(j));
      normj = logsumexp(normj, out_prob(j));
    }
    out_Prob.row(i) = exp(out_prob - normj );
  }
  return( out_Prob);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec mixture_likelihoods_univ_normal(const arma::vec y, const arma::vec theta, const arma::vec sigma,
                                   const arma::vec eta, const bool sum=true, const bool logt=true)
{
  int k = theta.n_elem;
  arma::vec out_vec(y.n_elem); out_vec.fill(-9999);
  arma::vec out_vec_exp(y.n_elem); out_vec_exp.fill(-9999);
  
  arma::vec out_sum(1); out_sum(0) = 0;
  if( theta.n_elem != eta.n_elem) {Rcout << "length mistmatch." << std::endl; }
  for( int i=0; i < y.n_elem; i++)
  {
    double normj = -INFINITY;
    arma::rowvec out_temp(k); 
    
    for( int j=0; j < k; j++)
    {
      out_temp(j) = log(eta(j)) + log_normal_pdf(y(i), theta(j),sigma(j));
      normj = logsumexp(normj, out_temp(j));
    }
    out_vec(i) = normj;
    out_vec_exp(i) = exp(normj);
    out_sum += normj;
  }
if(sum & logt)
{
 return(out_sum);
}
if(!sum & logt){
 return(out_vec);
}
if( sum & !logt)
{
  return(0);
}
if( !sum & !logt)
{
  return(out_vec_exp);
}
return(0);

}
