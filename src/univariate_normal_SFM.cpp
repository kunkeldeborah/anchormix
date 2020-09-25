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
double log_priordensity_gamma_indprior(const arma::vec gamma, const double kappa, const double mu0,
                                const double a0, const double b0 )
{ 
  double out = 0;
  out += log_normal_pdf(gamma(0), mu0, pow(kappa,-.5)) + log_gamma_pdf(gamma(1),a0,b0);
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec sample_gamma_indprior(const double kappa, const double mu0,
                                const double a0, const double b0 )
{
  int n_params = 2;
  arma::vec sample(n_params);
  arma::vec mu0vec(1); mu0vec(0) = mu0;
  arma::vec sigma(1); sigma(0) = pow(kappa,-.5);
  sample(0) = rnormArma(1,mu0vec,sigma)(0);
  sample(1) = R::rgamma(a0, 1/b0);
  return(sample);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sample_gamma_indprior_mat(const int nsamples, const double kappa, const double mu0,
                                const double a0, const double b0 )
{
  int n_params = 2;
  arma::mat sample(nsamples,n_params);
  arma::vec mu0vec(1); mu0vec(0) = mu0;
  arma::vec sigma(1); sigma(0) = pow(kappa,-.5);
  sample.col(0) = rnormArma(nsamples,mu0vec,sigma);
  for( int i =0 ; i < nsamples; i ++)
  {
  sample(i,1) = R::rgamma(a0, 1/b0);
  }
  return(sample);
}


 // [[Rcpp::depends(RcppArmadillo)]]
 // [[Rcpp::export]]
 arma::vec update_gammaj_indprior_SFM(const int component, const arma::vec y, const arma::vec s,
                                 const arma::vec gamma, const double kappa, const double mu0,
                                 const double a0, const double b0)
 {
   // Gamma is a vector with 1st element muj and 2nd element 1/sigmaj2
   // this is needed for the full conditional of mj | sigmaj
   int n = y.n_elem; int n_c = 0; double n_cdouble=0;
   arma::vec error(1); error.fill(-INFINITY);
   if( s.n_elem != n ) { return(error); }
   int n_params = 2;
   arma::vec gamma_out(n_params); arma::vec yj(n);
   arma::vec mean_muj(1); arma::vec std_muj(1);
   arma::vec new_muj(1); arma::vec new_sigma2inv(1);
   arma::vec sumj(1); sumj.fill(0);
   
   for(int i=0; i < n; i ++)
   {
     if( s(i)==component) 
       { 
        yj(n_c) = y(i); n_c += 1;  
        sumj(0) += y(i);
       }
   }
  
   n_cdouble = n_c;
   mean_muj(0) = (gamma(1)*sumj(0) +  kappa*mu0)/(gamma(1)*n_cdouble + kappa);
   std_muj(0) = pow(gamma(1)*n_cdouble + kappa, -.5);
   new_muj(0) = R::rnorm(mean_muj(0), std_muj(0));
   gamma_out(0) = new_muj(0);
   
   double sumsq = 0;
   
   for( int i =0; i < n_c; i ++) { sumsq += (yj(i)-new_muj(0))*(yj(i)-new_muj(0)); }
   
   new_sigma2inv = R::rgamma(a0 +.5*n_cdouble, 1.0/(b0 + .5*sumsq));
   gamma_out(1) = new_sigma2inv(0); 
   return(gamma_out);
   
 }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_gamma_indprior_SFM(const int K, const arma::vec y, const arma::vec s,
                                const arma::mat gamma, const double kappa, const double mu0,
                                const double a0, const double b0)
{
  // Gamma is a vector with 1st element muj and 2nd element 1/sigmaj2
  int n = y.n_elem; 
  arma::vec nvec(K); nvec = nvec_count(s, K);
  arma::vec error(1); error.fill(-INFINITY);
  if( s.n_elem != n ) { return(error); }
  int n_params = 2;
  arma::mat GammaMat_out(K,n_params);

    for(int j=0; j < K; j ++)
  {
      int nj = nvec(j);
    if(nj > 0 )
    {
      GammaMat_out.row(j) = update_gammaj_indprior_SFM(j+1,y, s, gamma.row(j).t(),kappa,mu0,a0,b0 ).t();
    }
    if( nj== 0)
    {
      GammaMat_out.row(j) = sample_gamma_indprior(kappa, mu0, a0, b0).t();
    
    }
  }
  
  return(GammaMat_out);
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double exp_inf(double x)
{
  return(exp(-INFINITY));
}

/* next functions are for conjugate prior */


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec sample_gamma_conjprior(const double kappa, const double mu0,
                                const double a0, const double b0 )
{
  int n_params = 2;
  arma::vec sample(n_params);
  arma::vec mu0vec(1); mu0vec(0) = mu0;
  sample(1) = R::rgamma(a0, 1/b0);
  sample(0) = R::rnorm(mu0vec(0),pow(kappa,-.5)/pow(sample(1),.5));
  return(sample);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sample_gamma_conjprior_mat(const int nsamples, const double kappa, const double mu0,
                                    const double a0, const double b0 )
{
  int n_params = 2;
  arma::mat sample(nsamples,n_params);
  arma::vec mu0vec(nsamples); mu0vec.fill(-9e3);
  arma::vec sigma(nsamples);
  double prior_std = pow(kappa,-.5);
  for( int i =0 ; i < nsamples; i ++)
  {
    sample(i,1) = R::rgamma(a0, 1/b0);
    mu0vec(i) = mu0;
    sigma(i) = prior_std*pow(sample(i,1),-.5);
  }
  sample.col(0) = rnormArma(nsamples,mu0vec,sigma);
  return(sample);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_gammaj_conjprior_SFM(const int component, const arma::vec y, const arma::vec s,
                                     const arma::vec gamma, const double kappa, const double mu0,
                                     const double a0, const double b0)
{
  // Gamma is a vector with 1st element muj and 2nd element 1/sigmaj2
  // this is needed for the full conditional of mj | sigmaj
  int n = y.n_elem; int n_c = 0; double n_cdouble=0;
  arma::vec error(1); error.fill(-INFINITY);
  if( s.n_elem != n ) { return(error); }
  int n_params = 2;
  arma::vec gamma_out(n_params); arma::vec yj(n);

  double sumj = 0; 
  double ssx = 0;
  
  for(int i=0; i < n; i ++)
  {
    if( s(i)==component) 
    { 
      yj(n_c) = y(i); n_c += 1;  
      sumj += y(i);
      ssx += y(i)*y(i);
    }
  }
  n_cdouble = n_c;
  
  
  double mun = (kappa*mu0+sumj)/(kappa + n_cdouble);
  double kappan = kappa + n_cdouble;
  double an = a0 + n_cdouble/2.0;
  double bn = b0 + .5*(ssx - sumj*sumj/n_cdouble) + (kappa*n_cdouble*(sumj/n_cdouble - mu0)*(sumj/n_cdouble - mu0))/(2*(kappa+ n_cdouble));
  if( n_cdouble == 0) { bn = b0; mun = mu0; kappan = kappa; }
  gamma_out = sample_gamma_conjprior(kappan,mun, an, bn );

  return(gamma_out);
  
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_gamma_conjprior_SFM(const int K, const arma::vec y, const arma::vec s,
                                    const arma::mat gamma, const double kappa, const double mu0,
                                    const double a0, const double b0)
{
  // Gamma is a vector with 1st element muj and 2nd element 1/sigmaj2
  int n = y.n_elem; 
  arma::vec nvec(K); nvec = nvec_count(s, K);
  arma::vec error(1); error.fill(-INFINITY);
  if( s.n_elem != n ) { return(error); }
  int n_params = 2;
  arma::mat GammaMat_out(K,n_params);
  
  for(int j=0; j < K; j ++)
  {
    int nj = nvec(j);
    if(nj > 0 )
    {
      GammaMat_out.row(j) = update_gammaj_conjprior_SFM(j+1,y, s, gamma.row(j).t(),kappa,mu0,a0,b0 ).t();
    }
    if( nj== 0)
    {
      GammaMat_out.row(j) = sample_gamma_conjprior_mat(1, kappa, mu0, a0, b0).row(0);
      
    }
  }
  
  return(GammaMat_out);
  
}

/* uniform Rlow,Rup prior on mean, gamma a0, b0 prior on precision. */

 // [[Rcpp::depends(RcppArmadillo)]]
 // [[Rcpp::export]]
 arma::vec sample_gamma_unifprior(const double Rlow, const double Rup,
                                 const double a0, const double b0 )
 {
   int n_params = 2;
   arma::vec sample(n_params);
   sample(0) = R::runif(Rlow,Rup);
   sample(1) = R::rgamma(a0, 1/b0);
   return(sample);
 }

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sample_gamma_unifprior_mat(const int nsamples, const double Rlow, const double Rup,
                                    const double a0, const double b0 )
{
  int n_params = 2;
  arma::mat sample(nsamples,n_params);
  for( int i =0 ; i < nsamples; i ++)
  {
    sample(i,0) = R::runif(Rlow,Rup);
    sample(i,1) = R::rgamma(a0, 1/b0);
  }
  return(sample);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_gammaj_unifprior_SFM(const int component, const arma::vec y, const arma::vec s,
                                     const arma::vec gamma, const double Rlow, const double Rup,
                                     const double a0, const double b0)
{
  // Gamma is a vector with 1st element muj and 2nd element 1/sigmaj2
  // this is needed for the full conditional of mj | sigmaj
  int n = y.n_elem; int n_c = 0; double n_cdouble=0;
  arma::vec error(1); error.fill(-INFINITY);
  if( s.n_elem != n ) { return(error); }
  int n_params = 2;
  arma::vec gamma_out(n_params); arma::vec yj(n);
  arma::vec mean_muj(1); arma::vec std_muj(1);
  arma::vec new_muj(1); arma::vec new_precision(1);
  arma::vec sumj(1); sumj.fill(0);
  
  for(int i=0; i < n; i ++)
  {
    if( s(i)==component) 
    { 
      yj(n_c) = y(i); n_c += 1;  
      sumj(0) += y(i);
    }
  }
  if( n_c > 0)
  {
  
  n_cdouble = n_c;
  mean_muj(0) = sumj(0)/n_cdouble;
  std_muj(0) = pow(gamma(1)*n_cdouble, -.5);
  new_muj(0) = R::rnorm(mean_muj(0), std_muj(0));
  gamma_out(0) = new_muj(0);
  
  double sumsqy = 0;
  double sq_ybarmu = 0;
  sq_ybarmu = pow(mean_muj(0) - new_muj(0), 2.0);
  for( int i =0; i < n_c; i ++) { sumsqy += (yj(i)-mean_muj(0))*(yj(i)-mean_muj(0)); }
  
  new_precision = R::rgamma(a0 +.5*n_cdouble, 1.0/(b0 + .5*sq_ybarmu*n_cdouble + sumsqy));
  gamma_out(1) = new_precision(0); 
  }
  if( n_c == 0)
  {
    gamma_out = sample_gamma_unifprior( Rlow, Rup, a0, b0 );
  }
  return(gamma_out);
  
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat update_gamma_unifprior_SFM(const int K, const arma::vec y, const arma::vec s,
                                    const arma::mat gamma, const double Rlow, const double Rup,
                                    const double a0, const double b0)
{
  // Gamma is a vector with 1st element muj and 2nd element 1/sigmaj2
  int n = y.n_elem; 
  arma::vec nvec(K); nvec = nvec_count(s, K);
  arma::vec error(1); error.fill(-INFINITY);
  if( s.n_elem != n ) { return(error); }
  int n_params = 2;
  arma::mat GammaMat_out(K,n_params);
  
  for(int j=0; j < K; j ++)
  {
    int nj = nvec(j);
    if(nj > 0 )
    {
      GammaMat_out.row(j) = update_gammaj_unifprior_SFM(j+1,y, s, gamma.row(j).t(),Rlow,Rup,a0,b0 ).t();
    }
    if( nj== 0)
    {
      GammaMat_out.row(j) = sample_gamma_unifprior(Rlow, Rup, a0, b0).t();
      
    }
  }
  
  return(GammaMat_out);
  
}
 