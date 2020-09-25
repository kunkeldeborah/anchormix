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
arma::vec update_theta_indprior(const arma::vec cvector, const arma::vec Y, const double precision, const double kappa, const double mu0)
{ 
  const int n = Y.n_elem;
  arma::vec out(1);
  double sumCY = 0; double sumC = 0;
  for( int i=0; i < n; i ++)
  {
    sumCY += Y(i)*cvector(i);
    sumC += cvector(i);
  }
  out(0) = (mu0*kappa + precision*sumCY)/(precision*sumC + kappa);
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_precision_indprior(const arma::vec cvector, const arma::vec Y, const double theta, const double a0, const double b0)
{ 
  const int n = Y.n_elem;
  arma::vec out(1);
  double sumCsq = 0; double sumC = 0;
  for( int i=0; i < n; i ++)
  {
    sumCsq += (Y(i)-theta)*(Y(i)-theta)*cvector(i);
    sumC += cvector(i);
  }
  out(0) = (a0+.5*sumC-1)/(b0 + .5*sumCsq);
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec update_b0_gammaprior(const double a0, const double g0, 
                    const double h0, const arma::vec precision)
{ 
  const int n_components = precision.n_elem;
  arma::vec out(1);
  double sum_p = 0;
  for( int j=0; j < n_components; j ++)
  {
    sum_p += precision(j);
  }
  out(0) = (a0*n_components+g0-1)/(h0 + sum_p);
  return(out);
}



// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec lb_ind_MAP_b(arma::vec Y, arma::vec theta,
                     arma::vec sigma, arma::vec eta, const arma::vec kappa, const arma::vec mu0,
                                const double a0, const double b0,
                                const double g0, const double h0)
{
  arma::vec out(1);
  int k = eta.n_elem;
  arma::vec prec(k);
  for( int j=0; j<k; j++)
  {
    prec(j) = pow(sigma(j),-2);
  }
  out(0) = mixture_likelihoods_univ_normal(Y,theta,sigma,eta, true)(0) + normal_pdf_arma(theta, mu0,pow(1.0/kappa,.5), true, true)(0) + 
  gamma_pdf_arma(prec, a0, b0, true,true)(0) + log_gamma_pdf(b0, g0, h0);
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec lb_ind_MAP_b_A(arma::vec Y, arma::vec theta,
                       arma::vec sigma, arma::vec eta, const arma::vec kappa, const arma::vec mu0,
                       const double a0, const double b0,
                       const double g0, const double h0, const arma::mat Cm, const arma::mat Am)
{
  arma::vec out(1);
  int k = eta.n_elem;
  arma::vec prec(k);
  double penalty = 0;
  int Acols = Am.n_cols;
  
  for( int j=0; j<k; j++)
  {
    prec(j) = pow(sigma(j),-2);
    for (int i=0; i<Acols; i++)
    {
      penalty += log(Cm(Am(j,i)-1,j));
    }
  }
  out(0) = mixture_likelihoods_univ_normal(Y,theta,sigma,eta, true)(0) + normal_pdf_arma(theta, mu0,pow(1.0/kappa,.5), true, true)(0) + 
    gamma_pdf_arma(prec, a0, b0, true,true)(0) + log_gamma_pdf(b0, g0, h0) + penalty;
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec lb_ind_MAP(arma::vec Y, arma::vec theta,
                       arma::vec sigma, arma::vec eta, const arma::vec kappa, const arma::vec mu0,
                       const double a0, const double b0)
{
  arma::vec out(1);
  out(0) = mixture_likelihoods_univ_normal(Y,theta,sigma,eta, true)(0) + normal_pdf_arma(theta, mu0,pow(1.0/kappa,.5), true, true)(0) + 
    gamma_pdf_arma(sigma, a0, b0, true,true)(0);
  return(out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec lb_ind_MAP_A(arma::vec Y, arma::vec theta,
                     arma::vec sigma, arma::vec eta, const arma::vec kappa, const arma::vec mu0,
                     const double a0, const double b0, const arma::mat Cm, const arma::mat Am)
{
  arma::vec out(1);
  double penalty = 0;
  int Acols = Am.n_cols;
  int k = eta.n_elem;
  
  for( int j=0; j<k; j++)
  {
  for (int i=0; i<Acols; i++)
  {
    penalty += log(Cm(Am(j,i)-1,j));
  }
  }
  out(0) = mixture_likelihoods_univ_normal(Y,theta,sigma,eta, true)(0) + normal_pdf_arma(theta, mu0,pow(1.0/kappa,.5), true, true)(0) + 
    gamma_pdf_arma(sigma, a0, b0, true,true)(0);
  return(out);
}
