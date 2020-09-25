#include <Rcpp.h>
#ifndef DISTRIBUTIONS
#define DISTRIBUTIONS

# include <RcppArmadillo.h>
# include <omp.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>


const double log2pi = std::log(2.0 * M_PI);

double std_normal_pdf(const double y);

double normal_pdf(const double y, const double theta,const double sigma);

double log_normal_pdf(const double y, const double theta, const double sigma);

double gamma_pdf(const double y, const double a, const double b );

double log_gamma_pdf(const double y, const double a, const double b );

arma::vec normal_pdf_arma(const arma::vec y, const arma::vec theta, const arma::vec sigma,
                          const bool sum=true, const bool log=true);

arma::vec gamma_pdf_arma(const arma::vec y, const double a, const double b,
                         const bool sum=true, const bool log=true);


arma::vec rnormArma(int n, arma::vec mu, arma::vec sigma); 

arma::mat allocationprobs_univ_normal(const arma::vec y, const arma::vec theta, const arma::vec sigma,
                                      const arma::vec eta);

arma::vec mixture_likelihoods_univ_normal(const arma::vec y, const arma::vec theta, const arma::vec sigma,
                                          const arma::vec eta, const bool sum=true, const bool logt=true);


#endif