#include <Rcpp.h>
#ifndef HELPERS
#define HELPERS

# include <RcppArmadillo.h>
# include <omp.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <iostream>

double logsumexp_vec(arma::vec avec);

double logsumexp(const double a, const double b);

arma::mat normalize_A(arma::mat A);

arma::mat normalize_vec(arma::vec v);

arma::mat normalize_logvec(arma::vec lv);

arma::vec nvec_count(arma::vec svec, int k);
 
arma::vec ordered_remove(arma::vec list, int index);

arma::vec ordered_insert_next(arma::vec list);

int ordered_next(arma::vec list);

bool in_vec_int(double n, arma::vec v);

int which_vec_int(int n, arma::vec v);

arma::vec int_seq(const int k);


#endif
