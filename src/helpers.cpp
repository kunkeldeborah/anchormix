# include <RcppArmadillo.h>
# include <omp.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <iostream>
using namespace Rcpp;



/* returns log(sum(exp(v))) */
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double logsumexp_vec(arma::vec avec)
{
  double m = -999;
  double expsum = 0;
  m = max(avec);
  int n = avec.n_elem;
  if( m == -INFINITY ) return(m);
  for( int i = 0; i < n; i ++)
  {
    expsum += exp(avec(i) - m);
  }
  return(log(expsum) + m );
}

/* returns log(exp(a) + exp(b)) */ 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double logsumexp(const double a, const double b)
{
  double m = -999;
  if( a < b ) m = b;
  if( b <= a ) m = a;
  if( m == -INFINITY ) return(m);
  return(log(exp(a-m) + exp(b-m)) + m );
}

/* normalizes matrix so that rows sum to 1. */
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat normalize_A(arma::mat A)
{
  int n = A.n_rows; int p = A.n_cols;
  arma::mat A2(n,p);
  for( int i = 0; i < n ; i ++)
  {
    double sumi = arma::sum(A.row(i));
    A2.row(i) = A.row(i)/sumi;
  }
  return(A2);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat normalize_A_reg(arma::mat A)
{
  int n = A.n_rows; int p = A.n_cols;
  arma::mat A2(n,p);
  for( int i = 0; i < n ; i ++)
  {
    double sumi = 1e-14;
    sumi += arma::sum(A.row(i));
    A2.row(i) = A.row(i)/sumi;
  }
  return(A2);
}


/* normalizes vector to sum to 1. */
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat normalize_vec(arma::vec v)
{
  int n = v.n_elem; 
  arma::vec v2(n);
  v2 = v/arma::sum(v);
  return(v2);
}

/* normalizes vector to sum to 1. */
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat normalize_logvec(arma::vec lv)
{
  int n = lv.n_elem; 
  arma::vec lv2(n);
  lv2 = lv - logsumexp_vec(lv);
  return(exp( lv2 ));
}

/* returns counts for each allocation. */
 // [[Rcpp::depends(RcppArmadillo)]]
 // [[Rcpp::export]]
 arma::vec nvec_count(arma::vec svec, int k)
 {
   arma::vec n(k); n.fill(0); int ntot = 0;
   int nt = svec.n_elem;
   for( int i=0; i<nt; i ++)
   {
     for( int j = 0 ; j < k; j ++)
     {
       if( svec(i) == j+1 ) { n(j) += 1; ntot += 1; }
     }
   }
   if( ntot < nt ){ Rcout << "Warning: sum(nvec) < n" << std::endl;}
   return(n);
 }

/* gives list of clusters */
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec ordered_remove(arma::vec list, int index)
{
  // index is location (using 1:n indexing) of element to remove.
  int t = list.n_elem;
  arma::vec new_list(t-1); new_list.fill(0);
  
  if( t < (index) ) { Rcout<< "error: index > t"<< std::endl; return(new_list); }
  for( int i = 0; i < t; i ++)
  {
    if( i < (index-1) ) { new_list(i) = list(i); }
    if( i > (index-1) ) { new_list(i-1) = list(i); }
  }
  return(new_list);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec ordered_insert_next(arma::vec list)
{
  // inserts cluster in lowest unused spot
  int t = list.n_elem;
  arma::vec new_list(t+1); new_list.fill(0);
  int i = t;
  while(  i > 0 && i < list(i-1) )
  {
    //if( i < list(i-1) )
    //{
    new_list(i) = list(i-1);  
    //}
    i += -1 ;
  }
  new_list(i) = i+1;
  
  if( i > 0 )
  {
    for( int h=0; h < i; h ++ )
    {
      new_list(h) = list(h);
    }
  }
   
  return(new_list);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int ordered_next(arma::vec list)
{
  // index of lowest empty cluster label
  int t = list.n_elem;
  arma::vec new_list(t+1); new_list.fill(0);
  int i = t;
  while(  i > 0 && i < list(i-1) )
  {
    //if( i < list(i-1) )
    //{
    new_list(i) = list(i-1);  
    //}
    i += -1 ;
  }
  //new_list(i) = i+1;
  
  
  return(i+1);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
bool in_vec_int(double n, arma::vec v)
{
  bool invec = false;
  int vlen = v.n_elem;
  for( int j=(vlen-1); j >= 0; j--)
  {
    if( v(j)==n ) { invec = true; } 
  }
  return( invec );
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int which_vec_int(int n, arma::vec v)
{
  int out = -99; int nmatch = 0;
  int vlen = v.n_elem;
  for( int j=(vlen-1); j >= 0; j--)
  {
   if( v(j)==n ) { nmatch += 1; out = j+1; } 
  }
  if( nmatch > 1 ) {Rcout<< "warning: which returns FIRST of " << nmatch  <<" matches." << std::endl;}
  return( out);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec int_seq(const int k)
{
  arma::vec out(k);
  for( int j=1; j < (k+1); j ++)
  {
    out(j-1) = j;
  }
  return(out);
}


