# include <RcppArmadillo.h>
# include <omp.h>
# include <math.h>
# include <Rmath.h>
# include <R.h>
# include <iostream>
# include "helpers.h"
using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat melt_probmatrix(const arma::mat probM)
{
  // inserts cluster in lowest unused spot
  const int n = probM.n_rows;
  const int k = probM.n_cols;

  arma::mat new_probM(n*k,3); new_probM.fill(0);
  int i=0;
    for( int ii=0; ii<n; ii++)
	{
	for( int j=0; j<k; j++)
	{
	new_probM(i,0) = probM(ii,j);
        new_probM(i,1) = j+1;
	new_probM(i,2) = ii+1; 
	i+=1;
	}
	}
  
  return(new_probM);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat sort_probmatrix(const arma::mat melt_probM, const int sort_col)
{
   int n = melt_probM.n_rows;
   int k = melt_probM.n_cols;
   arma::mat sorted(n,k);
   arma::uvec sorted_ind(n);
   sorted_ind = sort_index(melt_probM.col(sort_col));

   int q = n-1;
   for( int i=0; i<n; i++)
	{
	sorted.row(q) = melt_probM.row(sorted_ind(i));
    q += -1;
	}

  return(sorted);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
bool check_if_done(const arma::vec current_counts, const arma::vec max_counts)
{
  bool done=true;
  int k = max_counts.n_elem;
  for( int j=0; j<k; j++)
  {
    if(current_counts(j)<max_counts(j))
    {
      done = false;
    }
  }
  return(done);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
bool in_matrix(const int value, const arma::mat M)
{
  bool in_M=false;
  int m = M.n_rows;
  for( int j=0; j<m; j++)
  {
    if(in_vec_int(value, M.row(j).t()))
    {
      in_M = true;
    }
  }
  return(in_M);
}

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat find_anchored(const arma::mat P_sorted, const arma::vec m_vec, const int k)
{
  int max_m = max(m_vec);
  arma::mat  A_out(k, max_m); A_out.fill(0);
  
  arma::vec counter(k);  counter.fill(0);
  bool all_full=false;
  
  int i=0;
  const int P_nrows =P_sorted.n_rows;
  while(!all_full &(i < P_nrows) )
  {
    
    if( (counter(P_sorted(i,1)-1)<(m_vec(P_sorted(i,1)-1))) & !(in_matrix(P_sorted(i,2),A_out)))
    {
    A_out(P_sorted(i,1)-1,counter(P_sorted(i,1)-1)) = P_sorted(i,2);
    //A_out(P_sorted(i,1)-1,1) = P_sorted(i,2);
      
    counter(P_sorted(i,1)-1) +=1;
    }
     
    i +=1;
    all_full = check_if_done(counter, m_vec);
  }
  
  return(A_out);
}





