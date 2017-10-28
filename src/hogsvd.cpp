#include <RcppArmadillo.h>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate normalised S matrix
//' @param D list of data matrices
//' @param ncols number of columns of input matrices
//' @param nthreads number of omp threads, 0 for max
//' @export
// [[Rcpp::export]]
arma::mat calcNormS(const List& D, int& ncols, int nthreads = 0, bool verbose = true) {
  arma::uword N = D.size();
  std::vector<arma::mat> A(N);
  std::vector<arma::mat> Ainv(N);
  arma::mat S(ncols, ncols, arma::fill::zeros);
  
#ifdef _OPENMP
  if (nthreads == 0) {
    nthreads = omp_get_max_threads(); 
  }
  
  omp_set_num_threads(nthreads);
  
  if (verbose)
    Rcpp::Rcout << " using up to " << nthreads << " threads...";
#endif
  
#pragma omp parallel for
  for ( arma::uword i = 0; i < N; i++) {
    arma::mat Di = as<arma::mat>(D[i]); 
    arma::mat matA = Di.t() * Di;
    A[i] = matA;
    Ainv[i] = arma::inv(matA);
  }
  
  if (verbose)
    Rcpp::Rcout << " A, A inverse computation complete... ";
  
  // Force thread sync
  {
    #pragma omp barrier  
  }
  
  // Check for user interupt
  Rcpp::checkUserInterrupt();
  
  /* Old style double loop that is bad for omp
  
  for ( arma::uword i = 0; i < N; i++ ) {
    for ( arma::uword j = i + 1; j < N; j++) {
      arma::mat tmp = A[i] * Ainv[j] + A[j] * Ainv[i]; 
      #pragma omp critical
      {
        S = S + tmp;
      }
    }
  }

  */
  
  std::vector< std::pair<int,int> > pairs;
  for (int i = 0; i < N; i++) 
    for (int j = i + 1; j < N; j++)
      pairs.push_back(std::pair<int,int>(i,j));
  
#pragma omp parallel for
  for(int k = 0; k < N * (N+1) / 2; k++) {
    std::pair<int,int> p = pairs[k];
    int i = p.first;
    int j = p.second;

    arma::mat tmp = A[i] * Ainv[j] + A[j] * Ainv[i]; 
#pragma omp critical
    {
      S = S + tmp;
    }    
  }

  
  // Force thread sync
  {
    #pragma omp barrier  
  }
  
  // Check for user interupt
  Rcpp::checkUserInterrupt();
  
  S = S / (N * (N - 1));
  
  return S;
}

