#include <RcppArmadillo.h>
#include <vector>

#include <omp.h>
// [[Rcpp::plugins(openmp)]]

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate normalised S matrix
//' @param D list of data matrices
//' @param ncols number of columns of input matrices
//' @param nthreads number of omp threads, 0 for max
//' @export
// [[Rcpp::export]]
arma::mat calcNormS(const List& D, int& ncols, int nthreads = 1) {
  arma::uword N = D.size();
  std::vector<arma::mat> A(N);
  std::vector<arma::mat> Ainv(N);
  arma::mat S(ncols, ncols, arma::fill::zeros);
  
  if (nthreads == 0) {
   nthreads = omp_get_max_threads(); 
  }

  omp_set_num_threads(nthreads);
  #pragma omp parallel for
  for ( arma::uword i = 0; i < N; i++) {
    arma::mat Di = as<arma::mat>(D[i]); 
    arma::mat matA = Di.t() * Di;
    A[i] = matA;
    Ainv[i] = arma::inv(matA);
  }

  #pragma omp parallel for
  for ( arma::uword i = 0; i < N; i++ ) {
    for ( arma::uword j = i + 1; j < N; j++) {
      arma::mat tmp = A[i] * Ainv[j] + A[j] * Ainv[i]; 
      #pragma omp critical
      {
        S = S + tmp;
      }
    }
  }
    
  S = S / (N * (N - 1));
  
  return S;
}

