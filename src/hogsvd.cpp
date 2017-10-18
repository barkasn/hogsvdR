#include <RcppArmadillo.h>
#include <vector>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate normalised S matrix
//' @param D list of data matrices
//' @param ncols number of columns of input matrices
//' @export
// [[Rcpp::export]]
arma::mat calcNormS(const List& D, int ncols) {
  arma::uword N = D.size();
  
  std::vector<arma::mat> A(N);
  std::vector<arma::mat> Ainv(N);

  #pragma omp parallel for
  for ( arma::uword i = 0; i < N; i++) {
    arma::mat Di = as<arma::mat>(D[i]);
 
    arma::mat matA = Di.t() * Di;
    A[i] = matA;
    Ainv[i] = arma::inv(matA);
  }
  
  arma::mat S(ncols, ncols, arma::fill::zeros);

  #pragma omp parallel for
  for ( arma::uword i = 0; i < N; i++ ) {
    for ( arma::uword j = i + 1; j < N; j++) {
      arma::mat tmp = A[i] * Ainv[j] + A[j] * Ainv[i]; 
      #pragma omp critical
      S = S + tmp;
    }
  }
  
  
  S = S / (N * (N - 1));
  
  return S;
}

