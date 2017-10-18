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
  int N = D.size();
  
  std::vector<arma::mat> A(N);
  std::vector<arma::mat> Ainv(N);
  
  for (int i = 0; i < N; i++) {
    arma::mat Di = as<arma::mat>(D[i]);
 
    arma::mat matA = Di.t() * Di;
    A[i] = matA;
    Ainv[i] = arma::inv(matA);
  }
  
  arma::mat S(ncols, ncols, arma::fill::zeros);
  
  for ( arma::uword i = 0; i < N; i++ ) {
    for ( arma::uword j = i + 1; j < N; j++) {
      arma::mat Ai = A[i];
      arma::mat A_j = Ainv[j];
      arma::mat Aj = A[j];
      arma::mat A_i = Ainv[i];
      S = S + Ai * A_j + Aj * A_i;
    }
  }
  
  
  S = S / (N * (N - 1));
  
  return S;
}

