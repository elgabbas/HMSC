#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//' @export
// [[Rcpp::export]]
arma::mat Solve1(const arma::mat & x) {
    return arma::inv(x);
}

//' @export
// [[Rcpp::export]]
arma::mat Solve2(const arma::mat& A, const arma::mat& B) {
    return arma::solve(A, B);
}

//' @export
// [[Rcpp::export]]
arma::vec Solve2vect(const arma::mat& A, const arma::vec& B) {
    arma::vec x = arma::solve(A, B);
    // Ensure the result is a vector
    return x;
}
