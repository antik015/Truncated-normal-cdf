#include <RcppArmadillo.h>
#include <Rmath.h>
#include <RcppArmadilloExtensions/sample.h>
#include <omp.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

// [[Rcpp::export]]
double cpp_tmvnorm_prob(vec l, vec u, vec m, mat S, int N){
  int p = S.n_cols;
  mat L(p, p, fill::zeros);
  L = chol(S, "lower");
  vec pr(N, fill::zeros);
  for(int i =0; i<N; ++i){
    vec v(p, fill::zeros);
    double p_i = 1;
    for(int j = 0; j < p; ++j){
      if(j==0){
        double a_j = (l(j) - m(j))/L(j,j);
        double b_j = (u(j) - m(j))/L(j,j);
        v(j) = R::qnorm((R::pnorm(b_j, 0, 1, 1, 0) - R::pnorm(a_j, 0, 1, 1, 0))*randu() + R::pnorm(a_j, 0, 1, 1, 0), 0, 1, 1, 0);
        p_i = p_i*(R::pnorm(b_j, 0, 1, 1, 0) - R::pnorm(a_j, 0, 1, 1, 0));
      } else {
        vec x = (L.row(j)).t();
        double a_j = (l(j) - m(j) - sum(x.subvec(0, j-1)%v.subvec(0, j-1)))/L(j,j);
        double b_j = (u(j) - m(j) - sum(x.subvec(0, j-1)%v.subvec(0, j-1)))/L(j,j);
        v(j) = R::qnorm((R::pnorm(b_j, 0, 1, 1, 0) - R::pnorm(a_j, 0, 1, 1, 0))*randu() + R::pnorm(a_j, 0, 1, 1, 0), 0, 1, 1, 0);
        p_i = p_i*(R::pnorm(b_j, 0, 1, 1, 0) - R::pnorm(a_j, 0, 1, 1, 0));
      }
    }
    pr(i) = p_i;
  }
  return mean(pr);
}