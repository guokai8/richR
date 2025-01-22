#include <Rcpp.h>
using namespace Rcpp;

inline double hyp_test(double k, double M, double NM, double n) {
  return R::phyper(k, M, NM, n, false, false);
}

// [[Rcpp::export]]
NumericVector hyper_bench_vector(NumericVector xin,
                                  NumericVector yin,
                                  double N,
                                  double n) {
  int xsize = xin.size();
  NumericVector res(xsize);
  double* px = &xin[0];
  double* py = &yin[0];
  double* pr = &res[0];

  for (int i = 0; i < xsize; i++) {
    double x_val = px[i];
    double y_val = py[i];
    pr[i] = hyp_test(x_val - 1.0, y_val, N - y_val, n);
  }

  res.attr("names") = xin.names();
  return res;
}

