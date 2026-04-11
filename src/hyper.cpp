#include <Rcpp.h>
using namespace Rcpp;

// Vectorized hypergeometric p-values for ORA.
// xin = significant genes per term (k), yin = annotated genes per term (M),
// N = total background genes, n = total input genes found in annotation.
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
    double k = px[i];
    double M = py[i];
    // Skip phyper when no overlap or empty term (p = 1.0).
    if (k <= 0.0 || M <= 0.0) {
      pr[i] = 1.0;
    } else {
      pr[i] = R::phyper(k - 1.0, M, N - M, n, false, false);
    }
  }

  res.attr("names") = xin.names();
  return res;
}

