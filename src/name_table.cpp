#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Count unique elements per list entry.
// Strategy: copy CHARSXP pointers into a reusable buffer, sort, count
// distinct runs. Sorting raw pointers is cache-friendly and avoids all
// hash-table overhead. O(n log n) but with tiny constants.
// [[Rcpp::export]]
IntegerVector name_table(List lh) {
  int n = lh.size();
  IntegerVector res(n);

  SEXP nm = Rf_getAttrib(lh, R_NamesSymbol);

  // Reusable buffer — grows as needed, never shrinks.
  std::vector<SEXP> buf;

  for (int i = 0; i < n; i++) {
    SEXP elem = lh[i];
    if (TYPEOF(elem) != STRSXP) {
      res[i] = NA_INTEGER;
      continue;
    }

    R_xlen_t len = Rf_xlength(elem);

    if (len <= 1) {
      res[i] = (int)len;
      continue;
    }

    // Two elements: simple pointer comparison.
    if (len == 2) {
      res[i] = (STRING_ELT(elem, 0) == STRING_ELT(elem, 1)) ? 1 : 2;
      continue;
    }

    // General case: sort pointers, count distinct.
    buf.resize(len);
    for (R_xlen_t j = 0; j < len; j++) {
      buf[j] = STRING_ELT(elem, j);
    }
    std::sort(buf.begin(), buf.end());
    int nuniq = 1;
    for (R_xlen_t j = 1; j < len; j++) {
      nuniq += (buf[j] != buf[j - 1]);
    }
    res[i] = nuniq;
  }

  Rf_setAttrib(res, R_NamesSymbol, nm);
  return res;
}
