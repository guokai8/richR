#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// Count unique elements per list entry.
// Uses CHARSXP pointer comparison — R interns strings so equal strings
// share the same pointer, avoiding std::string copies entirely.
// [[Rcpp::export]]
IntegerVector name_table(List lh) {
  int n = lh.size();
  IntegerVector res(n);

  SEXP nm = Rf_getAttrib(lh, R_NamesSymbol);

  for (int i = 0; i < n; i++) {
    SEXP elem = lh[i];
    if (TYPEOF(elem) != STRSXP) {
      res[i] = NA_INTEGER;
      continue;
    }

    R_xlen_t len = Rf_xlength(elem);
    // SEXP pointers are unique per interned string — pointer hash is O(1).
    std::unordered_set<SEXP> uset;
    uset.reserve(len);

    for (R_xlen_t j = 0; j < len; j++) {
      uset.insert(STRING_ELT(elem, j));
    }
    res[i] = (int)uset.size();
  }

  Rf_setAttrib(res, R_NamesSymbol, nm);
  return res;
}
