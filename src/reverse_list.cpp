#include <Rcpp.h>
#include <unordered_map>
#include <vector>
using namespace Rcpp;

// Reverse a named list: input {A: [x,y], B: [y,z]} -> {x: [A], y: [A,B], z: [B]}
// Uses CHARSXP pointers as hash keys (R interns strings, so pointer == identity).
// [[Rcpp::export]]
List reverseList(const List& lhs) {

  SEXP names_sexp = Rf_getAttrib(lhs, R_NamesSymbol);
  int n = lhs.size();

  // Pass 1: build reversed mapping  value -> vector of source names.
  std::unordered_map<SEXP, std::vector<SEXP>> rev;
  std::vector<SEXP> order;  // insertion order of unique values

  for (int i = 0; i < n; i++) {
    SEXP nm_i = STRING_ELT(names_sexp, i);
    SEXP elem = lhs[i];
    if (TYPEOF(elem) != STRSXP) continue;
    R_xlen_t len = Rf_xlength(elem);
    for (R_xlen_t j = 0; j < len; j++) {
      SEXP val = STRING_ELT(elem, j);
      auto it = rev.find(val);
      if (it == rev.end()) {
        order.push_back(val);
        rev[val].push_back(nm_i);
      } else {
        it->second.push_back(nm_i);
      }
    }
  }

  // Pass 2: build output List.
  int ngroups = order.size();
  List res(ngroups);
  CharacterVector res_names(ngroups);

  for (int g = 0; g < ngroups; g++) {
    SEXP key = order[g];
    const std::vector<SEXP>& src = rev[key];
    int m = src.size();
    CharacterVector vals(m);
    for (int j = 0; j < m; j++) {
      SET_STRING_ELT(vals, j, src[j]);
    }
    res[g] = vals;
    SET_STRING_ELT(res_names, g, key);
  }

  res.attr("names") = res_names;
  return res;
}
