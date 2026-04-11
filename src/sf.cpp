#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <string>
using namespace Rcpp;

template <int RTYPE>
SEXP fast_factor_template( const Vector<RTYPE>& x ) {
  Vector<RTYPE> levs = sort_unique(x);
  IntegerVector out = match(x, levs);
  out.attr("levels") = as<CharacterVector>(levs);
  out.attr("class") = "factor";
  return out;
}

//[[Rcpp::export]]
SEXP fast_factor( SEXP x ) {
  switch( TYPEOF(x) ) {
  case INTSXP: return fast_factor_template<INTSXP>(x);
  case REALSXP: return fast_factor_template<REALSXP>(x);
  case STRSXP: return fast_factor_template<STRSXP>(x);
  }
  return R_NilValue;
}

// Pure C++ split: group col0 values by col1 keys using a hash map.
// Returns a named List where names = unique col1 values,
// each element = CharacterVector of col0 values in that group.
//[[Rcpp::export]]
List sf(DataFrame &x){
  SEXP col0 = VECTOR_ELT(x, 0);
  SEXP col1 = VECTOR_ELT(x, 1);
  R_xlen_t n = Rf_xlength(col0);

  // Pass 1: collect indices per group key (col1 value).
  // Use CHARSXP pointer as key for speed (R interns strings).
  std::unordered_map<SEXP, std::vector<R_xlen_t>> groups;
  groups.reserve(n / 4);
  // Maintain insertion order via a separate vector of unique keys.
  std::vector<SEXP> order;
  order.reserve(n / 4);

  for (R_xlen_t i = 0; i < n; i++) {
    SEXP key = STRING_ELT(col1, i);
    auto it = groups.find(key);
    if (it == groups.end()) {
      order.push_back(key);
      groups[key].push_back(i);
    } else {
      it->second.push_back(i);
    }
  }

  // Pass 2: build output List.
  int ngroups = order.size();
  List res(ngroups);
  CharacterVector res_names(ngroups);

  for (int g = 0; g < ngroups; g++) {
    SEXP key = order[g];
    const std::vector<R_xlen_t>& idx = groups[key];
    int m = idx.size();
    CharacterVector vals(m);
    for (int j = 0; j < m; j++) {
      SET_STRING_ELT(vals, j, STRING_ELT(col0, idx[j]));
    }
    res[g] = vals;
    SET_STRING_ELT(res_names, g, key);
  }

  res.attr("names") = res_names;
  return res;
}
