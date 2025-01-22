#include <Rcpp.h>
using namespace Rcpp;
using namespace std;
template <int RTYPE>
SEXP fast_factor_template( const Vector<RTYPE>& x ) {
  Vector<RTYPE> levs = sort_unique(x);
  IntegerVector out = match(x, levs);
  out.attr("levels") = as<CharacterVector>(levs);
  out.attr("class") = "factor";
  return out;
}

//[[Rcpp::export]]
SEXP fast_factor( SEXP x ) { //modified from RcppCore
  switch( TYPEOF(x) ) {
  case INTSXP: return fast_factor_template<INTSXP>(x);
  case REALSXP: return fast_factor_template<REALSXP>(x);
  case STRSXP: return fast_factor_template<STRSXP>(x);
  }
  return R_NilValue;
}

//[[Rcpp::export]]
List sf(DataFrame &x){
  StringVector Gene=x(0);
  StringVector Fa=x(1);
  StringVector tmp=unique(Fa);
  List res(Fa.size());
  Function sp("split");
  res=sp(Gene,fast_factor(Fa));
  return(res);
}
