#include <Rcpp.h>
using namespace Rcpp;
//[[Rcpp::export]]
StringVector uniq(StringVector& xa){
  StringVector rhs=unique(xa);
  return(rhs);
}
NumericVector uniq(NumericVector& xa){
  NumericVector rhs=unique(xa);
  return(rhs);
}
