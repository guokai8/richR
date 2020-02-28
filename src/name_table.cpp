#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector name_table(List& lh) {
  int n=lh.length();
  StringVector names=lh.names();
  NumericVector res(n);
  for(int i=0;i<n;++i){
    StringVector tmp=lh[i];
    //tmp=unique(tmp);
    res[i]=unique(tmp).size();
  }
  res.attr("names")=names;
  return res;
}
