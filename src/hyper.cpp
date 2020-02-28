#include <Rcpp.h>
using namespace Rcpp;
double hyp_test( std::vector<double>& xa) {
  double k=xa[0];
  double M=xa[1];
  double NM=xa[2];
  double n=xa[3];
  return R::phyper(k,M,NM,n,FALSE,FALSE);
}
//[[Rcpp::export]]
NumericVector hyper_bench_vector(NumericVector& xin,NumericVector& yin,double N,double n){
  int xsize=xin.size();
  //double xt;
  NumericVector res;
  StringVector xnames=xin.names();
  for(int i=0;i<xsize;++i){
    String tname=xnames[i];
    std::vector<double>xres;
    xres.push_back(xin[tname]-1.0);
    xres.push_back(yin[tname]);
    xres.push_back(N-yin[tname]);
    xres.push_back(n);
    double pp=hyp_test(xres);
    res.push_back(pp);
  }
  res.attr("names")=xnames;
  return res;
}
