#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector name_table(List lh) {
  // 长度
  int n = lh.size();
  IntegerVector res(n);

  // 提取 List 的 names
  SEXP nm = Rf_getAttrib(lh, R_NamesSymbol);

  for (int i = 0; i < n; i++) {
    // 第 i 个元素
    SEXP elem = lh[i];
    // 确认它是字符串向量
    if (TYPEOF(elem) != STRSXP) {
      // 如果不是字符串向量，可以根据需求决定抛错或标记 NA
      res[i] = NA_INTEGER;
      continue;
    }

    // 元素长度
    R_xlen_t len = Rf_xlength(elem);
    std::unordered_set<std::string> uset;
    uset.reserve(len);

    for (R_xlen_t j = 0; j < len; j++) {
      // 取得 R 的 CHARSXP
      SEXP rstr = STRING_ELT(elem, j);
      if (rstr == NA_STRING) {
        uset.insert("NA");
      } else {
        // 拷贝到 std::string
        uset.insert(std::string(CHAR(rstr)));
      }
    }
    res[i] = uset.size();
  }

  Rf_setAttrib(res, R_NamesSymbol, nm);
  return res;
}
