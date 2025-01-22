#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List reverseList(const List& lhs) {
  Function sf("sf");  // 从 R 全局环境中获取名为 "sf" 的函数
  
  // 1) 拿到 lhs 的名字向量
  CharacterVector nm = lhs.names();
  int n = lhs.size(); // 列表的长度
  
  // 2) 统计所有子向量元素总数，为后面预分配做准备
  R_xlen_t total_length = 0;
  std::vector<R_xlen_t> lens(n);
  for(int i = 0; i < n; i++){
    SEXP elem = lhs[i];
    R_xlen_t len = Rf_xlength(elem); 
    lens[i] = len;
    total_length += len;
  }
  
  CharacterVector lhs_n(total_length);
  CharacterVector unl(total_length);
  
  // 4) 逐个子向量收集数据
  R_xlen_t idx = 0;
  for(int i = 0; i < n; i++){
    R_xlen_t len = lens[i];
    if(len == 0) continue;
    
    // 取出第 i 个子向量 (假设是字符串向量，可按需改成别的类型)
    CharacterVector sub = lhs[i];
    
    for(R_xlen_t j = 0; j < len; j++){
      lhs_n[idx] = nm[i];  // 重复 names(lhs)[i]
      unl[idx]   = sub[j]; // 拿到子向量中的元素
      idx++;
    }
  }
  
  DataFrame df = DataFrame::create(
    _["V1"] = lhs_n,
    _["V2"] = unl,
    _["stringsAsFactors"] = false
  );
  
  // 6) 调用您现有的 sf() 函数
  //    如果 sf() 是 C++ 导出的函数，也可以直接调用；若它在 R 里，可用 Function("sf")
  List res = sf(df);
  
  // 返回结果
  return res;
}
