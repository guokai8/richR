##' calculate kappa cluster
##' @param x richResult object
##' @param gene (Optional).a vector of gene list
##' @param useTerm to use the term or not (TRUE/FALSE)
##' @param cutoff kappa score threshold for significant dispersion results
##' @param overlap
##' @param minSize minimal number of terms in the cluster
##' @param escore kappa enrichment score cutoff value (default: 3)
##' @author Kai Guo
.kappa_cluster_internal<-function(x,gene=NULL,useTerm=FALSE,cutoff=0.5,overlap=0.5,minSize=5,escore=3){
  if(isTRUE(useTerm)){
    rownames(x) <- x$Term
  }else{
    rownames(x) <- x$Annot
  }
  mat<-expand.grid(rownames(x),rownames(x),stringsAsFactors=F)
  mat <- mat[mat$Var1>mat$Var2,]
  mat$kappa <- apply(mat,1,function(y).kappa(x[y[1],"GeneID"],x[y[2],"GeneID"],gene))
  mat<-mat[mat$kappa > cutoff,]
  if(nrow(mat)==0){
    stop("kappa cutoff too high\n")
  }
  ml1<-split(mat$Var1,mat$Var2)
  ml2<-split(mat$Var2,mat$Var1)
  ###
  ml<-sapply(union(names(ml1),names(ml2)),function(x)c(ml1[[x]],ml2[[x]]))
  ###
  res<-.merge_term(ml,overlap)
  #### remove the class with smaller number
  res<-res[unlist(lapply(res, function(x)length(x)>=minSize))]
  ###find the smallest p value
  idx<-lapply(res,function(y)which.min(x[y,"Pvalue"]))
  ###use the term with smallest p value as name
  names(res)<-unlist(lapply(names(idx),function(x)res[[x]][idx[[x]]]))
  res<-sapply(unique(names(res)),function(x)unique(unlist(res[names(res)==x])))
  #####
  es <- lapply(res, function(y).calculate_Enrichment_Score(y,x))
  tl <- unlist(lapply(res,length))
  es <- unlist(es)/tl
  es <- es[es>escore]
  if(length(es)==0){
    stop("Enrichment score too high, No significant cluster\n")
  }
  dx <- x[names(es),]
  tl <- tl[names(es)]
  dy <-data.frame(AnnotationCluster=1:length(es),EnrichmentScore=es,filteredClusterSize=tl)
  rs<-cbind(dy,dx)
  rs<-rs[order(rs$EnrichmentScore,decreasing = T),]
  rs$Cluster<-unlist(lapply(res[rownames(rs)],function(x)paste(x,collapse=",",sep="")))[rownames(rs)]
  return(rs)
}
#' kappa cluster analysis
##' @param x richResult object or dataframe
##' @param gene (Optional).a vector of gene list
##' @param useTerm to use the term or not (TRUE/FALSE)
##' @param cutoff kappa score threshold for significant dispersion results
##' @param overlap cutoff value of the overlap between two Terms
##' @param minSize minimal number of terms in the cluster
##' @param escore kappa enrichment score cutoff value (default: 3)
#' @examples
#' \dontrun{
#'   hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#'   gene=sample(unique(hsago$GeneID),1000)
#'   res<-richGO(gene,godata = hsago,ontology ="BP")
#'   resc<-richCluster(res)
#' }
#' @export
#' @author Kai Guo
setMethod("richCluster", signature(x = "data.frame"),definition = function(x,gene=NULL,useTerm=FALSE,cutoff=0.5,overlap=0.5,minSize=5,escore=3) {
  if(is.null(gene)){
    gene=unique(unlist(strsplit(x$GeneID,",")))
  }
  .kappa_cluster_internal(x,gene=gene,useTerm=useTerm,cutoff=cutoff,overlap=overlap,minSize=minSize,escore=escore)
})
#' kappa cluster analysis
##' @param x richResult object or dataframe
##' @param gene (Optional).a vector of gene list
##' @param useTerm to use the term or not (TRUE/FALSE)
##' @param cutoff kappa score threshold for significant dispersion results
##' @param overlap cutoff value of the overlap between two Terms
##' @param minSize minimal number of terms in the cluster
##' @param escore kappa enrichment score cutoff value (default: 3)
#' @examples
#' \dontrun{
#'   hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#'   gene=sample(unique(hsago$GeneID),1000)
#'   res<-richGO(gene,godata = hsago,ontology ="BP")
#'   resc<-richCluster(res)
#' }
#' @export
#' @author Kai Guo
setMethod("richCluster", signature(x = "richResult"),definition = function(x,gene=NULL,useTerm=FALSE,cutoff=0.5,overlap=0.5,minSize=5,escore=3) {
     if(is.null(gene)){
       gene=x@gene
     }
    .kappa_cluster_internal(x@result,gene=gene,useTerm=useTerm,cutoff=cutoff,overlap=overlap,minSize=minSize,escore=escore)
})
