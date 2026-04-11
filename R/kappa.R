## ----------------------------------------------------------------
##  Kappa-related helper functions
## ----------------------------------------------------------------

#' Compute kappa statistic between two gene sets
#' @param x comma-separated gene string
#' @param y comma-separated gene string
#' @param geneall vector of all genes (background)
#' @return numeric kappa score
#' @keywords internal
.kappa <- function(x, y, geneall) {
  x <- unlist(strsplit(x, ","))
  y <- unlist(strsplit(y, ","))
  if (length(intersect(x, y)) == 0) {
    kab <- 0
  } else {
    tmp <- matrix(0, 2, 2)
    tmp[1, 1] <- length(intersect(x, y))
    tmp[2, 1] <- length(setdiff(x, y))
    tmp[1, 2] <- length(setdiff(y, x))
    tmp[2, 2] <- length(setdiff(geneall, union(x, y)))
    oab <- (tmp[1, 1] + tmp[2, 2]) / sum(tmp)
    aab <- ((tmp[1, 1] + tmp[2, 1]) * (tmp[1, 1] + tmp[1, 2]) +
              (tmp[1, 2] + tmp[2, 2]) * (tmp[2, 1] + tmp[2, 2])) / (sum(tmp) * sum(tmp))
    if (aab == 1) {
      kab <- 0
    } else {
      kab <- (oab - aab) / (1 - aab)
    }
  }
  return(kab)
}

#' Calculate enrichment score for a cluster
#' @param x character vector of term names
#' @param df data.frame with Pvalue column
#' @return numeric enrichment score
#' @keywords internal
.calculate_Enrichment_Score <- function(x, df) {
  pvalue <- df[x, "Pvalue"]
  esp <- ifelse(pvalue == 0, 16, -log10(pvalue))
  es <- sum(esp)
}

#' Merge overlapping term clusters
#' @param x named list of term vectors
#' @param overlap overlap threshold for merging
#' @return list of merged term clusters
#' @keywords internal
.merge_term <- function(x, overlap) {
  ml <- x
  res <- list()
  for (i in names(ml)) {
    lhs <- setdiff(names(ml), i)
    for (j in lhs) {
      ov <- intersect(ml[[i]], ml[[j]])
      un <- union(ml[[i]], ml[[j]])
      ovl <- length(ov) / length(un)
      if (ovl > overlap) {
        res[[i]] <- c(i, un)
        ml <- ml[setdiff(names(ml), j)]
      } else {
        res[[i]] <- c(i, ml[[i]])
      }
    }
  }
  return(res)
}

## ----------------------------------------------------------------
##  Kappa cluster main functions
## ----------------------------------------------------------------

##' calculate kappa cluster
##' @param x richResult object
##' @param gene (Optional).a vector of gene list
##' @param useTerm to use the term or not (TRUE/FALSE)
##' @param cutoff kappa score threshold for significant dispersion results
##' @param overlap overlap between clusters
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
