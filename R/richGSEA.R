#' Enrichment analysis for any type of annotation data
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param filename output filename
#' @param padj.method p value adjust method
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @importFrom fgsea fgsea
#' @export
#' @author Kai Guo
richGSEA_internal<-function(x,object,keytype="",pvalue=0.05,padj=NULL,minSize=15,maxSize=500,nperm=5000,filename=NULL,
                            padj.method="BH",organism=NULL,ontology=NULL,table=TRUE){
  x<-sort(x)
  if("Annot"%in%colnames(object)){
    object[,2]<-object$Annot
  }
  annod<-sf(object);
  res<-fgsea(pathways=annod,stats=x,minSize=minSize,maxSize=maxSize,nperm=nperm)
  if(isTRUE(table)){
    res$leadingEdge<-unlist(lapply(res$leadingEdge, function(x)paste(gsub(" ","",x),collapse = ",",sep="")))
  }
  if(is.null(padj)){
    res<-res[res$pval<pvalue,]
    padj=numeric()
  }else{
    res<-res[res$padj<padj,]
  }
  res<-res[order(res$pval),]
  if(is.null(organism)){
    organism=character()
  }
  if(is.null(ontology)){
    ontology=character()
  }
  result<-new("GSEAResult",
              result=res,
              pvalueCutoff   = pvalue,
              pAdjustMethod  = padj.method,
              padjCutoff   = padj,
              genenumber    = length(x),
              organism       = organism,
              gene           = names(x),
              ontology = ontology,
              keytype        = keytype
  )
  return(result)
}

#' GSEA Enrichment analysis function
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param filename output filename
#' @param padj.method p value adjust method
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @export
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   hsako<-as.data.frame(hsako)
#'   name=sample(unique(hsako$GeneID),1000)
#'   gene<-rnorm(1000)
#'   names(gene)<-name
#'   res<-richKEGG(gene,kodata = hsako)
#' }
#' @author Kai Guo
setMethod("richGSEA", signature(object = "data.frame"),definition = function(x,object,keytype="",pvalue=0.05,padj=NULL,minSize=15,ontology=ontology,
                                                                             maxSize=500,nperm=5000,filename=NULL,padj.method="BH",organism=NULL,table=TRUE) {
  richGSEA_internal(x,object,keytype=keytype,pvalue=pvalue,padj=padj,minSize=minSize,ontology=ontology,
                    maxSize=maxSize,nperm=nperm,filename=filename,padj.method=padj.method,organism=organism,table=table)
})
#' GSEA Enrichment analysis function
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param filename output filename
#' @param padj.method p value adjust method
#' @param nperm Number of permutations to do. Minimial possible nominal p-value is about 1/nperm
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   name=sample(unique(hsako$GeneID),1000)
#'   gene<-rnorm(1000)
#'   names(gene)<-name
#'   res<-richGSEA(gene,kodata = hsako)
#' }
#' @export
#' @author Kai Guo
setMethod("richGSEA", signature(object = "Annot"),definition = function(x,object,keytype="",pvalue=0.05,padj=NULL,minSize=15,ontology=ontology,
                                                                            maxSize=500,nperm=5000,filename=NULL,padj.method="BH",organism=NULL,table=TRUE) {
  richGSEA_internal(x,object@annot,keytype=object@keytype,pvalue=pvalue,padj=padj,minSize=minSize,ontology=object@anntype,
                    maxSize=maxSize,nperm=nperm,filename=filename,padj.method=padj.method,organism=object@species,table=table)
})


#' @name ggGSEA
#' @title plot the gsea result
#' @param x a vector include all log2FC with gene name
#' @param object Annot object
#' @param term the significant term
#' @param info Term with annotation details
#' @param gseaRes GSEAResult object
#' @importFrom fgsea plotEnrichment plotGseaTable
#' @examples
#' \dontrun{
#' set.seed(123)
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   name=sample(unique(hsako$GeneID),1000)
#'   gene<-rnorm(1000)
#'   names(gene)<-name
#'   res<-richGSEA(gene,kodata = hsako)
#'  # ggGSEA(gene,term=res$)
#' }
#' @export
#' @author Kai Guo
ggGSEA<-function(x,term,object,gseaRes=gseaRes){
  x<-sort(x)
  annot<-object@annot
  gseaRes<-gseaRes@result
  if(!is.null(annot$Annot)){
    annot[,2]<-annot$Annot
  }
  annod <- sf(annot)
  if(length(term)>1&!is.null(gseaRes)){
    plotGseaTable(annod[term],stats=x,gseaRes,gseaParam=0.5)
  }else{
    plotEnrichment(annod[[term]],stats=x)
  }
}

