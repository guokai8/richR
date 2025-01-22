#' Pathway Enrichment analysis for different level function
#' @importFrom dplyr filter left_join
#' @importFrom rlang sym
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata KEGG annotation data
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param ontology ontology type
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimal size of genes annotated by ontology term for testing.
#' @param maxGSSize maximal size of each geneset for analyzing
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @author Kai Guo
richLevel_internal<-function(x,kodata,level="Level2",pvalue =0.05, padj=NULL,
                             organism=NULL,keytype="SYMBOL",ontology="",minSize=2,maxSize=500,
                             minGSSize = 10, maxGSSize = 500,
                             keepRich=TRUE, filename=NULL,padj.method="BH",sep=","){
  data(path)
  annot<-left_join(kodata,path,by=c('PATH'="ko"))
  annot<-annot[,c("GeneID",level)]
  annot<-annot[!is.na(annot[,2]),]
  result<-enrich(x = x,object = annot,pvalue=pvalue,padj=padj,organism=organism,minSize=minSize,
              maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,keytype=keytype,ontology=ontology,filename=filename,
              padj.method=padj.method,sep = sep)
  return(result)
}
#' Pathway Enrichment analysis for different level function
#' @importFrom dplyr filter left_join
#' @importFrom rlang sym
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata KEGG annotation data
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param ontology ontology type
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimal size of genes annotated by ontology term for testing.
#' @param maxGSSize maximal size of each geneset for analyzing
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   hsako<-as.data.frame(hsako)
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richLevel(gene,kodata = hsako,level="Level2")
#' }
#' @export
#' @author Kai Guo
setMethod("richLevel", signature(kodata = "data.frame"),definition = function(x,kodata,level="Level2",pvalue =0.05, padj=NULL,
                                                                              organism=NULL,keytype="SYMBOL",ontology="KEGG",minSize=2,maxSize=500,
                                                                              minGSSize = 10, maxGSSize = 500,
                                                                              keepRich=TRUE, filename=NULL,padj.method="BH",sep=",") {
  richLevel_internal(x,kodata=kodata,level=level, pvalue=pvalue,padj=padj,
                    organism=organism,keytype=keytype,ontology=ontology,minSize=minSize,
                    maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                    padj.method=padj.method,sep=sep)
})

#' Pathway Enrichment analysis for different level function
#' @importFrom dplyr filter left_join
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata KEGG annotation data
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param ontology ontology type
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimal size of genes annotated by ontology term for testing.
#' @param maxGSSize maximal size of each geneset for analyzing
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richLevel(gene,kodata = hsako,level="Level2")
#' }
#' @export
#' @author Kai Guo
setMethod("richLevel", signature(kodata = "Annot"),definition = function(x,kodata,level="Level2",pvalue =0.05, padj=NULL,
                                                                         organism=NULL,keytype="SYMBOL",ontology="",minSize=2,maxSize=500,
                                                                         minGSSize = 10, maxGSSize = 500,
                                                                         keepRich=TRUE, filename=NULL,padj.method="BH",sep=",") {
  richLevel_internal(x,kodata@annot,level=level,pvalue=pvalue,padj=padj,
                    organism=kodata@species,keytype=kodata@keytype,ontology=kodata@anntype,minSize=minSize,
                    maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                    padj.method=padj.method,sep=sep)
})
