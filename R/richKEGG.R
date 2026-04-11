#' KEGG Pathway Enrichment analysis function
#' @importFrom dplyr filter
#' @importFrom rlang sym
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata GO annotation data
#' @param ontology KEGG
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimal size of genes annotated by ontology term for testing.
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param builtin use KEGG built-in annotation or not (set FALSE if you want to use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @author Kai Guo
richKEGG_internal<-function(x,kodata,pvalue=0.05,padj=NULL,ontology="KEGG",
                            organism=NULL,keytype="SYMBOL",minSize=2,maxSize=500,
                            minGSSize = 10, maxGSSize = 500,
                            keepRich=TRUE, filename=NULL,padj.method="BH",builtin=TRUE,sep=","){
    .validateParams(pvalue=pvalue, padj=padj, minSize=minSize, maxSize=maxSize,
                    minGSSize=minGSSize, maxGSSize=maxGSSize, func_name="richKEGG")
    .validateGeneInput(x, annotation=kodata, func_name="richKEGG")
  ## KEGG-specific: term name lookup + pathway level annotation
  if (ontology == "KEGGM") {
    all_ko <- .get_kgm.data()
  } else {
    all_ko <- .get_kg_dat(builtin = builtin)
  }
  term_names <- setNames(all_ko[, 1], rownames(all_ko))
  path_env <- new.env(parent = emptyenv())
  data("path", envir = path_env)
  .run_ora(x, annot = kodata, term_names = term_names, extra_cols = path_env$path,
           pvalue = pvalue, padj = padj, padj.method = padj.method,
           minSize = minSize, maxSize = maxSize,
           minGSSize = minGSSize, maxGSSize = maxGSSize,
           keepRich = keepRich, organism = organism, ontology = ontology,
           keytype = keytype, filename = filename, sep = sep)
}
#' KEGG Pathway Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata KEGG annotation data
#' @param ontology KEGG
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimal size of genes annotated by ontology term for testing.
#' @param maxGSSize maximal size of each geneset for analyzing
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param builtin use KEGG built-in annotation or not (set FALSE if you want to use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   hsako<-as.data.frame(hsako)
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#' }
#' @export
#' @author Kai Guo
setMethod("richKEGG", signature(kodata = "data.frame"),definition = function(x,kodata,pvalue=0.05,padj=NULL,organism=NULL,ontology="KEGG",
                                                                             keytype=NULL,minSize=2,maxSize=500,minGSSize = 10, maxGSSize = 500,
                                                                             keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",") {
  richKEGG_internal(x,kodata=kodata,pvalue=pvalue,padj=padj,
                    organism=organism,ontology=ontology,keytype=keytype,minSize=minSize,
                    maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                    padj.method=padj.method,builtin=builtin,sep=sep)
})

#' KEGG Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param kodata Annot object with KEGG annotation
#' @param ontology KEGG
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimal size of genes annotated by ontology term for testing.
#' @param maxGSSize maximal size of each geneset for analyzing
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param builtin use KEGG built-in annotation or not (set FALSE if you want to use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#' }
#' @export
#' @author Kai Guo
setMethod("richKEGG", signature(kodata = "Annot"),definition = function(x,kodata,pvalue=0.05,padj=NULL,organism=NULL,ontology="KEGG",
                                                                        keytype=NULL,minSize=2,maxSize=500,minGSSize = 10, maxGSSize = 500,
                                                                        keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",") {
  richKEGG_internal(x,kodata@annot,pvalue=pvalue,padj=padj,
                  organism=kodata@species,ontology=kodata@anntype,keytype=kodata@keytype,minSize=minSize,
                  maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                  padj.method=padj.method,builtin=builtin,sep=sep)
})



