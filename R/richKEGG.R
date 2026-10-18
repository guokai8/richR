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
#' @param minGSSize minimum annotated term size for reporting; positive-overlap terms are adjusted before this reporting filter.
#' @param maxGSSize maximal size of each geneSet for analyzing
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param builtin use KEGG built-in annotation or not (set FALSE if you want to use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @param universe optional character vector defining the ORA background
#' @param annotation_provenance optional annotation provenance captured from an Annot object
#' @author Kai Guo
richKEGG_internal<-function(x,kodata,pvalue=0.05,padj=NULL,ontology="KEGG",
                            organism=NULL,keytype="SYMBOL",minSize=2,maxSize=500,
                            minGSSize = 10, maxGSSize = 500,
                            keepRich=TRUE, filename=NULL,padj.method="BH",builtin=TRUE,sep=",",universe=NULL,
                            annotation_provenance=NULL){
    .validateParams(pvalue=pvalue, padj=padj, minSize=minSize, maxSize=maxSize,
                    minGSSize=minGSSize, maxGSSize=maxGSSize, func_name="richKEGG")
    normalization <- .normalize_annotation(kodata, func_name = "richKEGG")
  ## KEGG-specific: term name lookup + pathway level annotation
  if (ontology == "KEGGM") {
    all_ko <- .get_kgm.data()
  } else {
    all_ko <- .get_kg_dat(builtin = builtin)
  }
  term_names <- setNames(all_ko[, 1], rownames(all_ko))
  path_env <- new.env(parent = emptyenv())
  data("path", envir = path_env)
    .run_ora(x, annot = normalization, term_names = term_names, extra_cols = path_env$path,
           pvalue = pvalue, padj = padj, padj.method = padj.method,
           minSize = minSize, maxSize = maxSize,
           minGSSize = minGSSize, maxGSSize = maxGSSize,
           keepRich = keepRich, organism = organism, ontology = ontology,
           keytype = keytype, filename = filename, sep = sep, universe = universe,
           annotation_provenance = annotation_provenance,
           metadata_ontology = ontology)
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
#' @param minGSSize minimum annotated term size for reporting; positive-overlap terms are adjusted before this reporting filter.
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
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setMethod("richKEGG", signature(kodata = "data.frame"),definition = function(x,kodata,pvalue=0.05,padj=NULL,organism=NULL,ontology="KEGG",
                                                                             keytype=NULL,minSize=2,maxSize=500,minGSSize = 10, maxGSSize = 500,
                                                                             keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",",universe=NULL) {
  richKEGG_internal(x,kodata=kodata,pvalue=pvalue,padj=padj,
                    organism=organism,ontology=ontology,keytype=keytype,minSize=minSize,
                    maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                    padj.method=padj.method,builtin=builtin,sep=sep,universe=universe)
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
#' @param minGSSize minimum annotated term size for reporting; positive-overlap terms are adjusted before this reporting filter.
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
#' @param universe optional character vector intersected with annotated genes to define the ORA background
setMethod("richKEGG", signature(kodata = "Annot"),definition = function(x,kodata,pvalue=0.05,padj=NULL,organism=NULL,ontology="KEGG",
                                                                        keytype=NULL,minSize=2,maxSize=500,minGSSize = 10, maxGSSize = 500,
                                                                        keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",",universe=NULL) {
  provenance <- .annotation_provenance(kodata)
  richKEGG_internal(x,kodata@annot,pvalue=pvalue,padj=padj,
                  organism=kodata@species,ontology=kodata@anntype,keytype=kodata@keytype,minSize=minSize,
                  maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,keepRich=keepRich,filename=filename,
                  padj.method=padj.method,builtin=builtin,sep=sep,universe=universe,
                  annotation_provenance=provenance)
})



