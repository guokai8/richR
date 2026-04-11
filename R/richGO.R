#' GO Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param godata GO annotation data
#' @param ontology BP,MF or CC
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param minGSSize minimal size of genes annotated by ontology term for testing.
#' @param maxGSSize maximal size of each geneset for analyzing
#' @param keepRich keep terms with rich factor value equal 1 or not (default: FALSE)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @author Kai Guo
richGO_internal<-function(x,godata,ontology="BP",pvalue=0.05,padj=NULL,
                 organism=NULL,keytype="SYMBOL",minSize=2,maxSize=500,
                 minGSSize = 10, maxGSSize = 500,
                 keepRich=FALSE, filename=NULL,padj.method="BH",sep=","){
  .validateParams(pvalue=pvalue, padj=padj, minSize=minSize, maxSize=maxSize,
                  minGSSize=minGSSize, maxGSSize=maxGSSize, func_name="richGO")
  .validateGeneInput(x, annotation=godata, func_name="richGO")
  ## GO-specific: filter annotation to terms in this ontology
  if ("ONTOLOGYALL" %in% colnames(godata)) {
    godata_filt <- godata[which(godata$ONTOLOGYALL == ontology), ]
  } else {
    all_go <- .get_go_dat(ont = ontology)
    godata_filt <- godata[godata[, 2] %in% rownames(all_go), ]
  }
  if (nrow(godata_filt) == 0) {
    stop("richGO: no GO terms found for ontology '", ontology,
         "'. Check that your annotation contains this ontology.", call. = FALSE)
  }
  ## term_names: named vector  GO_ID -> readable term
  if ("Annot" %in% colnames(godata_filt) && "GOALL" %in% colnames(godata_filt)) {
    ## Build from existing annotation columns (avoids requiring GO.db)
    tm <- unique(godata_filt[, c("GOALL", "Annot")])
    term_names <- setNames(tm$Annot, tm$GOALL)
    annot_2col <- godata_filt[, c("GeneID", "GOALL")]
  } else {
    all_go <- .get_go_dat(ont = ontology)
    term_names <- setNames(all_go[, 1], rownames(all_go))
    annot_2col <- godata_filt[, 1:2]
  }
  .run_ora(x, annot = annot_2col, term_names = term_names,
           pvalue = pvalue, padj = padj, padj.method = padj.method,
           minSize = minSize, maxSize = maxSize,
           minGSSize = minGSSize, maxGSSize = maxGSSize,
           keepRich = keepRich, organism = organism, ontology = ontology,
           keytype = keytype, filename = filename, sep = sep)
}
#' GO Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param godata GO annotation data
#' @param ontology BP,MF or CC
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
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
#'   hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#'   hsago <- as.data.frame(hsago)
#'   gene <- sample(unique(hsago$GeneID),1000)
#'   res<-richGO(gene,godata = hsago,ontology ="BP")
#' }
#' @export
#' @author Kai Guo
setMethod("richGO", signature(godata = "data.frame"),definition = function(x,godata,ontology="BP",pvalue=0.05,padj=NULL,
                                                                           organism=NULL,keytype=NULL,minSize=2,maxSize=500,
                                                                           minGSSize = 10, maxGSSize = 500,
                                                                           keepRich=TRUE, filename=NULL,padj.method="BH",sep=",") {
 # godata<-as(godata,"Annot")
  richGO_internal(x,godata,ontology=ontology,pvalue=pvalue,padj=padj,
                  organism=organism,keytype=keytype,minSize=minSize,maxSize=maxSize,
                  minGSSize = minGSSize, maxGSSize = maxGSSize,
                  keepRich=keepRich,filename=filename,padj.method=padj.method,sep=sep)
})

#' GO Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param godata GO annotation data
#' @param ontology BP,MF or CC
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
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
#'   hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#'   gene=sample(unique(hsago$GeneID),1000)
#'   res<-richGO(gene,godata = hsago,ontology ="BP")
#' }
#' @export
#' @author Kai Guo
setMethod("richGO", signature(godata = "Annot"),definition = function(x,godata,ontology="BP",pvalue=0.05,padj=NULL,minSize=2,maxSize=500,
                                                                      minGSSize = 10, maxGSSize = 500,
                                                                      keepRich=TRUE,filename=NULL,padj.method="BH",sep=",") {
  richGO_internal(x,godata@annot,ontology=ontology,pvalue=pvalue,padj=padj,
                  organism=godata@species,keytype=godata@keytype,minSize=minSize,maxSize=maxSize,minGSSize = minGSSize, maxGSSize = maxGSSize,
                  keepRich=keepRich,filename=filename,padj.method=padj.method,sep=sep)
})
