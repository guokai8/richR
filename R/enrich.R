#' Enrichment analysis for any type of annotation data
#' @param x vector contains gene names or dataframe with DEGs information
#' @param annot ontology type
#' @param ontology ontology type
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @export
#' @author Kai Guo
enrich_internal<-function(x,annot,ontology= "",pvalue=0.05,padj=NULL,organism=NULL,minSize=1,maxSize=500,
                          keepRich=TRUE,keytype="",filename=NULL,padj.method="BH",sep = ","){
  ontology=""
  ao2gene<-sf(annot)
  ao2gene_num<-name_table(ao2gene)
  gene2ao<-sf(annot[,c(2,1)])
  if(is.data.frame(x)){
    input<-rownames(x)
  }else{
    input=as.vector(x)
  }
  fgene2ao=gene2ao[input]
  fao2gene=reverseList(fgene2ao)
  k=name_table(fao2gene)
  n=length(unique(unlist(fao2gene)))
  M=ao2gene_num[names(k)]
  N=length(unique(annot[,1]))
  rhs<-hyper_bench_vector(k,M,N,n)
  lhs<-p.adjust(rhs,method=padj.method)
  rhs_gene<-unlist(lapply(fao2gene, function(x)paste(unique(x),sep="",collapse = sep)))
  resultFis<-data.frame("Annot"=names(rhs),"Term"=names(rhs),"Annotated"=M[names(rhs)],
                        "Significant"=k[names(rhs)],"Pvalue"=as.vector(rhs),"Padj"=lhs,
                        "GeneID"=rhs_gene[names(rhs)])
  resultFis<-resultFis[order(resultFis$Pvalue),]
  if(is.null(padj)){
    resultFis<-resultFis[resultFis$Pvalue<pvalue,]
    padj=numeric()
  }else{
    resultFis<-resultFis[resultFis$Padj<padj,]
  }
  colnames(resultFis)[2]="Term"
  resultFis<-subset(resultFis,Significant<=maxSize)
  if(keepRich==FALSE){
    resultFis<-subset(resultFis,Significant>=minSize)
  }else{
    resultFis<-subset(resultFis,Significant>=minSize|(Significant/Annotated)==1)
  }
  rownames(resultFis)<-resultFis$Annot
  if(!is.null(filename)){
    write.table(resultFis,file=paste(filename,".txt",sep=""),sep="\t",quote=F,row.names=F)
  }
  if(is.data.frame(x)){
    detail<-getdetail(resultFis,x)
  }else{
    gene<-strsplit(as.vector(resultFis$GeneID),split=sep)
    names(gene)<-resultFis$Annot
    gened<-data.frame("TERM"=rep(names(gene),times=unlist(lapply(gene,length))),
                      "Annot"=rep(resultFis$Term,times=unlist(lapply(gene,length))),
                      "GeneID"=unlist(gene),row.names=NULL,
                      "Pvalue"=rep(resultFis$Pvalue,times=unlist(lapply(gene,length))),
                      "Padj"=rep(resultFis$Padj,times=unlist(lapply(gene,length)))
                      )
    gened$GeneID<-as.character(gened$GeneID)
    detail<-gened
  }
  if(is.null(organism)){
    organism=character()
  }
  result<-new("richResult",
              result=resultFis,
              detail=detail,
              pvalueCutoff   = pvalue,
              pAdjustMethod  = padj.method,
              padjCutoff   = padj,
              genenumber    = length(input),
              organism       = organism,
              ontology       = ontology,
              gene           = input,
              keytype        = keytype,
              sep = sep
  )
  return(result)
}
#' Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param annot ontology type
#' @param ontology ontology type
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @export
#' @author Kai Guo
setMethod("enrich", signature(annot = "data.frame"),definition = function(x,annot,ontology="",pvalue=0.05,padj=NULL,organism=NULL,
                                                                             keytype="",filename=NULL,minSize=2,maxSize=500,
                                                                             keepRich=TRUE,padj.method="BH",sep=",") {
  enrich_internal(x,annot=annot,ontology=ontology,pvalue=pvalue,padj=padj,
                    organism=organism,keytype=keytype,minSize=minSize,maxSize=maxSize,keepRich=keepRich,
                    filename=filename,padj.method=padj.method,sep=sep)
})

#' Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param annot ontology type
#' @param ontology ontology type
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param organism organism
#' @param keytype keytype for input genes
#' @param minSize minimal number of genes included in significant terms
#' @param maxSize maximum number of genes included in significant terms
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param bulitin use KEGG bulit in KEGG annotation or not(set FALSE if you want use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @export
#' @author Kai Guo
setMethod("enrich", signature(annot = "Annot"),definition = function(x,annot,pvalue=0.05,padj=NULL,organism=NULL,
                                                                        keytype="",filename=NULL,minSize=2,maxSize=500,
                                                                        keepRich=TRUE,padj.method="BH",builtin=TRUE,sep=",") {
  enrich_internal(x=x,annot=annot@annot,ontology=annot@anntype,pvalue=pvalue,padj=padj,
                    organism=annot@species,keytype=annot@keytype,minSize=minSize,maxSize=maxSize,keepRich=keepRich,
                    filename=filename,padj.method=padj.method,sep=sep)
})


#' Functional enrichment analysis with DAVID
#'
#' @param gene vector contains gene names
#' @param keytype key type export
#' @param anntype GOTERM_BP_FAT, KEGG_PATHWAY,GOTERM_CC_FAT
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param david.user richR@und.edu
#' @return Annot object
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="ENTREZID",anntype = "KEGG")
#'   hsako<-as.data.frame(hsako)
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richDAVID(gene,keytype="ENTREZID",species="human")
#' }
#' @export
#' @author Kai Guo
richDAVID <- function(gene,keytype="ENTREZ_GENE_ID",species="human",anntype="GOTERM_BP_FAT",
                      universe,pvalue=0.05,padj=NULL,minSize=2,maxSize=500,
                      keepRich=TRUE, filename=NULL,padj.method="BH",sep=",",
                      david.user="richR@und.edu"){
  pkg <- "RDAVIDWebService"
  if (!require(pkg,character.only=TRUE)){
    if(!require("BiocManager",character.only=TRUE)){
      install.packages("BiocManager")
    }else{
      BiocManager::install(pkg)
    }
  }else{
    suppressMessages(requireNamespace(pkg))
  }
  #require(pkg, character.only=TRUE)
  idtype=keytype
  if(keytype=="SYMBOL"){
    gene<-idconvert(gene,species=species,fkeytype = "SYMBOL",tkeytype = "ENTREZID")
    keytype="ENTREZ_GENE_ID"
  }else if(keytype=="ENTREZID"){
    keytype="ENTREZ_GENE_ID"
  }else if(keytype=="ENSEMBL"){
    keytype="ENSEMBL_GENE_ID"
  }else if(keytype=="REFSEQ"){
    keytype="REFSEQ_MRNA"
  }else{
    keytype=keytype
  }
  david <- DAVIDWebService$new(email=david.user,
                               url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  keytype <- match.arg(keytype, getIdTypes(david))
  rh <- addList(david, gene, idType=keytype,
                listName="richR",
                listType="Gene")
  if (rh$inDavid == 0) {
    stop("Please check 'keytype' parameter...")
  }
  if (!missing(universe)) {
    rh <- addList(david, universe, idType=keytype,
                  listName="universe",
                  listType="Background")
  }
  #####
  setAnnotationCategories(david, anntype)
  res <- getFunctionalAnnotationChart(david, threshold=1, count=minSize)
  if (length(res@.Data) == 0) {
    warning("No significant enrichment found...")
    return(NULL)
  }
  terms <- res$Term
  if (length(grep("~", terms[1])) == 0) {
    sep <- ":"
  } else {
    sep <- "~"
  }
  terml <- sapply(terms, function(y) strsplit(y, split=sep))
  term <- do.call("rbind", terml)
  Annot <- term[,1]
  Term <- term[,2]
  Annotated <- res$Pop.Hits
  Significant <- res$Count
  resultFis <- data.frame(Annot          = Annot,
                          Term = Term,
                          Annotated   = Annotated,
                          Significant     = Significant,
                          Pvalue      = res$PValue,
                          stringsAsFactors = FALSE)
  row.names(resultFis) <- Annot
  if (padj.method == "bonferroni") {
    resultFis$Padj <- res$Bonferroni
  } else if(padj.method == "FDR") {
    resultFis$Padj <- res$FDR
  }else{
    resultFis$Padj <- res$Benjamini
  }
  resultFis$GeneID <- gsub(' ','',res$Genes)
  resultFis<-resultFis[order(resultFis$Pvalue),]
  if(is.null(padj)){
    resultFis<-resultFis[resultFis$Pvalue<pvalue,]
    padj=numeric()
  }else{
    resultFis<-resultFis[resultFis$Padj<padj,]
  }
  resultFis<-subset(resultFis, Significant<=maxSize)
  if(keepRich==FALSE){
    resultFis<-subset(resultFis, Significant>=minSize)
  }else{
    resultFis<-subset(resultFis, Significant>=minSize|(Annotated/Significant)==1)
  }
  rownames(resultFis)<-resultFis$Annot
  if(idtype=="SYMBOL"){
    xx<-suppressMessages(lapply(strsplit(resultFis$GeneID,","),function(x)idconvert(x,species=species,fkeytype = "ENTREZID",tkeytype = "SYMBOL")))
    xx<-unlist(lapply(xx,function(x)paste(x,sep="",collapse = ",")))
    resultFis$GeneID <- xx
  }
  if(!is.null(filename)){
    write.table(resultFis,file=paste(filename,ontology,"res.txt",sep="_"),sep="\t",quote=F,row.names=F)
  }
  genes<-strsplit(as.vector(resultFis$GeneID),split=sep)
  names(genes)<-resultFis$Annot
  gened<-data.frame("TERM"=rep(names(genes),times=unlist(lapply(genes,length))),
                    "Annot"=rep(resultFis$Term,times=unlist(lapply(genes,length))),
                    "GeneID"=unlist(genes),row.names=NULL,
                    "Pvalue"=rep(resultFis$Pvalue,times=unlist(lapply(genes,length))),
                    "Padj"=rep(resultFis$Padj,times=unlist(lapply(genes,length)))
  )
  gened$GeneID<-as.character(gened$GeneID)
  detail<-gened
  species <- getSpecieNames(david)
  species <- gsub("\\(.*\\)", "", species)
  new("richResult",
      result=resultFis,
      detail=detail,
      pvalueCutoff   = pvalue,
      pAdjustMethod  = padj.method,
      padjCutoff   = padj,
      genenumber    = length(gene),
      organism       = species,
      ontology       = anntype,
      gene           = gene,
      keytype        = keytype,
      sep=",")
}

