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
#' @param keepRich keep terms with rich factor value equal 1 or not (default: TRUE)
#' @param bulitin use KEGG bulit in KEGG annotation or not(set FALSE if you want use newest KEGG data)
#' @param filename output filename
#' @param padj.method pvalue adjust method(default:"BH")
#' @param sep character string used to separate the genes when concatenating
#' @author Kai Guo
richKEGG_internal<-function(x,kodata,pvalue=0.05,padj=NULL,ontology="KEGG",
                            organism=NULL,keytype="SYMBOL",minSize=2,maxSize=500,
                            keepRich=TRUE, filename=NULL,padj.method="BH",builtin=TRUE,sep=","){

    ko2gene<-sf(kodata)
    ko2gene_num<-name_table(ko2gene)
    gene2ko<-sf(kodata[,c(2,1)])
    if(is.data.frame(x)){
      input=rownames(x)
    }else{
      input=as.vector(x)
    }
    fgene2ko=gene2ko[input]
    fko2gene=reverseList(fgene2ko)
    k=name_table(fko2gene)
    n=length(unique(unlist(fko2gene)))
    M=ko2gene_num[names(k)]
    N=length(unique(kodata[,1]))
    rhs<-hyper_bench_vector(k,M,N,n)
    lhs<-p.adjust(rhs,method=padj.method)
    if(ontology=="KEGGM"){
      all_ko<-.get_kgm.data()
    }else{
      all_ko<-.get_kg_dat(builtin=builtin)
    }
    rhs_an<-all_ko[names(rhs),]
    rhs_gene<-unlist(lapply(fko2gene, function(x)paste(unique(x),sep="",collapse = sep)))
    resultFis<-data.frame("Annot"=names(rhs),"Term"=rhs_an,"Annotated"=M[names(rhs)],
                          "Significant"=k[names(rhs)],"Pvalue"=as.vector(rhs),"Padj"=lhs,
                          "GeneID"=rhs_gene[names(rhs)])
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
      resultFis<-subset(resultFis, Significant>=minSize|(Significant/Annotated)==1)
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
  if(is.null(keytype)){
    keytype=character()
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
#' KEGG Pathway Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param ontology KEGG
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
                                                                             keytype=NULL,minSize=2,maxSize=500,
                                                                             keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",") {
  richKEGG_internal(x,kodata=kodata,pvalue=pvalue,padj=padj,
                    organism=organism,ontology=ontology,keytype=keytype,minSize=minSize,
                    maxSize=maxSize,keepRich=keepRich,filename=filename,
                    padj.method=padj.method,builtin=builtin,sep=sep)
})

#' KEGG Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
#' @param ontology KEGG
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
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#' }
#' @export
#' @author Kai Guo
setMethod("richKEGG", signature(kodata = "Annot"),definition = function(x,kodata,pvalue=0.05,padj=NULL,organism=NULL,ontology="KEGG",
                                                                        keytype=NULL,minSize=2,maxSize=500,
                                                                        keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,sep=",") {
  richKEGG_internal(x,kodata@annot,pvalue=pvalue,padj=padj,
                  organism=kodata@species,ontology=kodata@anntype,keytype=kodata@keytype,minSize=minSize,
                  maxSize=maxSize,keepRich=keepRich,filename=filename,
                  padj.method=padj.method,builtin=builtin,sep=sep)
})



