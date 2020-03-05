#' Enrichment analysis for any type of annotation data
#' @importFrom dplyr filter_
#' @importFrom magrittr %>%
#' @param x vector contains gene names or dataframe with DEGs information
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
enrich_internal<-function(x,annot,pvalue=0.05,padj=NULL,organism=NULL,ontology="",minSize=1,maxSize=500,
                          keepRich=TRUE,keytype="",filename=NULL,padj.method="BH",sep = ","){
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
  resultFis<-resultFis%>%filter_(~Significant<=maxSize)
  if(keepRich==FALSE){
    resultFis<-resultFis%>%filter_(~Significant>=minSize)
  }else{
    resultFis<-resultFis%>%filter_(~Significant>=minSize|(~Annotated/~Significant)==1)
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
#' KEGG Pathway Enrichment analysis function
#' @param x vector contains gene names or dataframe with DEGs information
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
setMethod("enrich", signature(annot = "data.frame"),definition = function(x,annot,pvalue=0.05,padj=NULL,organism=NULL,ontology="",
                                                                             keytype="",filename=NULL,minSize=2,maxSize=500,
                                                                             keepRich=TRUE,padj.method="BH",sep=",") {
  enrich_internal(x,annot=annot,ontology=ontology,pvalue=pvalue,padj=padj,
                    organism=organism,keytype=keytype,minSize=minSize,maxSize=maxSize,keepRich=keepRich,
                    filename=filename,padj.method=padj.method,sep=sep)
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
#' @export
#' @author Kai Guo
setMethod("enrich", signature(annot = "Annot"),definition = function(x,annot,pvalue=0.05,padj=NULL,organism=NULL,ontology="",
                                                                        keytype="",filename=NULL,minSize=2,maxSize=500,
                                                                        keepRich=TRUE,padj.method="BH",builtin=TRUE,sep=",") {
  enrich_internal(x=x,annot=annot@annot,ontology=annot@anntype,pvalue=pvalue,padj=padj,
                    organism=annot@species,keytype=annot@keytype,minSize=minSize,maxSize=maxSize,keepRich=keepRich,
                    filename=filename,padj.method=padj.method,sep=sep)
})
