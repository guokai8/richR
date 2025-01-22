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
  go2gene<-sf(godata)
  all_go<-.get_go_dat(ont=ontology)
  go2gene<-go2gene[names(go2gene)%in%rownames(all_go)];
  gene2go<-reverseList(go2gene)
  if(is.data.frame(x)){
    input=rownames(x)
  }else{
    input=as.vector(x)
  }
  fgene2go<-gene2go[input];
  fgo2gene<-reverseList(fgene2go)
  k=name_table(fgo2gene);
  n=sum(!is.na(names(fgene2go)))
  IGO<-names(fgo2gene);
  N <- length(unique(unlist(go2gene)));
  M <- name_table(go2gene[IGO])
  rhs<-hyper_bench_vector(k,M,N,n)
  lhs<-p.adjust(rhs,method=padj.method)
  rhs_an<-all_go[names(rhs),]
  rhs_gene<-unlist(lapply(fgo2gene, function(x)paste(unique(x),sep="",collapse = sep)))
  Annotated=M[names(rhs)]
  Significant=k[names(rhs)]
  RichFactor <- Significant / Annotated
  FoldEnrichment <- RichFactor * N / n
  Pvalue=as.vector(rhs)
  GeneID=rhs_gene[names(rhs)]
  # mu and sigma are the mean and standard deviation of the hypergeometric distribution
  ## https://en.wikipedia.org/wiki/Hypergeometric_distribution
  mu <- M * n / N
  sigma <- mu * (N - n) * (N - M) / N / (N-1)
  zscore <- (k - mu)/sqrt(sigma)

  resultFis<-data.frame("Annot"=names(rhs),"Term"=rhs_an,"Annotated"=Annotated,
                        "Significant"=Significant,"RichFactor" = RichFactor,"FoldEnrichment"= FoldEnrichment,
                        "zscore"=zscore,"Pvalue"=Pvalue,"Padj"=lhs,
                        "GeneID"=GeneID)
  resultFis<-resultFis[order(resultFis$Pvalue),]
  resultFis<-subset(resultFis, Significant<=maxSize)
  ## remove gene Set with too much gene annotated
  resultFis<-subset(resultFis, Annotated<=maxGSSize)
  if(keepRich==FALSE){
    resultFis<-subset(resultFis, Significant>=minSize)
    resultFis<-subset(resultFis, Annotated>=minGSSize)
  }else{
    resultFis<-subset(resultFis, Significant>=minSize|RichFactor==1|Annotated >=minGSSize)
  }
  if(is.null(padj)){
    resultFis<-resultFis[resultFis$Pvalue<pvalue,]
    padj=numeric()
  }else{
    resultFis<-resultFis[resultFis$Padj<padj,]
  }
  rownames(resultFis)<-resultFis$Annot
  if(!is.null(filename)){
    write.table(resultFis,file=paste(filename,".txt",sep=""),sep="\t",quote=F,row.names=F)
  }
  if(is.data.frame(x)){
    detail<-getdetail(resultFis,x,sep=sep)
  }else{
    if(length(as.vector(resultFis$GeneID)>=1)){
      gene<-strsplit(as.vector(resultFis$GeneID),split=sep)
      names(gene)<-resultFis$Annot
      gened<-data.frame("TERM"=rep(names(gene),times=unlist(lapply(gene,length))),
                        "Annot"=rep(resultFis$Term,times=unlist(lapply(gene,length))),
                        "GeneID"=unlist(gene),row.names=NULL,
                        "Pvalue"=rep(resultFis$Pvalue,times=unlist(lapply(gene,length))),
                        "Padj"=rep(resultFis$Padj,times=unlist(lapply(gene,length)))
      )
    }else{
      gene = x
      names(gene)<-resultFis$Annot
      gened<-data.frame("TERM"="",
                        "Annot"="",
                        "GeneID"=x,row.names=NULL,
                        "Pvalue"=1,
                        "Padj"=1)
    }
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
              sep=sep
  )
  return(result);
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
