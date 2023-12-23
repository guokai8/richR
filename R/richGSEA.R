#' Enrichment analysis for any type of annotation data
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off method
#' @param padj.method p value adjust method
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
#' @importFrom fgsea fgseaMultilevel
#' @export
#' @author Kai Guo
richGSEA_internal<-function(x,object,keytype="",pvalue=0.05,padj=NULL,minSize=15,maxSize=500,
                            padj.method="BH",organism=NULL,ontology=NULL,table=TRUE,sep=","){
  x<-sort(x)
  if("Annot"%in%colnames(object)){
    object[,2]<-object$Annot
  }
  annod<-sf(object);
  res<-fgseaMultilevel(pathways=annod,stats=x,minSize=minSize,maxSize=maxSize)
  if(isTRUE(table)){
    res$leadingEdge<-unlist(lapply(res$leadingEdge, function(x)paste(gsub(" ","",x),collapse = sep,sep="")))
  }
  if(is.null(padj)){
    res<-res[res$pval<pvalue,]
    padj=numeric()
  }else{
    res<-res[res$padj<padj,]
  }
  res<-res[order(res$pval),]
  res <- res[!is.na(res$pathway),]
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
              input = x,
              ontology = ontology,
              keytype        = keytype,
              sep = sep
  )
  return(result)
}

#' GSEA Enrichment analysis function
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off method
#' @param padj.method p value adjust method
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
#' @export
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsako<-as.data.frame(hsako)
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' }
#' @author Kai Guo
setMethod("richGSEA", signature(object = "data.frame"),definition = function(x,object,keytype="",pvalue=0.05,padj=NULL,minSize=15,ontology=ontology,
                                                                             maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=",") {
  richGSEA_internal(x,object,keytype=keytype,pvalue=pvalue,padj=padj,minSize=minSize,ontology=ontology,
                    maxSize=maxSize,padj.method=padj.method,organism=organism,table=table,sep=sep)
})
#' GSEA Enrichment analysis function
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off method
#' @param padj.method p value adjust method
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' }
#' @export
#' @author Kai Guo
setMethod("richGSEA", signature(object = "Annot"),definition = function(x,object,keytype="",pvalue=0.05,padj=NULL,minSize=15,ontology=ontology,
                                                                            maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=",") {
  richGSEA_internal(x,object@annot,keytype=object@keytype,pvalue=pvalue,padj=padj,minSize=minSize,ontology=object@anntype,
                    maxSize=maxSize,padj.method=padj.method,organism=object@species,table=table,sep=sep)
})


#' @name ggGSEA
#' @title plot the gsea result
#' @param x a vector include all log2FC with gene name
#' @param object Annot object
#' @param term the significant term
#' @param info Term with annotation details
#' @param gseaRes GSEAResult object
#' @importFrom fgsea plotEnrichment plotGseaTable
#' @importFrom ggplot2 ggtitle
#' @importFrom cowplot plot_grid
#' @examples
#' \dontrun{
#' set.seed(123)
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' ggGSEA(gene,term = res$pathway,object = hsako,gseaRes = res,default = F)
#' }
#' @export
#' @author Kai Guo
ggGSEA<-function(x,term,object,gseaRes=gseaRes,top=10, default = TRUE){
  x<-sort(x)
  annot<-object@annot
  gseaRes<-gseaRes@result
  if(nrow(gseaRes)<top){
    top <- nrow(gseaRes)
  }
  if(!is.null(annot$Annot)){
    annot[,2] <- annot$Annot
  }
  annod <- sf(annot)
  if(length(term)>1&!is.null(gseaRes)){
    if(isTRUE(default)){
       plotGseaTable(annod[term],stats=x,gseaRes,gseaParam=0.5)
    }else{
      res<-lapply(gseaRes$pathway, function(y)plotEnrichment(annod[[y]],stats=x)+ggtitle(y))
      plot_grid(plotlist = res[1:top],ncol=5)
    }
  }else{
    plotEnrichment(annod[[term]],stats=x)
  }
}

#' plot multiple significant pathways
#' @importFrom ggplot2 ggplot geom_hline aes geom_point geom_segment
#' @importFrom ggplot2 geom_line theme_bw element_blank theme scale_color_manual
#' @param x a vector include all log2FC with gene name
#' @param object Annot object
#' @param gseaRes GSEAResult object
#' @param mycol a vector indicate the colors used for the figure
#' @param top number of terms you want to display,
#' @param pvalue cutoff value of pvalue (if padj set as NULL)
#' @param padj cutoff value of p adjust value
#' @examples
#' \dontrun{
#' set.seed(123)
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' plotGSEA(hsako,res,object = hsako)
#' }
#' @export
#' @author Kai Guo
plotGSEA<-function(x, object,gseaRes,mycol=NULL,top=10,pvalue=0.05,padj=NULL,
                   gseaParam = 1, ticksSize = 0.2){
  if(!is.null(padj)){
    cutoff <- padj
    sigpathway<-gseaRes$pathway[gseaRes$padj<cutoff]
  }else{
    cutoff<-pvalue
    sigpathway<-gseaRes$pathway[gseaRes$pval<cutoff]
  }
  if(top>length(sigpathway)){
    top<-length(sigpathway)
  }
  if(is.null(mycol)){
    mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9",
                          "#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D",
                          "#7CC767")
  }
  sigpathway<-sigpathway[1:top]
  fc <- x
  res<-lapply(sigpathway,function(x).calGSEA(object,x,fc,gseaParam=gseaParam,ticksSize=ticksSize))
  toPlot<-do.call(rbind,lapply(res,'[[','toPlot'))
  pathway<-do.call(rbind,lapply(res,'[[','pathway'))
  tops<-unlist(lapply(res,'[[','tops'))
  bottoms<-unlist(lapply(res,'[[','bottoms'))
  diff <- (max(tops) - min(bottoms))/8
  p<-ggplot(toPlot, aes(x = x, y = y,color=Group)) + geom_point(size = 0.1) +
    geom_hline(yintercept = max(tops), colour = "red",
               linetype = "dashed") + geom_hline(yintercept = min(bottoms),
                                                 colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "black") + geom_line() + theme_bw()
  p <- p+geom_segment(data =pathway, mapping = aes(x = x,
                                              y = -diff/4, xend = x,
                                              yend = diff/4,color=Group),
                 size = ticksSize) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(values=mycol)+labs(x = "rank", y = "Enrichment score")
  return(p)
}


