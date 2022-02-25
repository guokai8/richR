##' richplot
##' @description plot the sigificant terms and shared genes with network format
##' @importFrom GGally ggnet2
##' @importFrom igraph graph_from_data_frame
##' @importFrom igraph simplify
##' @importFrom ggrepel geom_text_repel
##' @importFrom igraph V
##' @importFrom igraph V<-
##' @importFrom igraph degree
##' @importFrom ggplot2 geom_text
##' @param object richResult or dataframe
##' @param top number of terms to show (default: 50)
##' @param pvalue cutoff p value for enrichment result
##' @param padj cutoff p adjust value for enrichment result
##' @param usePadj use adjust p value for the color or not
##' @param low color used for small value
##' @param high color used for large value
##' @param useTerm use terms for nodes (default: TRUE)
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param label.color label color
##' @param label.size label size
##' @param node.shape vector of shape and names of the vector should be the terms (default: 20)
##' @param layout layout method
##' @param savefig save figures or not
##' @param filename output figure name
##' @param width width for output figure
##' @param height height for output figure
##' @param node.alpha alpha-transparency scales
##' @param node.shape shape of the node
##' @param repel use ggrepel text function or not
##' @param segment.size segment size for ggrepel text
##' @param sep character string used to separate the genes when concatenating
ggrich_internal <- function(object,top=50, pvalue=0.05, padj=NULL,
                            usePadj =TRUE, useTerm=TRUE,low="orange",high="red",
                            writeCyt=FALSE, cytoscapeFile="network-file-for-cytoscape.txt",
                            label.color = "black", label.size = 2, node.shape=NULL,
                            layout = layout.fruchterman.reingold,savefig=FALSE,filename="network",
                            width=7,height=7,node.alpha=0.7,repel=TRUE,segment.size=0.2,sep=","){
  if(!is.null(padj)){
    object<-object[object$Padj<padj,]
  }else{
    object<-object[object$Pvalue<pvalue,]
  }
  if(nrow(object)>=top){
    object<-object[1:top,]
  }
  object$GeneID<-as.vector(object$GeneID)
  lhs <- strsplit(object$GeneID,sep)
  if(isTRUE(useTerm)){
    rhs <- data.frame(from=unlist(lhs),to=rep(object$Term,lapply(lhs,length)))
    rownames(object) <- object$Term
  }else{
    rhs <- data.frame(from=unlist(lhs),to=rep(object$Annot,lapply(lhs,length)))
    rownames(object) <- object$Annot
  }
  if(isTRUE(usePadj)){
    rhs$value <- -log10(object[rhs$to,"Padj"])

  }else{
    rhs$value <- -log10(object[rhs$to,"Pvalue"])
  }
  if (writeCyt == TRUE) {
    write.table(rhs, file = cytoscapeFile, sep = "\t",
                row.names = F, quote = F)
  }
  pvalue1 <- rhs$value
  names(pvalue1) <- rhs$from
  pvalue1 <- unlist(lapply(split(pvalue1,names(pvalue1)),mean))
  pvalue2 <- rhs$value
  names(pvalue2)<-rhs$to
  pvalue2 <- unlist(lapply(split(pvalue2,names(pvalue2)),mean))
  pvalue<-c(pvalue1,pvalue2)
  g <- graph_from_data_frame(rhs, directed = F)
  pvalue = pvalue[V(g)$name]
  cols <- .color_scale(high, low)
  V(g)$color <- cols[sapply(pvalue, .getIdx, min(pvalue), max(pvalue))]
  if(!is.null(node.shape)){
    shape=rep(20,length(V(g)$name))
    names(shape)<-V(g)$name
    shape[names(node.shape)]<-node.shape
    node.shape <- shape
  }else{
    node.shape=rep(20,length(V(g)$name))
    names(node.shape)<-V(g)$name
  }
  p <- ggnet2(g, node.size = degree(g), node.color = V(g)$color,
              node.shape=node.shape,legend.position = "none",
              node.alpha=node.alpha,edge.size=0.05)
  if(isTRUE(repel)){
    p <- p + geom_text_repel(label = V(g)$name,size=label.size, color=label.color,
                             segment.size=segment.size)
  }else{
    p <- p +geom_text(label = V(g)$name,size=label.size, color=label.color)
  }
  if(savefig==TRUE){
    ggsave(p,file=paste(filename,"pdf",sep="."),width=width,height = height)
  }
  p
}
##' richplot for Enrichment results
##' @param object richResult or dataframe
##' @param top number of terms to show (default: 50)
##' @param pvalue cutoff p value for enrichment result
##' @param padj cutoff p adjust value for enrichment result
##' @param usePadj use adjust p value for the color or not
##' @param low color used for small value
##' @param high color used for large value
##' @param useTerm use terms for nodes (default: TRUE)
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param label.color label color
##' @param label.size label size
##' @param node.shape vector of shape and names of the vector should be the terms (default: 20)
##' @param layout layout method
##' @param savefig save figures or not
##' @param filename output figure name
##' @param width width for output figure
##' @param height height for output figure
##' @param node.alpha alpha-transparency scales
##' @param node.shape shape of the node
##' @param repel use ggrepel text function or not
##' @param segment.size segment size for ggrepel text
##' @param sep character string used to separate the genes when concatenating
##' @examples
##' \dontrun{
##'   hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
##'   gene=sample(unique(hsago$GeneID),1000)
##'   res<-richKEGG(gene,kodata = hsako)
##'   ggnetplot(res)
##' }
##' @export
##' @author Kai Guo
setMethod("ggnetplot", signature(object = "richResult"),definition = function(object,top=50, pvalue=0.05, padj=NULL,
                                                                           usePadj =TRUE, useTerm=TRUE,low="orange",high="red",
                                                                           writeCyt=FALSE, cytoscapeFile="network-file-for-cytoscape.txt",
                                                                           label.color = "black", label.size = 2, node.shape=NULL,
                                                                           layout = layout.fruchterman.reingold,savefig=FALSE,filename="network",
                                                                           width=7,height=7,node.alpha=0.7,repel=TRUE,segment.size=0.2,sep=",") {
  ggrich_internal(object@result,top=top,pvalue=pvalue,padj=padj,
                 usePadj=usePadj,useTerm=useTerm,low=low,high=high,
                 writeCyt=writeCyt, cytoscapeFile=cytoscapeFile,
                 label.color = label.color, label.size = label.size, node.shape=node.shape,
                 layout = layout,savefig=savefig,filename=filename,
                 width=width,height=height,node.alpha=node.alpha,repel=repel,segment.size=segment.size,sep=object@sep)
})
##' richplot for Enrichment result
##' @rdname ggnetplot
##' @param object richResult or dataframe
##' @param top number of terms to show (default: 50)
##' @param pvalue cutoff p value for enrichment result
##' @param padj cutoff p adjust value for enrichment result
##' @param usePadj use adjust p value for the color or not
##' @param low color used for small value
##' @param high color used for large value
##' @param useTerm use terms for nodes (default: TRUE)
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param label.color label color
##' @param label.size label size
##' @param node.shape vector of shape and names of the vector should be the terms (default: 20)
##' @param layout layout method
##' @param savefig save figures or not
##' @param filename output figure name
##' @param width width for output figure
##' @param height height for output figure
##' @param node.alpha alpha-transparency scales
##' @param node.shape shape of the node
##' @param repel use ggrepel text function or not
##' @param segment.size segment size for ggrepel text
##' @param sep character string used to separate the genes when concatenating
##' @examples
##' \dontrun{
##'   hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
##'   gene=sample(unique(hsago$GeneID),1000)
##'   res<-richKEGG(gene,kodata = hsako)
##'   ggnetplot(result(res))
##' }
##' @export
##' @author Kai Guo
setMethod("ggnetplot", signature(object = "data.frame"),definition = function(object,top=50, pvalue=0.05, padj=NULL,
                                                                           usePadj =TRUE, useTerm=TRUE,low="orange",high="red",
                                                                           writeCyt=FALSE, cytoscapeFile="network-file-for-cytoscape.txt",
                                                                           label.color = "black", label.size = 2, node.shape=NULL,
                                                                           layout = layout.fruchterman.reingold,savefig=FALSE,filename="network",
                                                                           width=7,height=7,node.alpha=0.7,repel=TRUE,segment.size=0.2,sep=",") {
  ggrich_internal(object,top=top,pvalue=pvalue,padj=padj,
                  usePadj=usePadj,useTerm=useTerm,low=low,high=high,
                  writeCyt=writeCyt, cytoscapeFile=cytoscapeFile,
                  label.color = label.color, label.size = label.size, node.shape=node.shape,
                  layout = layout,savefig=savefig,filename=filename,
                  width=width,height=height,node.alpha=node.alpha,repel=repel,segment.size=segment.size,sep=sep)
})

##' richplot for Enrichment result
##' @rdname ggnetplot
##' @param top number of terms to show (default: 50)
##' @param pvalue cutoff p value for enrichment result
##' @param padj cutoff p adjust value for enrichment result
##' @param usePadj use adjust p value for the color or not
##' @param low color used for small value
##' @param high color used for large value
##' @param useTerm use terms for nodes (default: TRUE)
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param label.color label color
##' @param label.size label size
##' @param node.shape vector of shape and names of the vector should be the terms (default: 20)
##' @param layout layout method
##' @param savefig save figures or not
##' @param filename output figure name
##' @param width width for output figure
##' @param height height for output figure
##' @param node.alpha alpha-transparency scales
##' @param node.shape shape of the node
##' @param repel use ggrepel text function or not
##' @param segment.size segment size for ggrepel text
##' @examples
##' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsako<-as.data.frame(hsako)
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richKEGG(gene,object = hsako)
#' ggnetplot(res)
##' }
setMethod("ggnetplot", signature(object = "GSEAResult"),definition = function(object,top=50, pvalue=0.05, padj=NULL,
                                                                              usePadj =TRUE, useTerm=TRUE,low="orange",high="red",
                                                                              writeCyt=FALSE, cytoscapeFile="network-file-for-cytoscape.txt",
                                                                              label.color = "black", label.size = 2, node.shape=NULL,
                                                                              layout = layout.fruchterman.reingold,savefig=FALSE,filename="network",
                                                                              width=7,height=7,node.alpha=0.7,repel=TRUE,segment.size=0.2,sep=",") {
  object<-object@result
  sep <- object@sep
  if(is.list(object$leadingEdge)){
    object$leadingEdge<-unlist(lapply(object$leadingEdge, function(x)paste(x,collapse = sep,sep="")))
  }
  object$Annot<-object$pathway
  object<-object[,c(9,1:3,8)]
  colnames(object)<-c("Annot","Term","Pvalue","Padj","GeneID")
  ggrich_internal(object,top=top,pvalue=pvalue,padj=padj,
                  usePadj=usePadj,useTerm=useTerm,low=low,high=high,
                  writeCyt=writeCyt, cytoscapeFile=cytoscapeFile,
                  label.color = label.color, label.size = label.size, node.shape=node.shape,
                  layout = layout,savefig=savefig,filename=filename,
                  width=width,height=height,node.alpha=node.alpha,repel=repel,segment.size=segment.size,sep=sep)
})




