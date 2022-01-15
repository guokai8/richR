##' generate network based on Enrichment results
##' @rdname ggnetwork
##' @param object richResult,GSEAResult object or dataframe
##' @param gene vector contains gene names or dataframe with DEGs information
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param weightcut cutoff valule for edge
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param layout layout for the network (layout.fruchterman.reingold)
##' @param low color used for small value
##' @param high color used for large value
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param cytoscapeFormat Character string giving the output file format
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param label.font label font
##' @param label.color label color
##' @param label.size label size
##' @param filename figure output name
##' @param savefig save the figure or not
##' @param width figure width
##' @param height figure height
##' @param sep character string used to separate the genes when concatenating
##' @importFrom reshape2 melt
##' @importFrom igraph graph.data.frame
##' @importFrom ggrepel geom_text_repel
##' @importFrom igraph delete.edges
##' @importFrom ggplot2 theme
##' @importFrom GGally ggnet2
##' @importFrom ggplot2 ggsave
##' @importFrom intergraph asNetwork
##' @importFrom igraph E
##' @importFrom igraph E<-
##' @importFrom igraph V
##' @importFrom igraph V<-
##' @importFrom igraph write_graph
##' @importFrom visNetwork visIgraph
##' @importFrom visNetwork visSave
##' @importFrom visNetwork visInteraction
##' @importFrom visNetwork visOptions
##' @importFrom magrittr %>%
##' @author Kai Guo
ggnetwork_internal<-function (object=object,gene=gene,top = 50, pvalue = 0.05, padj = NULL,usePadj=TRUE,low = "orange",high = "red",
                              weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,cytoscapeFormat="graphml",
                              label.color = "black", label.size = 2,node.shape=NULL, layout = layout.fruchterman.reingold,savefig=FALSE,
                              visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename=NULL,
                              width=7,height=7,segment.size=0.2,node.alpha=0.7,sep=",",...){
  if (!is.null(padj)) {
    object <- object[object$Padj < padj, ]
  }
  else {
    object <- object[object$Pvalue < pvalue, ]
  }
  if (nrow(object) <= top) {
    object <- object
  }
  else {
    object <- object[1:top, ]
  }
  if (is.data.frame(gene)) {
    gene_p <- -log10(gene$padj)
    names(gene_p) <- rownames(gene)
  }
  else {
    gene_p <- rep(1, length(gene))
    names(gene_p) <- gene
  }
  if(isTRUE(usePadj)){
    pvalue<- -log10(object$Padj)

  }else{
    pvalue <- -log10(object$Pvalue)
  }
  names(pvalue) <- rownames(object)
  go2gen <- strsplit(x = as.vector(object$GeneID), split = sep)
  names(go2gen) <- rownames(object)
  gen2go <- reverseList(go2gen)
  golen <- object$Significant
  names(golen) <- rownames(object)
  gen2golen <- lapply(gen2go, function(x) golen[x])
  gen2gosum <- lapply(gen2golen, function(x) sum(x)/x)
  gen2res <- lapply(gen2gosum, function(x) x/sum(x))
  id <- rownames(object)
  n = nrow(object)
  w <- matrix(NA, nrow = n, ncol = n)
  colnames(w) <- rownames(w) <- rownames(object)
  for (i in 1:n) {
    ni <- id[i]
    for (j in i:n) {
      nj <- id[j]
      genein = intersect(go2gen[[ni]], go2gen[[nj]])
      geneup <- sum(gene_p[genein] * unlist(lapply(lapply(gen2res[genein],
                                                          "[", c(ni, nj)), sum)))
      genei <- setdiff(go2gen[[ni]], go2gen[[nj]])
      genej <- setdiff(go2gen[[nj]], go2gen[[ni]])
      geneid <- sum(gene_p[genei] * unlist(lapply(lapply(gen2res[genei],
                                                         "[", ni), sum)))
      genejd <- sum(gene_p[genej] * unlist(lapply(lapply(gen2res[genej],
                                                         "[", nj), sum)))
      gened <- geneup + geneid + genejd
      w[i, j] <- geneup/gened
    }
  }
  if (isTRUE(useTerm)) {
    colnames(w) <- rownames(w) <- object$Term
    names(pvalue) = object$Term
    rownames(object)<-object$Term
  }else{
    colnames(w) <- rownames(w) <- object$Annot
    names(pvalue) = object$Annot
    rownames(object)<-object$Annot
  }
  wn <- melt(w, as.is = TRUE)
  wn <- wn[wn[, 1] != wn[, 2], ]
  wn <- wn[!is.na(wn[, 3]), ]
  wn <- wn[wn[, 3] > 0, ]
  if(nrow(wn)==0){
    stop("All the terms share no overlap!")
  }
  g <- graph.data.frame(wn[, -3], directed = F)
  E(g)$width = sqrt(wn[, 3] * 5)
  pvalue = pvalue[V(g)$name]
  if (isTRUE(useTerm)) {
    idx <- unlist(sapply(V(g)$name, function(x) which(x == object$Term)))
  }else {
    idx <- unlist(sapply(V(g)$name, function(x) which(x == object$Annot)))
  }
  cols <- .color_scale(high, low)
  V(g)$color <- cols[sapply(pvalue, .getIdx, min(pvalue), max(pvalue))]
  if(sum(wn[,3]>weightcut)>=1){
     g <- delete.edges(g, E(g)[wn[, 3] < weightcut])
  }
  gs <- object$Significant
  if (isTRUE(useTerm)) {
    names(gs) <- object$Term
  }else {
    names(gs) <- object$Annot
  }
  V(g)$size <- log(gs[V(g)$name], base = 10) * 10
  if(!is.null(node.shape)){
    node.shape=rep(20,length(V(g)$name))
    names(node.shape)<-V(g)$name
    node.shape[names(node.shape)]<-node.shape
  }else{
    node.shape=rep(20,length(V(g)$name))
    names(node.shape)<-V(g)$name
  }
  if (isTRUE(writeCyt)) {
    if(is.null(cytoscapeFile)){
      cytoscapeFile="cytoscape-network.txt"
    }else{
      cytoscapeFile=paste0(cytoscapeFile,cytoscapeFormat)
    }
    write_graph(g, file = cytoscapeFile, format=cytoscapeFormat)
  }
  if (isTRUE(visNet)) {
    graph <- visIgraph(g, smooth = smooth)
    if(nodeselect==TRUE){
      graph<-graph%>%visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE)%>%
        visInteraction(navigationButtons = TRUE)
    }
    if(edit==TRUE){
      graph<-graph%>%visInteraction(navigationButtons = TRUE)%>%visOptions(manipulation = TRUE)
    }
    if(savehtml==TRUE){
      if(is.null(filename)){
        filename="network"
      }
      visSave(graph,file = paste0(filename,".html"))
    }
    graph
  }else{
  p<-ggnet2(g, node.size = V(g)$size, node.color = V(g)$color,
            edge.size = E(g)$width/20,node.shape=node.shape,node.alpha=node.alpha) +
    geom_text_repel(label = V(g)$name,
                    size=label.size,segment.size=segment.size,color=label.color)+
    theme(legend.position = "none")
  if(savefig==TRUE){
    if(is.null(filename)){
      filename="network"
    }
    ggsave(p,file=paste0(filename,".pdf"),width=width,height = height)
  }
  print(p)
  }
}

##' network for Enrichment results
##' @rdname ggnetwork
##' @param object richResult,GSEAResult object or dataframe
##' @param gene vector contains gene names or dataframe with DEGs information
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param weightcut cutoff valule for edge
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param layout layout for the network (layout.fruchterman.reingold)
##' @param low color used for small value
##' @param high color used for large value
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param cytoscapeFormat Character string giving the output file format
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param label.font label font
##' @param label.color label color
##' @param label.size label size
##' @param filename figure output name
##' @param savefig save the figure or not
##' @param width figure width
##' @param height figure height
##' @export
##' @author Kai Guo
setMethod("ggnetwork", signature(object = "richResult"),definition = function(object,top = 50, pvalue = 0.05, padj = NULL,usePadj=TRUE,low = "orange",high = "red",
                                                                              weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,cytoscapeFormat="graphml",
                                                                              label.color = "black", label.size = 2,node.shape=NULL, layout = layout.fruchterman.reingold,savefig=FALSE,
                                                                              visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename=NULL,
                                                                              width=7,height=7,segment.size=0.2,node.alpha=0.7,...) {
  ggnetwork_internal(object@result,gene=object@gene,top=top,pvalue=pvalue,padj=padj,usePadj=usePadj,weightcut=weightcut,useTerm=useTerm,writeCyt=writeCyt,cytoscapeFile=cytoscapeFile,cytoscapeFormat=cytoscapeFormat,
                     label.font=label.font,label.color=label.color,label.size=label.size,node.shape=node.shape,
                     layout=layout,savefig=savefig,width=width,height=height,
                     visNet=visNet,smooth=smooth,nodeselect=nodeselect,edit=edit,
                     filename=filename,node.alpha=node.alpha,sep=object@sep,...)
})
##' network for Enrichment results
##' @rdname ggnetwork
##' @param object richResult,GSEAResult object or dataframe
##' @param gene vector contains gene names or dataframe with DEGs information
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param weightcut cutoff valule for edge
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param layout layout for the network (layout.fruchterman.reingold)
##' @param low color used for small value
##' @param high color used for large value
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param cytoscapeFormat Character string giving the output file format
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param label.font label font
##' @param label.color label color
##' @param label.size label size
##' @param filename figure output name
##' @param savefig save the figure or not
##' @param width figure width
##' @param height figure height
##' @export
##' @author Kai Guo
setMethod("ggnetwork", signature(object = "data.frame"),definition = function(object,gene=NULL,top = 50, pvalue = 0.05, padj = NULL,usePadj=TRUE,low = "orange",high = "red",
                                                                              weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,cytoscapeFormat="graphml",
                                                                              label.color = "black", label.size = 2,node.shape=NULL, layout = layout.fruchterman.reingold,savefig=FALSE,
                                                                              visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename=NULL,
                                                                              width=7,height=7,segment.size=0.2,node.alpha=0.7,sep=",",...) {
  ggnetwork_internal(object,gene=gene,top=top,pvalue=pvalue,padj=padj,usePadj=usePadj,weightcut=weightcut,useTerm=useTerm,writeCyt=writeCyt,cytoscapeFile=cytoscapeFile,cytoscapeFormat=cytoscapeFormat,
                     label.font=label.font,label.color=label.color,label.size=label.size,node.shape=node.shape,
                     layout=layout,savefig=savefig,width=width,height=height,
                     visNet=visNet,smooth=smooth,nodeselect=nodeselect,edit=edit,
                     filename=filename,node.alpha=node.alpha,sep=sep,...)
})
#
##' network for Enrichment results
##' @rdname ggnetwork
##' @param object data.frame object
##' @param gene vector contains gene names or dataframe with DEGs information
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param weightcut cutoff valule for edge
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param layout layout for the network (layout.fruchterman.reingold)
##' @param low color used for small value
##' @param high color used for large value
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param cytoscapeFormat Character string giving the output file format
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param label.font label font
##' @param label.color label color
##' @param label.size label size
##' @param filename figure output name
##' @param savefig save the figure or not
##' @param width figure width
##' @param height figure height
##' @export
##' @author Kai Guo
setMethod("ggnetwork", signature(object = "GSEAResult"),definition = function(object,gene,top = 50, pvalue = 0.05, padj = NULL,usePadj=TRUE,low = "orange",high = "red",
                                                                              weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,cytoscapeFormat="graphml",
                                                                              label.color = "black", label.size = 2,node.shape=NULL, layout = layout.fruchterman.reingold,savefig=FALSE,
                                                                              visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename="network",
                                                                              width=7,height=7,segment.size=0.2,node.alpha=0.7,...) {
  sep=object@sep
  if(is.list(object@result$leadingEdge)){
    object@result$leadingEdge<-unlist(lapply(object@result$leadingEdge, function(x)paste(x,collapse = sep, sep="")))
  }
  object@result$Annot<-object@result$pathway
  object@result<-object@result[,c(9,1:3,7,8)]
  colnames(object@result)<-c("Annot","Term","Pvalue","Padj","Significant","GeneID")
  ggnetwork_internal(object@result,gene=object@gene,top=top,pvalue=pvalue,padj=padj,usePadj=usePadj,weightcut=weightcut,useTerm=useTerm,writeCyt=writeCyt,cytoscapeFile=cytoscapeFile,cytoscapeFormat=cytoscapeFormat,
                     label.font=label.font,label.color=label.color,label.size=label.size,node.shape=node.shape,
                     layout=layout,savefig=savefig,width=width,height=height,
                     visNet=visNet,smooth=smooth,nodeselect=nodeselect,edit=edit,
                     filename=filename,node.alpha=node.alpha,sep=sep,...)
})



