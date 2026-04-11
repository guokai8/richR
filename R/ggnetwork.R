##' generate network based on Enrichment results
##' @rdname richNetwork
##' @param object richResult,GSEAResult object or dataframe
##' @param gene vector contains gene names or dataframe with DEGs information
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param low color used for small value
##' @param high color used for large value
##' @param weightcut cutoff value for edge weight
##' @param useTerm use Term name instead of Annotation ID
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param cytoscapeFormat Character string giving the output file format
##' @param label.color label color
##' @param label.size label size
##' @param node.shape shape of the nodes
##' @param layout layout method ('fruchtermanreingold','kamadakawai','target','circle')
##' @param savefig save the figure or not
##' @param visNet use visNetwork for interactive plot
##' @param smooth smooth edges in visNetwork
##' @param nodeselect enable node selection in visNetwork
##' @param edit enable editing in visNetwork
##' @param savehtml save visNetwork as HTML file
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param sep character string used to separate the genes when concatenating
##' @param ... additional arguments
##' @importFrom igraph graph.data.frame
##' @importFrom igraph delete.edges
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 ggsave
##' @importFrom igraph E
##' @importFrom igraph E<-
##' @importFrom igraph V
##' @importFrom igraph V<-
##' @importFrom igraph write_graph
##' @author Kai Guo
ggnetwork_internal<-function (object=object,gene=gene,top = 50, pvalue = 0.05, padj = NULL,usePadj=TRUE,low = "orange",high = "red",
                              weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,cytoscapeFormat="graphml",
                              label.color = "black", label.size = 2,node.shape=NULL, layout = "fruchtermanreingold",savefig=FALSE,
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
  wn <- as.data.frame(as.table(w), stringsAsFactors = FALSE)
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
    shape=rep(20,length(V(g)$name))
    names(shape)<-V(g)$name
    shape[names(node.shape)]<-node.shape
    node.shape <- shape
  }else{
    node.shape=rep(20,length(V(g)$name))
    names(node.shape)<-V(g)$name
  }
  if (isTRUE(writeCyt)) {
    if(is.null(cytoscapeFile)){
      cytoscapeFile="cytoscape-network.txt"
    }else{
      cytoscapeFile=paste0(cytoscapeFile,".",cytoscapeFormat)
    }
    if(cytoscapeFormat=="edgelist"){
      write.table(igraph::get.data.frame(g),file=cytoscapeFile,sep="\t",row.names=FALSE)
    }else{
      write_graph(g, file = cytoscapeFile, format=cytoscapeFormat)
    }
  }
  if (isTRUE(visNet)) {
    if (!requireNamespace("visNetwork", quietly = TRUE))
      stop("Package 'visNetwork' is required for interactive network plots. Install it with install.packages('visNetwork').")
    graph <- visNetwork::visIgraph(g, smooth = smooth)
    if(nodeselect==TRUE){
      graph <- visNetwork::visOptions(graph, highlightNearest = TRUE, nodesIdSelection = TRUE)
      graph <- visNetwork::visInteraction(graph, navigationButtons = TRUE)
    }
    if(edit==TRUE){
      graph <- visNetwork::visInteraction(graph, navigationButtons = TRUE)
      graph <- visNetwork::visOptions(graph, manipulation = TRUE)
    }
    if(savehtml==TRUE){
      if(is.null(filename)){
        filename="network"
      }
      visNetwork::visSave(graph,file = paste0(filename,".html"))
    }
    graph
  }else{
  if (!requireNamespace("GGally", quietly = TRUE))
    stop("Package 'GGally' is required for network plots. Install it with install.packages('GGally').")
  if (!requireNamespace("intergraph", quietly = TRUE))
    stop("Package 'intergraph' is required for network plots. Install it with install.packages('intergraph').")
  if (!requireNamespace("ggrepel", quietly = TRUE))
    stop("Package 'ggrepel' is required for network plots. Install it with install.packages('ggrepel').")
  p<-GGally::ggnet2(g, node.size = V(g)$size, node.color = V(g)$color,mode=layout,
            edge.size = E(g)$width/20,node.shape=node.shape,node.alpha=node.alpha) +
    ggrepel::geom_text_repel(label = V(g)$name,
                    size=label.size,segment.size=segment.size,color=label.color)+
    theme(legend.position = "none")
  if(savefig==TRUE){
    if(is.null(filename)){
      filename="network"
    }
    ggsave(p,filename=paste0(filename,".pdf"),width=width,height = height)
  }
  p
  }
}

##' network for Enrichment results
##' @rdname richNetwork-methods
##' @param object richResult object
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not
##' @param low color used for small value
##' @param high color used for large value
##' @param weightcut cutoff value for edge weight
##' @param useTerm use Term name instead of Annotation ID
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param cytoscapeFormat output file format
##' @param label.color label color
##' @param label.size label size
##' @param node.shape shape of the nodes
##' @param layout layout method
##' @param savefig save the figure or not
##' @param visNet use visNetwork for interactive plot
##' @param smooth smooth edges in visNetwork
##' @param nodeselect enable node selection in visNetwork
##' @param edit enable editing in visNetwork
##' @param savehtml save visNetwork as HTML file
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param ... additional arguments
##' @export
##' @author Kai Guo
setMethod("richNetwork", signature(object = "richResult"),definition = function(object,top = 50, pvalue = 0.05, padj = NULL,usePadj=TRUE,low = "orange",high = "red",
                                                                              weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,cytoscapeFormat="graphml",
                                                                              label.color = "black", label.size = 2,node.shape=NULL, layout = "fruchtermanreingold",savefig=FALSE,
                                                                              visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename=NULL,
                                                                              width=7,height=7,segment.size=0.2,node.alpha=0.7,...) {
  ggnetwork_internal(object@result,gene=object@gene,top=top,pvalue=pvalue,padj=padj,usePadj=usePadj,weightcut=weightcut,useTerm=useTerm,
                     writeCyt=writeCyt,cytoscapeFile=cytoscapeFile,cytoscapeFormat=cytoscapeFormat,
                     label.color=label.color,label.size=label.size,node.shape=node.shape,
                     layout=layout,savefig=savefig,width=width,height=height,
                     visNet=visNet,smooth=smooth,nodeselect=nodeselect,edit=edit,
                     filename=filename,node.alpha=node.alpha,sep=object@sep,...)
})
##' network for Enrichment results
##' @rdname richNetwork-methods
##' @param object data.frame of enrichment results
##' @param gene vector contains gene names
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not
##' @param low color used for small value
##' @param high color used for large value
##' @param weightcut cutoff value for edge weight
##' @param useTerm use Term name instead of Annotation ID
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param cytoscapeFormat output file format
##' @param label.color label color
##' @param label.size label size
##' @param node.shape shape of the nodes
##' @param layout layout method
##' @param savefig save the figure or not
##' @param visNet use visNetwork for interactive plot
##' @param smooth smooth edges in visNetwork
##' @param nodeselect enable node selection in visNetwork
##' @param edit enable editing in visNetwork
##' @param savehtml save visNetwork as HTML file
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param sep character string used to separate the genes
##' @param ... additional arguments
##' @export
##' @author Kai Guo
setMethod("richNetwork", signature(object = "data.frame"),definition = function(object,gene=NULL,top = 50, pvalue = 0.05, padj = NULL,usePadj=TRUE,low = "orange",high = "red",
                                                                              weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,cytoscapeFormat="graphml",
                                                                              label.color = "black", label.size = 2,node.shape=NULL, layout = "fruchtermanreingold",savefig=FALSE,
                                                                              visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename=NULL,
                                                                              width=7,height=7,segment.size=0.2,node.alpha=0.7,sep=",",...) {
  ggnetwork_internal(object,gene=gene,top=top,pvalue=pvalue,padj=padj,usePadj=usePadj,weightcut=weightcut,useTerm=useTerm,writeCyt=writeCyt,cytoscapeFile=cytoscapeFile,
                     cytoscapeFormat=cytoscapeFormat,
                     label.color=label.color,label.size=label.size,node.shape=node.shape,
                     layout=layout,savefig=savefig,width=width,height=height,
                     visNet=visNet,smooth=smooth,nodeselect=nodeselect,edit=edit,
                     filename=filename,node.alpha=node.alpha,sep=sep,...)
})
#
##' network for Enrichment results
##' @rdname richNetwork-methods
##' @param object GSEAResult object
##' @param gene vector contains gene names
##' @param top number of terms to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not
##' @param low color used for small value
##' @param high color used for large value
##' @param weightcut cutoff value for edge weight
##' @param useTerm use Term name instead of Annotation ID
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
##' @param cytoscapeFormat output file format
##' @param label.color label color
##' @param label.size label size
##' @param node.shape shape of the nodes
##' @param layout layout method
##' @param savefig save the figure or not
##' @param visNet use visNetwork for interactive plot
##' @param smooth smooth edges in visNetwork
##' @param nodeselect enable node selection in visNetwork
##' @param edit enable editing in visNetwork
##' @param savehtml save visNetwork as HTML file
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @param segment.size size for label segment
##' @param node.alpha alpha-transparency scales
##' @param ... additional arguments
##' @export
##' @author Kai Guo
setMethod("richNetwork", signature(object = "GSEAResult"),definition = function(object,gene,top = 50, pvalue = 0.05, padj = NULL,usePadj=TRUE,low = "orange",high = "red",
                                                                              weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = NULL,cytoscapeFormat="graphml",
                                                                              label.color = "black", label.size = 2,node.shape=NULL, layout = "fruchtermanreingold",savefig=FALSE,
                                                                              visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename="network",
                                                                              width=7,height=7,segment.size=0.2,node.alpha=0.7,...) {
  sep=object@sep
  if(is.list(object@result$leadingEdge)){
    object@result$leadingEdge<-unlist(lapply(object@result$leadingEdge, function(x)paste(x,collapse = sep, sep="")))
  }
  object@result$Annot<-object@result$pathway
  object@result<-object@result[,c(9,1:3,7,8)]
  colnames(object@result)<-c("Annot","Term","Pvalue","Padj","Significant","GeneID")
  ggnetwork_internal(object@result,gene=object@gene,top=top,pvalue=pvalue,padj=padj,usePadj=usePadj,weightcut=weightcut,useTerm=useTerm,
                     writeCyt=writeCyt,cytoscapeFile=cytoscapeFile,cytoscapeFormat=cytoscapeFormat,
                     label.color=label.color,label.size=label.size,node.shape=node.shape,
                     layout=layout,savefig=savefig,width=width,height=height,
                     visNet=visNet,smooth=smooth,nodeselect=nodeselect,edit=edit,
                     filename=filename,node.alpha=node.alpha,sep=sep,...)
})



