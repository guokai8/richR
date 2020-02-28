##' generate network based on Enrichment results
##' @rdname ggnetwork
##' @param richRes list of enrichment object
##' @param gene vector contains gene names or dataframe with DEGs information
##' @param top number of terms to display
##' @param top.display top number to display
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param weightcut cutoff valule for edge
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param layout layout for the network (layout.fruchterman.reingold)
##' @param low color used for small value
##' @param high color used for large value
##' @param writeCyt write out the cytoscape file
##' @param cytoscapeFile output cytoscape File
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
ggnetmap<-function(richRes,gene=NULL,top=50,top.display=NULL,pvalue = 0.05, padj = NULL,low = "orange",high = "red",
                   weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = "network-file-for-cytoscape.txt",
                   label.color = "black", label.size = 2,node.shape=NULL, layout = layout.fruchterman.reingold,savefig=FALSE,
                   visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename="network",
                   width=7,height=7,segment.size=0.2,node.alpha=0.7,...){
  if(!is.null(top)){
    object<-Reduce(function(x, y) rbind(x, y), lapply(richRes, function(x)x@result[1:top,]))
  }else{
    object<-Reduce(function(x, y) rbind(x, y), lapply(richRes, function(x)x@result))
  }
  object <- object[order(object$Pvalue),]
  if(is.null(gene)){
    gene <- unique(unlist(lapply(richRes, function(x) x@gene)))
  }
  if(is.null(top.display)){
    top.display=nrow(object)
  }
  ggnetwork(object,gene=gene,top=top,pvalue=pvalue,padj=padj,weightcut=weightcut,useTerm=useTerm,writeCyt=writeCyt,
            label.font=label.font,label.color=label.color,label.size=label.size,node.shape=node.shape,
            layout=layout,savefig=savefig,width=width,height=height,
            visNet=visNet,smooth=smooth,nodeselect=nodeselect,edit=edit,
            filename=filename,node.alpha=node.alpha,...)
}
