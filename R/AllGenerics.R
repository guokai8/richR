##' richGO
##'
##' @name richGO
##' @rdname richGO-methods
##' @title richGO method
##' @param x vector contains gene names
##' @param godata Annotation object or dataframe
##' @param ... additional parameters
##' @return richResult
##' @export
##' @author Kai Guo
setGeneric("richGO", function(x,godata,ontology="BP",pvalue=0.05,padj=NULL,minSize=2,maxSize=500,
                              keepRich=TRUE,filename=NULL,padj.method="BH",...)
  standardGeneric("richGO"))

##' richKEGG
##'
##' @name richKEGG
##' @rdname richKEGG-methods
##' @title richKEGG method
##' @param x vector contains gene names
##' @param kodata Annotation object or dataframe
##' @param ... additional parameters
##' @return richResult
##' @export
##' @author Kai Guo
setGeneric("richKEGG", function(x,kodata,pvalue=0.05,padj=NULL,organism=NULL,ontology="KEGG",
                                keytype=NULL,minSize=2,maxSize=500,
                                keepRich=TRUE,filename=NULL,padj.method="BH",builtin=TRUE,...)
  standardGeneric("richKEGG"))

##' richGSEA
##'
##' @name richGSEA
##' @rdname richGSEA-methods
##' @title richGSEA method
##' @param x vector contains gene names
##' @param Annot object or dataframe
##' @param ... additional parameters
##' @return GSEAResult
##' @export
##' @author Kai Guo
setGeneric("richGSEA", function(x,object,keytype="",pvalue=0.05,padj=NULL,minSize=15,ontology="",
                                maxSize=500,nperm=5000,filename=NULL,padj.method="BH",organism=NULL,...)
  standardGeneric("richGSEA"))
##' enrich
##'
##' @name enrich
##' @rdname enrich-methods
##' @title enrich method
##' @param x vector contains gene names
##' @param annot Annotation object or dataframe
##' @param ... additional parameters
##' @return richResult
##' @export
##' @author Kai Guo
setGeneric("enrich", function(x,annot,pvalue=0.05,...)
  standardGeneric("enrich"))


##' ggdot
##'
##' @name ggdot
##' @rdname ggdot-methods
##' @title ggdot method
##' @param object richResult or dataframe
##' @param ... additional parameters
##' @return plot
##' @export
##' @author Kai Guo
setGeneric("ggdot", function(object,top=50,pvalue=0.05,order=FALSE,
                             low="lightpink",high="red",alpha=0.7,
                             font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                             padj=NULL,usePadj=TRUE,filename=NULL,width=10,height=8,...)
    standardGeneric("ggdot"))

##' ggbar
##'
##' @name ggbar
##' @rdname ggbar-methods
##' @title ggbar method
##' @param object richResult or dataframe
##' @param ... additional parameters
##' @return plot
##' @export
##' @author Kai Guo
setGeneric("ggbar", function(object,top=50,pvalue=0.05,padj=NULL,order=FALSE,
                             usePadj=TRUE,...)
    standardGeneric("ggbar"))
##' ggnetplot
##' @rdname ggnetplot-method
##' @title ggnetplot method
##' @param object richResult or dataframe
##' @param ... additional parameters
##' @return plot
##' @export
##' @author Kai Guo
setGeneric("ggnetplot",function(object,top=50, pvalue=0.05, padj=NULL,
                                usePadj =TRUE, useTerm=TRUE,low="orange",high="red",
                                writeCyt=FALSE, cytoscapeFile="network-file-for-cytoscape.txt",
                                label.color = "black", label.size = 2, node.shape=NULL,
                                layout = layout.fruchterman.reingold,savefig=FALSE,filename="network",
                                width=7,height=7,node.alpha=0.7,repel=TRUE,segment.size=0.2,...)
  standardGeneric("ggnetplot"))


##' ggnetwork
##'
##' @name ggnetwork
##' @rdname ggnetwork-methods
##' @title ggnetwork method
##' @param object richResult object or dataframe
##' @param ... additional paramters
##' @return plot
##' @export
setGeneric("ggnetwork", function(object,gene,top = 50, pvalue = 0.05, padj = NULL,low = "orange",high = "red",
                                 weightcut = 0.2, useTerm = TRUE, writeCyt = FALSE,cytoscapeFile = "network-file-for-cytoscape.txt",
                                 label.color = "black", label.size = 2,node.shape=NULL, layout = layout.fruchterman.reingold,savefig=FALSE,
                                 visNet=FALSE,smooth=TRUE,nodeselect=FALSE,edit=FALSE,savehtml=FALSE,filename="network",
                                 width=7,height=7,segment.size=0.2,node.alpha=0.7,...)
  standardGeneric("ggnetwork"))


##' result generic
##' @param x richResult object
##' @return result return dataframe and print summary
##' @export
result<-function(x){
  UseMethod("result",x)
}
##' detail generic
##' @param x richResult object
##' @return detail return detial for these significant genes
##' @export
detail<-function(x){
  UseMethod("detail",x)
}
