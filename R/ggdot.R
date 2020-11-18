##' Dotplot for enrichment results
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_point
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 geom_text
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 scale_color_gradient
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ggsave
##' @importFrom ggplot2 theme_minimal
##' @importFrom ggplot2 labs
##' @importFrom ggplot2 guides
##' @importFrom ggplot2 guide_colourbar
##' @importFrom ggplot2 guide_legend
##' @param object richResult object
##' @param top number of terms you want to display,
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param low low color
##' @param high high color
##' @param alpha transparency alpha
##' @param font.x font of x axis
##' @param font y font of y axis
##' @param fontsize.x fontsize of x axis
##' @param fontsize.y fontsize of y axis
##' @param short automatic short name or not
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param font.size font size for xlim or ylim
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
ggdot_internal<-function(object,top=50,pvalue=0.05,order=FALSE,
                         low="lightpink",high="red",alpha=0.7,
                         font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                         short=FALSE,
                         padj=NULL,usePadj=TRUE,filename=NULL,width=10,height=8){
  if(!is.null(padj)){
    object<-object[object$Padj<padj,]
  }else{
    object<-object[object$Pvalue<pvalue,]
  }
  if(nrow(object)>=top){
    dd<-object[1:top,]
  }else{
    dd<-object
  }
  if(nrow(dd)>=1){
    dd[,3]<-dd[,4]/dd[,3]
    colnames(dd)[3]<-"rich";
    if(isTRUE(short)){
      dd$Term<-unlist(lapply(dd$Term,function(x).paste.char(x)))
    }
    if(order==TRUE){
      dd$Term<-factor(dd$Term,levels=dd$Term[order(dd$rich)])
    }
    if(usePadj==FALSE){
      p<-ggplot(dd,aes(x=rich,y=Term))+geom_point(aes(size=Significant,color=-log10(Pvalue)),alpha=alpha)+theme_minimal()+
        theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x))+
        scale_color_gradient(low=low,high=high)+ylab("Pathway name")+
        xlab("Rich factor")+labs(size="Gene number")+guides(color=guide_colourbar(order = 1),size=guide_legend(order = 2))
    }else{
      p<-ggplot(dd,aes(x=rich,y=Term))+geom_point(aes(size=Significant,color=-log10(Padj)),alpha=alpha)+theme_minimal()+
        theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x))+
        scale_color_gradient(low=low,high=high)+ylab("Pathway name")+
        xlab("Rich factor")+labs(size="Gene number")+guides(color=guide_colourbar(order = 1),size=guide_legend(order = 2))
    }}else{
      cat("No Pathway enrichment results were found!\n")
    }
  if(!is.null(filename)){
    ggsave(p,file=paste(filename,"dot.pdf",sep="_"),width=width,height=height)
  }
  p
}
##' dotplot for Enrichment results
##' @rdname ggdot
##' @param object richResult object
##' @param top number of terms you want to display,
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param low low color
##' @param high high color
##' @param alpha transparency alpha
##' @param font.x font of x axis
##' @param font y font of y axis
##' @param fontsize.x fontsize of x axis
##' @param fontsize.y fontsize of y axis
##' @param short automatic short name or not
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param font.size font size for xlim or ylim
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#'   ggdot(res)
#' }
##' @exportMethod ggdot
##' @author Kai Guo
setMethod("ggdot", signature(object = "richResult"),definition = function(object,top=50,pvalue=0.05,order=FALSE,
                                                                          low="lightpink",high="red",alpha=0.7,
                                                                          font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                                                                          short=FALSE,
                                                                          padj=NULL,usePadj=TRUE,filename=NULL,width=10,height=8) {
            ggdot_internal(object@result,top=top,pvalue=pvalue,order=order,
                           low=low,high=high,alpha=alpha,
                           font.x=font.x,font.y=font.y,fontsize.x=fontsize.x,fontsize.y=fontsize.y,
                           padj=padj,usePadj=usePadj,filename=filename,width=width,height=height)
          })
##' dotplot for Enrichment results
##' @rdname ggdot
##' @param object richResult object
##' @param top number of terms you want to display,
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param low low color
##' @param high high color
##' @param alpha transparency alpha
##' @param font.x font of x axis
##' @param font y font of y axis
##' @param fontsize.x fontsize of x axis
##' @param fontsize.y fontsize of y axis
##' @param short automatic short name or not
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param font.size font size for xlim or ylim
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @examples
##' \dontrun{
##'   hsago<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
##'   gene=sample(unique(hsago$GeneID),1000)
##'   res<-richKEGG(gene,kodata = hsako)
##'   ggdot(result(res))
##' }
##' @exportMethod ggdot
##' @author Kai Guo
setMethod("ggdot", signature(object = "data.frame"),definition = function(object,top=50,pvalue=0.05,order=FALSE,
                                                                          low="lightpink",high="red",alpha=0.7,
                                                                          font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                                                                          short=FALSE,
                                                                          padj=NULL,usePadj=TRUE,filename=NULL,width=10,height=8) {
            ggdot_internal(object,top=top,pvalue=pvalue,order=order,
                           low=low,high=high,alpha=alpha,
                           font.x=font.x,font.y=font.y,fontsize.x=fontsize.x,fontsize.y=fontsize.y,
                           padj=padj,usePadj=usePadj,filename=filename,width=width,height=height)
          })
