##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 geom_bar
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 geom_text
##' @importFrom ggplot2 theme
##' @importFrom ggplot2 scale_fill_gradient
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 ylim
##' @importFrom ggplot2 ggsave
##' @importFrom ggplot2 theme_light
##' @importFrom ggplot2 labs
##' @rdname richBar
##' @param resultFis data frame of enrichment results
##' @param top number of terms you want to display,
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param order order by Term or richFactor
##' @param horiz horiz or not
##' @param low low color
##' @param high high color
##' @param font.x font of x axis
##' @param font.y font of y axis
##' @param fontsize.x fontsize of x axis
##' @param fontsize.y fontsize of y axis
##' @param short automatic short name or not
##' @param fontsize.text fontsize for bar text labels
##' @param angle angle for x axis labels
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param orderp order by p value(adjusted p value)
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
richBar_internal<-function(resultFis,top=50,pvalue=0.05,order=FALSE,horiz=TRUE,
                         low="lightpink",high="red",
                         font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                         short=FALSE,
                         fontsize.text=3,angle=75,padj=NULL,usePadj=TRUE,orderp=FALSE,
                         filename=NULL,width=10,height=8){
  if(!is.null(padj)){
    resultFis<-resultFis[resultFis$Padj<padj,]
  }else{
    resultFis<-resultFis[resultFis$Pvalue<pvalue,]
  }
  if(nrow(resultFis)>=top){
    resultFis<-resultFis[1:top,]
  }
  if(max(resultFis$Significant/(resultFis$Annotated+0.1))<=1){
    yheight=max(resultFis$Significant/resultFis$Annotated)+0.1
  }else{
    yheight=1
  }
  if(isTRUE(short)){
    resultFis$Term<-unlist(lapply(resultFis$Term,function(x).paste.char(x,n=6)))
  }
  if(isTRUE(order)){
    resultFis$rich<-as.numeric(resultFis$Significant)/as.numeric(resultFis$Annotated)
    if(isTRUE(orderp)){
      resultFis$Term<-factor(resultFis$Term,levels=resultFis$Term[order(resultFis$Pvalue)])
    }else{
      resultFis$Term<-factor(resultFis$Term,levels=resultFis$Term[order(resultFis$rich)])
    }
  }
  if(isTRUE(usePadj)){
    resultFis$neg_log_p <- -log10(as.numeric(resultFis$Padj))
    fill_label <- "-log10(Padj)"
  }else{
    resultFis$neg_log_p <- -log10(as.numeric(resultFis$Pvalue))
    fill_label <- "-log10(Pvalue)"
  }
  p<-ggplot(resultFis,aes(x=Term,y=round(as.numeric(Significant/Annotated),2)))+
    geom_bar(stat="identity",aes(fill=neg_log_p))+
    scale_fill_gradient(low=low,high=high)+theme_light()
  if(isTRUE(horiz)){
    angle<-0
    p<-p+theme(axis.text.y=element_text(face=font.y,size=fontsize.y),
               axis.text.x=element_text(face=font.x,color="black",size=fontsize.x,angle=angle))+
      labs(fill=fill_label)+coord_flip()+
      geom_text(aes(label=Significant),hjust=-0.3,size=fontsize.text)+
      xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)
  }else{
    p<-p+theme(axis.text.y=element_text(face=font.y,size=fontsize.y),
               axis.text.x=element_text(face=font.x,color="black",size=fontsize.x,angle=angle,vjust=1,hjust=1))+
      labs(fill=fill_label)+
      geom_text(aes(label=Significant),vjust=-0.3,size=fontsize.text)+
      xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)
  }
  if(!is.null(filename)){
    ggsave(p,filename=paste(filename,"enrich.pdf",sep="_"),width=width,height=height)
  }
  p
}
##' barplot for Enrichment results
##' @rdname richBar
##' @param object richResult object
##' @param ... additional arguments
##' @examples
##' \dontrun{
##'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
##'   gene=sample(unique(hsako$GeneID),1000)
##'   res<-richKEGG(gene,kodata = hsako)
##'   richBar(res)
##' }
##' @export
##' @author Kai Guo
setMethod("richBar", signature(object = "richResult"),definition = function(object,top=50,pvalue=0.05,padj=NULL,order=FALSE,
                   usePadj=TRUE,fontsize.x=10,fontsize.y=10,short=FALSE,fontsize.text=3,angle=75,orderp=FALSE,filename=NULL,
                   width=10,height=8,horiz=TRUE,...) {
            richBar_internal(object@result,top=top,pvalue=pvalue,padj=padj,order=order,
                           usePadj=usePadj,fontsize.x=fontsize.x,fontsize.y=fontsize.y,short=short,fontsize.text = fontsize.text,angle=angle,
                           orderp=orderp,filename=filename,horiz=horiz, ...)
          })
##' barplot for Enrichment result
##' @rdname richBar
##' @examples
##' \dontrun{
##'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
##'   gene=sample(unique(hsako$GeneID),1000)
##'   res<-richKEGG(gene,kodata = hsako)
##'   richBar(result(res))
##' }
##' @export
##' @author Kai Guo
setMethod("richBar", signature(object = "data.frame"),definition = function(object,top=50,pvalue=0.05,padj=NULL,order=FALSE,
                                                                          usePadj=TRUE,fontsize.x=10,fontsize.y=10,short=FALSE,fontsize.text=3,angle=75,orderp=FALSE,filename=NULL,
                                                                          width=10,height=8,horiz=TRUE,...) {
  richBar_internal(object,top=top,pvalue=pvalue,padj=padj,order=order,
                 usePadj=usePadj,fontsize.x=fontsize.x,fontsize.y=fontsize.y,short=short,fontsize.text = fontsize.text,angle=angle,
                 orderp=orderp,filename=filename,width=width,height=height,horiz=horiz,...)
          })
