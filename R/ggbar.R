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
ggbar_internal<-function(resultFis,top=50,pvalue=0.05,order=FALSE,horiz=FALSE,
                         low="lightpink",high="red",
                         font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                         fontsize.text=3,angle=75,padj=NULL,usePadj=TRUE,
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
  if(order==TRUE){
    resultFis$rich<-as.numeric(resultFis$Significant)/as.numeric(resultFis$Annotated)
    resultFis$Term<-factor(resultFis$Term,levels=resultFis$Term[order(resultFis$rich)])
  }
  if(usePadj==FALSE){
    p<-ggplot(resultFis,aes(x=Term,y=round(as.numeric(Significant/Annotated),2)))+geom_bar(stat="identity",aes(fill=-log10(as.numeric(Pvalue))))
    p<-p+scale_fill_gradient(low=low,high=high)+theme_light()
    if(horiz==TRUE){
      p<-p+theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x,angle=angle))+labs(fill="-log10(Pvalue)")
      p<-p+coord_flip()
      p<-p+geom_text(aes(label=Significant),hjust=-0.3,size=fontsize.text)+xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)
    }else{
      p<-p+theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x,angle=angle,vjust=1,hjust=1))+labs(fill="-log10(Pvalue)")
      p<-p+geom_text(aes(label=Significant),vjust=-0.3,size=fontsize.text)+xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)
    }
  #  print(p)
  }else{
    p<-ggplot(resultFis,aes(x=Term,y=round(as.numeric(Significant/Annotated),2)))+geom_bar(stat="identity",aes(fill=-log10(as.numeric(Padj))))
    p<-p+scale_fill_gradient(low=low,high=high)+theme_light()
    if(horiz==TRUE){
      p<-p+theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x,angle=angle))+labs(fill="-log10(Padj)")
      p<-p+coord_flip()
      p<-p+geom_text(aes(label=Significant),hjust=-0.3,size=fontsize.text)+xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)

    }else{
      p<-p+theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x,angle=angle,vjust=1,hjust=1))+labs(fill="-log10(Padj)")
      p<-p+geom_text(aes(label=Significant),vjust=-0.3,size=fontsize.text)+xlab("Annotation")+ylab("Rich Factor")+ylim(0,yheight)
    }
  }
  if(!is.null(filename)){
    ggsave(p,file=paste(filename,OP,"enrich.pdf",sep="_"),width=width,height=height)
  }
  p
}
##' barplot for Enrichment results
##' @rdname ggbar
##' @param object richResult object
##' @param top number of terms you want to display,
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param font.size font size for xlim or ylim
##' @param angle label angle of xlim
##' @param filename figure output name
##' @exportMethod ggbar
##' @author Kai Guo
setMethod("ggbar", signature(object = "richResult"),definition = function(object,top=50,pvalue=0.05,padj=NULL,order=FALSE,
                   usePadj=TRUE,fontsize.x=10,fontsize.y=10,fontsize.text=3,angle=75,filename=NULL,
                   width=10,height=8,...) {
            ggbar_internal(object@result,top=top,pvalue=pvalue,padj=padj,order=order,
                           usePadj=usePadj,fontsize.x=fontsize.x,fontsize.y=fontsize.y,fontsize.text = fontsize.text,angle=angle,filename=filename, ...)
          })
##' barplot for Enrichment result
##' @rdname ggbar
##' @param object dataframe of enrichment results
##' @param top number of terms you want to display,
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param font.size font size for xlim or ylim
##' @param angle label angle of xlim
##' @param filename figure output name
setMethod("ggbar", signature(object = "data.frame"),definition = function(object,top=50,pvalue=0.05,padj=NULL,order=FALSE,
                                                                          usePadj=TRUE,fontsize.x=10,fontsize.y=10,fontsize.text=3,angle=75,filename=NULL,
                                                                          width=10,height=8,...) {
  ggbar_internal(object,top=top,pvalue=pvalue,padj=padj,order=order,
                 usePadj=usePadj,fontsize.x=fontsize.x,fontsize.y=fontsize.y,fontsize.text = fontsize.text,angle=angle,filename=filename, ...)
          })
