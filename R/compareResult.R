#' compare enrichment results across different samples
#' @param x list of richResults
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @examples
#' \dontrun{
#' hsako <- buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#' gene1 <- sample(unique(hsako$GeneID),1000)
#' gene2 <- sample(unique(hsako$GeneID),1000)
#' resko1 <-richKEGG(gene1,kodata = hsako)
#' resko2 <-richKEGG(gene2,kodata = hsako)
#' res<-compareResult(list(S1=resko1,S2=resko2))
#' }
#' @author Kai Guo
#' @export
compareResult<-function(x,pvalue=0.05,padj=NULL){
  if(!is.null(padj)){
    pvalue <- 1
  }else{
    padj <- 1
  }
  if(is.null(names(x))){
    names(x) <- paste("Group",  seq_along(x))
  }
  tmp <- lapply(x, function(x)filter(x,Padj<padj,Pvalue<pvalue))
  tmp <- lapply(names(x), function(x)mutate(tmp[[x]],group=x))
  dx<-do.call(rbind,tmp)
  return(dx)
}
##' draw dotplot for multiple enrichment results
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
##' @param x dataframe of enrichment result
##' @param pvalue cutoff value of pvalue (if padj set as NULL)
##' @param low low color
##' @param high high color
##' @param alpha transparency alpha
##' @param font.x font of x axis
##' @param font.y font of y axis
##' @param fontsize.x fontsize of x axis
##' @param fontsize.y fontsize of y axis
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @examples
#' \dontrun{
#' hsako <- buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#' gene1 <- sample(unique(hsako$GeneID),1000)
#' gene2 <- sample(unique(hsako$GeneID),1000)
#' resko1 <-richKEGG(gene,kodata = hsako)
#' resko2 <-richKEGG(gene,kodata = hsako)
#' res<-compareResult(list(S1=resko1,S2=resko2))
#' comparedot(res,pvalue=0.05)
#' }
#' @author Kai Guo
#' @export
comparedot <- function(x,pvalue=0.05,
                       low="lightpink",high="red",alpha=0.7,
                       font.x="bold",font.y="bold",fontsize.x=10,fontsize.y=10,
                       padj=NULL,usePadj=TRUE,filename=NULL,width=10,height=8){
  if(!is.null(padj)){
    x<-x[x$Padj<padj,]
  }else{
    x<-x[x$Pvalue<pvalue,]
  }
  x$Term <- unlist(lapply(x$Term,function(x).paste.char(x)))
  if(isTRUE(usePadj)){
    p<-ggplot(x,aes(x=group,y=Term))+geom_point(aes(size=Significant,color=-log10(Padj)),alpha=alpha)+theme_minimal()+
      theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x))+
      scale_color_gradient(low=low,high=high)+ylab("Pathway name")+
      xlab("")+labs(size="Gene number")+guides(color=guide_colourbar(order = 1),size=guide_legend(order = 2))
  }else{
    p<-ggplot(x,aes(x=group,y=Term))+geom_point(aes(size=Significant,color=-log10(Pvalue)),alpha=alpha)+theme_minimal()+
      theme(axis.text.y=element_text(face=font.y,size=fontsize.y),axis.text.x=element_text(face=font.x,color="black",size=fontsize.x))+
      scale_color_gradient(low=low,high=high)+ylab("Pathway name")+
      xlab("")+labs(size="Gene number")+guides(color=guide_colourbar(order = 1),size=guide_legend(order = 2))
  }
  if(!is.null(filename)){
    ggsave(p,file=paste(filename,"KEGG.pdf",sep="_"),width=width,height=height)
  }
  p
}
