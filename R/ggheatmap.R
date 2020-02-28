#' Plot compare heatmap of Enrichment result among DEG groups
#' @importFrom dplyr full_join
#' @importFrom dplyr filter_
#' @importFrom dplyr arrange
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 geom_tile
#' @importFrom ggplot2 scale_fill_gradient2
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 theme_minimal
#' @importFrom ggplot2 coord_equal
#' @importFrom ggplot2 coord_flip
#' @importFrom reshape2 melt
#' @importFrom magrittr %>%
#' @param richRes list of enrichment object
#' @param top the number of Terms you want to display
#' @param colnames the compare DEG group names
#' @param xsize size of group name
#' @param ysize size of Terms name
#' @author Kai Guo
ggheatmap<-function(richRes, top = 50, colnames = NULL, xsize = 6, ysize = 6,usePadj=FALSE,
                     horizontal=FALSE,returnData=FALSE,...)
{
  object<-Reduce(function(x, y) rbind(x, y), lapply(richRes, function(x)x@result[1:top,]))
  object<-as.data.frame(na.omit(object))
  sel<-as.vector(unique(object$Term))
  if(isTRUE(usePadj)){
    res <- Reduce(function(x, y) full_join(x, y, by = "Term"),lapply(richRes,function(x)x@result[,c("Term","Padj")]))
  }else{
    res <- Reduce(function(x, y) full_join(x, y, by = "Term"),lapply(richRes,function(x)x@result[,c("Term","Pvalue")]))
  }
  if (!is.null(colnames)) {
    colnames(res)[2:ncol(res)] <- colnames
  }
  else {
    colnames(res)[2:ncol(res)] <- paste("Group", 1:(ncol(res)-1), sep = "_")
  }
  res[is.na(res)] <- 1
  res<-res%>%filter_(~Term%in%sel)
  res<-as.data.frame(res)
  rownames(res)<-res$Term
  cor_mat<-cor(t(res[,2:ncol(res)]))
  dd <- as.dist((1-cor_mat)/2);
  hc <- hclust(dd);
  melted <- melt(res[hc$order,])
  melted$Term<-factor(melted$Term,levels=res$Term[hc$order])
  maxp = max(-log10(melted[, 3])) + 0.5
  if(!isTRUE(usePadj)){
    colnames(melted)[3] <- "Padj"
    p<-ggplot(melted, aes(x = variable, y = Term, fill = -log10(Padj))) +coord_equal(ratio = 0.8)+
      geom_tile(color = "white") + scale_fill_gradient2(low = "white",high = "red", midpoint = 0, limit = c(0, maxp)) +
      theme_minimal()+theme(axis.text.y = element_text(size = ysize), axis.text.x = element_text(angle = 70,
      vjust = 1, size = xsize, hjust = 1,face = "bold")) +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }else{
    colnames(melted)[3] <- "Pvalue"
    p<-ggplot(melted, aes(x = variable, y = Term, fill = -log10(Pvalue))) +coord_equal(ratio = 0.8)+
      geom_tile(color = "white") + scale_fill_gradient2(low = "white",high = "red", midpoint = 0, limit = c(0, maxp)) +
      theme_minimal()+theme(axis.text.y = element_text(size = ysize), axis.text.x = element_text(angle = 70,
      vjust = 1, size = xsize, hjust = 1,face = "bold")) +
      theme(axis.title.x = element_blank(),axis.title.y = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }
  if(isTRUE(horizontal)){
    p<-p+coord_flip()
  }
  if(isTRUE(returnData)){
    return(res)
  }else{
    p
  }
}
