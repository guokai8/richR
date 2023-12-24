#' compare enrichment results across different samples
#' @param x list of richResults
#' @param pvalue cutoff pvalue
#' @param padj cutoff p adjust value
#' @param include.all include all richResults even empty
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
compareResult<-function (x, pvalue = 0.05, padj = NULL,include.all=FALSE)
{
  if (!is.null(padj)) {
    pvalue <- 1
  }
  else {
    padj <- 1
  }
  if (is.null(names(x))) {
    names(x) <- paste("Group", seq_along(x))
  }
  tmp <- lapply(x, function(x) filter(x, Padj < padj, Pvalue <
                                        pvalue))
  tmp <- lapply(names(x), function(x) mutate(tmp[[x]], group = x))
  dx <- do.call(rbind, tmp)
  if(isTRUE(include.all)){
    id<-setdiff(names(x),unique(dx$group))
    dd<-data.frame(Annot=dx$Annot[1],Term=dx$Term[1],
                   Annotated=min(dx$Annotated),Significant=min(dx$Significant),Pvalue=1,
                   Padj=1,GeneID="",group=id)
    dx<-rbind(dx,dd)
  }
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
##' @param short automatic short name or not
##' @param padj cutoff value of p adjust value
##' @param usePadj use p adjust value as color or not (should use with padj)
##' @param filename figure output name
##' @param width figure width
##' @param height figure height
##' @param include.all include all richResults even empty
##' @examples
#' \dontrun{
#' hsako <- buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsago <- buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
#' gene1 <- sample(unique(hsako$GeneID),1000)
#' gene2 <- sample(unique(hsako$GeneID),1000)
#' resko1 <-richKEGG(gene1,kodata = hsako)
#' resko2 <-richKEGG(gene2,kodata = hsako)
#' res<-compareResult(list(S1=resko1,S2=resko2))
#' comparedot(res,pvalue=0.05)
#' }
#' @author Kai Guo
#' @export
comparedot<-function (x, pvalue = 0.05, low = "lightpink", high = "red",
                      alpha = 0.7, font.x = "bold", font.y = "bold", fontsize.x = 10,
                      fontsize.y = 10, short = FALSE, padj = NULL, usePadj = TRUE,
                      filename = NULL, width = 10, height = 8,include.all=FALSE)
{
  if(isTRUE(include.all)){
    padj<-pvalue<-1.1
    low<-"white"
  }
  if (!is.null(padj)) {
    x <- x[x$Padj < padj, ]
  }
  else {
    x <- x[x$Pvalue < pvalue, ]
  }
  if (isTRUE(short)) {
    x$Term <- unlist(lapply(x$Term, function(x) .paste.char(x)))
  }
  if (isTRUE(usePadj)) {
    p <- ggplot(x, aes(x = group, y = Term)) + geom_point(aes(size = Significant,
                                                              color = -log10(Padj)), alpha = alpha) + theme_minimal() +
      theme(axis.text.y = element_text(face = font.y, size = fontsize.y),
            axis.text.x = element_text(face = font.x, color = "black",
                                       size = fontsize.x)) + scale_color_gradient(low = low,
                                                                                  high = high) + ylab("Pathway name") + xlab("") +
      labs(size = "Gene number") + guides(color = guide_colourbar(order = 1),
                                          size = guide_legend(order = 2))
  }
  else {
    p <- ggplot(x, aes(x = group, y = Term)) + geom_point(aes(size = Significant,
                                                              color = -log10(Pvalue)), alpha = alpha) + theme_minimal() +
      theme(axis.text.y = element_text(face = font.y, size = fontsize.y),
            axis.text.x = element_text(face = font.x, color = "black",
                                       size = fontsize.x)) + scale_color_gradient(low = low,
                                                                                  high = high) + ylab("Pathway name") + xlab("") +
      labs(size = "Gene number") + guides(color = guide_colourbar(order = 1),
                                          size = guide_legend(order = 2))
  }
  if (!is.null(filename)) {
    ggsave(p, file = paste(filename, "KEGG.pdf", sep = "_"),
           width = width, height = height)
  }
  p
}

#'
#' compare GSEA enrichment results across different samples and generate figure
#' @importFrom ggplot2 facet_wrap
#' @importFrom ggplot2 ggplot geom_hline aes geom_point geom_segment
#' @importFrom ggplot2 geom_line theme_bw element_blank theme scale_color_manual
#' @param x list of GSEAResult
#' @param object Annot object
#' @param gene list of a vector include all log2FC with gene name (optional)
#' @param pathway pathways you want to display (optional)
#' @param mycol a vector indicate the colors used for the figure
#' @param top number of terms you want to display,
#' @param pvalue cutoff value of pvalue (if padj set as NULL)
#' @param padj cutoff value of p adjust value
#' @param scales Should scales be fixed ("fixed", the default), free ("free"), or free in one dimension ("free_x", "free_y")?
#' @examples
#' \dontrun{
#' hsako <- buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' gene1 <- sample(unique(hsako$GeneID),1000)
#' gene2 <- sample(unique(hsako$GeneID),1000)
#' fc1<-rnorm(1000,11,2)
#' names(fc1)<-gene1
#' fc2<-rnorm(1000,11,2)
#' names(fc2)<-gene2
#' resko1 <-richGSEA(fc1,kodata = hsako)
#' resko2 <-richGSEA(fc2,kodata = hsako)
#' res<-compareGSEA(list(S1=resko1,S2=resko2),hsako)
#' }
#' @author Kai Guo
#' @export
compareGSEA<-function(x,object,gene=NULL,pathway=NULL,
                      mycol=NULL,pvalue = 0.05, padj = NULL,
                      gseaParam = 1, ticksSize = 0.2,ncol=2,scales="fixed"){
  if (!is.null(padj)) {
    Pvalue <- 1
  }
  else {
    Padj <- 1
  }
  Pvalue=pvalue
  if (is.null(names(x))) {
    names(x) <- paste("Group", seq_along(x))
  }
  if(is.null(mycol)){
    mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9",
               "#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D",
               "#7CC767")
  }
  tmp <- lapply(x, function(x) filter(x, padj < Padj, pval <
                                        Pvalue))
  sigpathway <- lapply(tmp, function(x) x$pathway)
 # names(tmp)<-names(x)
  if(is.null(gene)){
    fc <- lapply(x,function(x){
      fc<-x@input
      names(fc)<-x@gene
      return(fc)
    })}else{
    fc<-gene
    names(fc)<-names(tmp)
  }
#  names(fc)<-names(x)
  res<-list()
  for(i in names(tmp)){
    sigp<-sigpathway[[i]]
    fct<-fc[[i]]
    sigt<-lapply(sigp,function(x).calGSEA(object,x,fct,gseaParam=gseaParam,ticksSize=ticksSize))
    res[[i]]<-sigt
  }
  ####
  toPlot<-lapply(res,function(x)lapply(x,'[[','toPlot'))
  path<-lapply(res,function(x)lapply(x,'[[','pathway'))
  tops<-lapply(res,function(x)lapply(x,'[[','tops'))
  bottoms<-lapply(res,function(x)lapply(x,'[[','bottoms'))
  ######
  toPlot<-do.call(rbind,lapply(toPlot, function(x)do.call(rbind,x)))
  path<-do.call(rbind,lapply(path, function(x)do.call(rbind,x)))
 # tops<-do.call(rbind,lapply(tops, function(x)do.call(rbind,x)))
 # bottoms<-do.call(rbind,lapply(bottoms, function(x)do.call(rbind,x)))
  ####
  toPlot$group<-sub("(\\.)\\d+$", "", rownames(toPlot))
  path$group<-sub("(\\.)\\d+$", "", rownames(path))
  if(!is.null(pathway)){
    toPlot<-subset(toPlot,Group%in%pathway)
    path<-subset(path,Group%in%pathway)
  }
#  tops$group<-sub("(\\.)\\d+$", "", rownames(tops))
#  bottoms$group<-sub("(\\.)\\d+$", "", rownames(bottoms))
  tops<-max(toPlot$y)
  bottoms<-min(toPlot$y)
  ####
  diff <- (max(tops) - min(bottoms))/8
  p<-ggplot(toPlot, aes(x = x, y = y,color=group)) + geom_point(size = 0.1) +
    geom_hline(yintercept = tops, colour = "red",
               linetype = "dashed") + geom_hline(yintercept = bottoms,
                                                 colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "black") + geom_line() + theme_bw()
  p <- p+geom_segment(data =path, mapping = aes(x = x,
                                                   y = -diff/4, xend = x,
                                                   yend = diff/4,color=group),
                      size = ticksSize) +
    theme(panel.border = element_blank(), panel.grid.minor = element_blank()) +
    scale_color_manual(values=mycol)+labs(x = "rank", y = "Enrichment score")
  p<-p+facet_wrap(.~Group,ncol=ncol,scales = scales)
  p
}

