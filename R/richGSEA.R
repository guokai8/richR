#' Enrichment analysis for any type of annotation data
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off method
#' @param KEGG a logical evaluating to TRUE or FALSE indicating whether KEGG GSEA were peformed or not.
#' @param padj.method p value adjust method
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param minGSSize minimal size of genes annotated by ontology term for testing.
#' @param maxGSSize maximal size of each geneset for analyzing
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
#' @importFrom fgsea fgseaMultilevel
#' @export
#' @author Kai Guo
richGSEA_internal<-function(x,object,keytype="",pvalue=0.05,padj=NULL,KEGG=FALSE,minSize=15,maxSize=500,
                            minGSSize = 10, maxGSSize = 500,
                            padj.method="BH",organism=NULL,ontology=NULL,table=TRUE,sep=","){
  x<-sort(x)
  if("Annot"%in%colnames(object)){
    object[,2]<-object$Annot
  }
  annod<-sf(object);
  res<-fgseaMultilevel(pathways=annod,stats=x,minSize=minSize,maxSize=maxSize)
  if(isTRUE(table)){
    res$leadingEdge<-unlist(lapply(res$leadingEdge, function(x)paste(gsub(" ","",x),collapse = sep,sep="")))
  }
  if(is.null(padj)){
    res<-res[res$pval<pvalue,]
    padj=numeric()
  }else{
    res<-res[res$padj<padj,]
  }
  res<-res[order(res$pval),]
  res <- res[!is.na(res$pathway),]
  if(isTRUE(KEGG)){
    data("path")
    rownames(path)<-path$Level3
    res<-cbind(res,path[res$pathway,])
  }
  if(is.null(organism)){
    organism=character()
  }
  if(is.null(ontology)){
    ontology=character()
  }
  result<-new("GSEAResult",
              result=res,
              pvalueCutoff   = pvalue,
              pAdjustMethod  = padj.method,
              padjCutoff   = padj,
              genenumber    = length(x),
              organism       = organism,
              gene           = names(x),
              input = x,
              ontology = ontology,
              keytype        = keytype,
              sep = sep
  )
  return(result)
}

#' GSEA Enrichment analysis function
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off method
#' @param KEGG a logical evaluating to TRUE or FALSE indicating whether KEGG GSEA were peformed or not.
#' @param organism organism name
#' @param keytype keytype for input genes
#' @param padj.method pvalue adjust method(default:"BH")
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
#' @export
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsako<-as.data.frame(hsako)
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' }
#' @author Kai Guo
setMethod("richGSEA", signature(object = "data.frame"),definition = function(x,object,keytype="",pvalue=0.05,padj=NULL,KEGG=FALSE,minSize=15,ontology=ontology,
                                                                             maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=",") {
  richGSEA_internal(x,object,keytype=keytype,pvalue=pvalue,padj=padj,KEGG=KEGG,minSize=minSize,ontology=ontology,
                    maxSize=maxSize,padj.method=padj.method,organism=organism,table=table,sep=sep)
})
#' GSEA Enrichment analysis function
#' @param x a vector include all log2FC with gene name
#' @param object annotation file for all genes
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off method
#' @param KEGG a logical evaluating to TRUE or FALSE indicating whether KEGG GSEA were peformed or not.
#' @param padj.method p value adjust method
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' }
#' @export
#' @author Kai Guo
setMethod("richGSEA", signature(object = "Annot"),definition = function(x,object,keytype="",pvalue=0.05,padj=NULL,KEGG=FALSE,minSize=15,ontology=ontology,
                                                                        maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=",") {
  richGSEA_internal(x,object@annot,keytype=object@keytype,pvalue=pvalue,padj=padj,KEGG=KEGG,minSize=minSize,ontology=object@anntype,
                    maxSize=maxSize,padj.method=padj.method,organism=object@species,table=table,sep=sep)
})


#' @name ggGSEA
#' @title plot the gsea result
#' @param object Annot object
#' @param term the significant term
#' @param gseaRes GSEAResult object
#' @param top number of top pathways to display
#' @param default whether to use default table plot or individual enrichment plots
#' @importFrom fgsea plotEnrichment plotGseaTable
#' @importFrom ggplot2 ggtitle
#' @importFrom cowplot plot_grid
#' @examples
#' \dontrun{
#' set.seed(123)
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' ggGSEA(term = res$pathway, object = hsako, gseaRes = res, default = FALSE)
#' }
#' @export
#' @author Kai Guo
ggGSEA <- function(term, object, gseaRes, top = 10, default = TRUE) {
  # Input validation
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  # Extract fold change vector from GSEAResult object
  x <- gseaRes@input
  x <- sort(x)

  annot <- object@annot
  gseaRes_result <- gseaRes@result

  if(nrow(gseaRes_result) == 0) {
    stop("No pathways found in GSEA results")
  }

  if(nrow(gseaRes_result) < top) {
    top <- nrow(gseaRes_result)
  }

  if(!is.null(annot$Annot)) {
    annot[,2] <- annot$Annot
  }

  annod <- sf(annot)

  if(length(term) > 1 & !is.null(gseaRes_result)) {
    if(isTRUE(default)) {
      plotGseaTable(annod[term], stats = x, gseaRes_result, gseaParam = 0.5)
    } else {
      res <- lapply(gseaRes_result$pathway[1:top], function(y) {
        plotEnrichment(annod[[y]], stats = x) + ggtitle(y)
      })
      plot_grid(plotlist = res, ncol = 5)
    }
  } else {
    plotEnrichment(annod[[term]], stats = x)
  }
}

#' plot multiple significant pathways with enhanced customization and directional coloring
#' @importFrom ggplot2 ggplot geom_hline aes geom_point geom_segment geom_line theme_bw element_blank theme scale_color_manual ggtitle element_text annotate ggsave guides guide_legend
#' @param object Annot object
#' @param gseaRes GSEAResult object
#' @param mycol a vector indicate the colors used for the figure (used when show_direction = FALSE)
#' @param up_color color for up-regulated pathways (positive NES)
#' @param down_color color for down-regulated pathways (negative NES)
#' @param show_direction logical, whether to use directional coloring based on NES
#' @param pathways specific pathway names to display (if NULL, use top/pvalue/padj filtering)
#' @param pathway_pattern regex pattern for pathway matching
#' @param category pathway category filter
#' @param min_genes minimum gene set size filter
#' @param max_genes maximum gene set size filter
#' @param top number of terms you want to display (used when pathways is NULL)
#' @param pvalue cutoff value of pvalue (if padj set as NULL, used when pathways is NULL)
#' @param padj cutoff value of p adjust value (used when pathways is NULL)
#' @param gseaParam GSEA parameter for enrichment score calculation
#' @param ticksSize size of the tick marks
#' @param plot_title title for the plot
#' @param show_pathway_names whether to show pathway names on plot
#' @param pathway_name_size size of pathway name text
#' @param line_alpha transparency of lines
#' @param point_alpha transparency of points
#' @param legend_position position of legend ("top", "bottom", "left", "right", "none")
#' @param legend_ncol number of columns in legend
#' @param return_data whether to return plot data along with plot
#' @param interactive whether to create interactive plot using plotly
#' @examples
#' \dontrun{
#' set.seed(123)
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' res<-richGSEA(gene,object = hsako)
#' # Plot with directional coloring (default)
#' plotGSEA(object = hsako, gseaRes = res, show_direction = TRUE)
#' # Plot specific pathways
#' plotGSEA(object = hsako, gseaRes = res,
#'          pathways = c("MAPK signaling pathway", "PI3K-Akt signaling pathway"))
#' # Plot pathways matching pattern
#' plotGSEA(object = hsako, gseaRes = res, pathway_pattern = "signaling")
#' # Traditional multi-color approach
#' plotGSEA(object = hsako, gseaRes = res, show_direction = FALSE)
#' }
#' @export
#' @author Kai Guo
plotGSEA <- function(object, gseaRes, mycol = NULL,
                     up_color = "#E31A1C", down_color = "#1F78B4", show_direction = TRUE,
                     pathways = NULL, pathway_pattern = NULL, category = NULL,
                     min_genes = NULL, max_genes = NULL, top = 10,
                     pvalue = 0.05, padj = NULL, gseaParam = 1, ticksSize = 0.2,
                     plot_title = "GSEA Enrichment Plot", show_pathway_names = FALSE,
                     pathway_name_size = 10, line_alpha = 0.8, point_alpha = 0.6,
                     legend_position = "bottom", legend_ncol = 2,
                     return_data = FALSE, interactive = FALSE) {

  # Input validation
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  # Extract fold change vector from GSEAResult object
  x <- gseaRes@input

  # Extract results from GSEAResult object
  gseaRes_result <- gseaRes@result

  if(nrow(gseaRes_result) == 0) {
    stop("No pathways found in GSEA results")
  }

  if(!is.null(pathways) && !is.character(pathways)) {
    stop("pathways must be a character vector")
  }

  # Select pathways to display
  if(!is.null(pathways)) {
    # Use user-specified pathways
    sigpathway <- pathways[pathways %in% gseaRes_result$pathway]
    if(length(sigpathway) == 0) {
      stop("None of the specified pathways found in GSEA results")
    }
    # Get corresponding NES values
    pathway_indices <- match(sigpathway, gseaRes_result$pathway)
    sig_nes <- gseaRes_result$NES[pathway_indices]
  } else if(!is.null(pathway_pattern)) {
    # Use pattern-based selection
    pattern_indices <- grepl(pathway_pattern, gseaRes_result$pathway, ignore.case = TRUE)
    sigpathway <- gseaRes_result$pathway[pattern_indices]
    sig_nes <- gseaRes_result$NES[pattern_indices]
    if(length(sigpathway) == 0) {
      stop("No pathways match the specified pattern")
    }
  } else {
    # Use filtering criteria
    if(!is.null(padj)) {
      cutoff <- padj
      sig_indices <- gseaRes_result$padj < cutoff
    } else {
      cutoff <- pvalue
      sig_indices <- gseaRes_result$pval < cutoff
    }

    sigpathway <- gseaRes_result$pathway[sig_indices]
    sig_nes <- gseaRes_result$NES[sig_indices]

    if(top > length(sigpathway)) {
      top <- length(sigpathway)
    }
    sigpathway <- sigpathway[1:top]
    sig_nes <- sig_nes[1:top]
  }

  # Apply gene set size filtering
  if(!is.null(min_genes) || !is.null(max_genes)) {
    pathway_indices <- match(sigpathway, gseaRes_result$pathway)
    valid_mask <- rep(TRUE, length(pathway_indices))

    if(!is.null(min_genes)) {
      valid_mask <- valid_mask & (gseaRes_result$size[pathway_indices] >= min_genes)
    }
    if(!is.null(max_genes)) {
      valid_mask <- valid_mask & (gseaRes_result$size[pathway_indices] <= max_genes)
    }

    sigpathway <- sigpathway[valid_mask]
    sig_nes <- sig_nes[valid_mask]
  }

  if(length(sigpathway) == 0) {
    stop("No pathways remaining after filtering")
  }

  # Set up colors based on direction preference
  if(show_direction) {
    # Create direction-based colors for individual pathways
    pathway_direction <- ifelse(sig_nes > 0, "Up-regulated", "Down-regulated")

    # Generate different shades/variants for each pathway within direction
    up_indices <- which(sig_nes > 0)
    down_indices <- which(sig_nes < 0)

    # Create color variations for up-regulated pathways
    if(length(up_indices) > 0) {
      up_colors <- colorRampPalette(c(up_color,
                                      adjustcolor(up_color, alpha.f = 0.7)))(length(up_indices))
    } else {
      up_colors <- character(0)
    }

    # Create color variations for down-regulated pathways
    if(length(down_indices) > 0) {
      down_colors <- colorRampPalette(c(down_color,
                                        adjustcolor(down_color, alpha.f = 0.7)))(length(down_indices))
    } else {
      down_colors <- character(0)
    }

    # Assign colors to pathways
    pathway_colors <- character(length(sigpathway))
    pathway_colors[up_indices] <- up_colors
    pathway_colors[down_indices] <- down_colors
    names(pathway_colors) <- sigpathway

    legend_title <- "Pathways (by Direction)"

  } else {
    # Use traditional multi-color approach
    if(is.null(mycol)) {
      mycol <- c("darkgreen", "chocolate4", "blueviolet", "#223D6C", "#D20A13", "#088247", "#58CDD9",
                 "#7A142C", "#5D90BA", "#431A3D", "#91612D", "#6E568C", "#E0367A", "#D8D155", "#64495D",
                 "#7CC767")
    }
    pathway_colors <- rep(mycol, length.out = length(sigpathway))
    names(pathway_colors) <- sigpathway
    legend_title <- "Pathway"
    pathway_direction <- rep("Individual", length(sigpathway))
  }

  # Calculate GSEA curves
  fc <- x
  res <- lapply(sigpathway, function(pathway) {
    .calGSEA(object, pathway, fc, gseaParam = gseaParam, ticksSize = ticksSize)
  })

  toPlot <- do.call(rbind, lapply(res, '[[', 'toPlot'))
  pathway_segments <- do.call(rbind, lapply(res, '[[', 'pathway'))
  tops <- unlist(lapply(res, '[[', 'tops'))
  bottoms <- unlist(lapply(res, '[[', 'bottoms'))

  # Keep pathways separate - each pathway keeps its individual identity
  # Add direction information for potential grouping in legend
  if(show_direction) {
    direction_mapping <- data.frame(
      Group = sigpathway,
      Direction = pathway_direction,
      stringsAsFactors = FALSE
    )
    toPlot <- merge(toPlot, direction_mapping, by = "Group", all.x = TRUE)
    pathway_segments <- merge(pathway_segments, direction_mapping, by = "Group", all.x = TRUE)
  }

  diff <- (max(tops) - min(bottoms))/8

  # Create the plot - each pathway maintains its individual line
  p <- ggplot(toPlot, aes(x = x, y = y, color = Group)) +
    geom_point(size = 0.1, alpha = point_alpha) +
    geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "black") +
    geom_line(alpha = line_alpha) +
    theme_bw() +
    geom_segment(data = pathway_segments,
                 mapping = aes(x = x, y = -diff/4, xend = x, yend = diff/4, color = Group),
                 size = ticksSize) +
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = legend_position,
          plot.title = element_text(hjust = 0.5, size = 14)) +
    labs(x = "Rank", y = "Enrichment Score") +
    ggtitle(plot_title)

  # Customize legend based on approach
  if(legend_position != "none") {
    if(show_direction) {
      # Create pathway labels with NES values and group by direction
      pathway_indices <- match(sigpathway, gseaRes_result$pathway)
      nes_values <- round(gseaRes_result$NES[pathway_indices], 2)

      # Create labels with pathway names and NES values
      pathway_labels_with_nes <- paste0(sigpathway, " (", nes_values, ")")

      # Truncate long pathway names but keep NES values
      max_pathway_length <- 40  # Adjust as needed
      truncated_pathways <- ifelse(nchar(sigpathway) > max_pathway_length,
                                   paste0(substr(sigpathway, 1, max_pathway_length-3), "..."),
                                   sigpathway)
      pathway_labels_with_nes <- paste0(truncated_pathways, " (", nes_values, ")")

      # Group pathways by direction
      up_indices <- which(sig_nes > 0)
      down_indices <- which(sig_nes < 0)

      # Create ordered labels (up-regulated first, then down-regulated)
      ordered_pathways <- character(0)
      ordered_colors <- character(0)
      ordered_labels <- character(0)

      if(length(up_indices) > 0) {
        # Sort up-regulated by NES (highest first)
        up_order <- up_indices[order(sig_nes[up_indices], decreasing = TRUE)]
        ordered_pathways <- c(ordered_pathways, sigpathway[up_order])
        ordered_colors <- c(ordered_colors, pathway_colors[sigpathway[up_order]])
        ordered_labels <- c(ordered_labels, pathway_labels_with_nes[up_order])
      }

      if(length(down_indices) > 0) {
        # Sort down-regulated by NES (least negative first, i.e., highest to lowest)
        down_order <- down_indices[order(sig_nes[down_indices], decreasing = TRUE)]
        ordered_pathways <- c(ordered_pathways, sigpathway[down_order])
        ordered_colors <- c(ordered_colors, pathway_colors[sigpathway[down_order]])
        ordered_labels <- c(ordered_labels, pathway_labels_with_nes[down_order])
      }

      # Add group separators in legend
      if(length(up_indices) > 0 && length(down_indices) > 0) {
        legend_title_with_groups <- paste0(legend_title, "\n",
                                           "↑ Up-regulated (", length(up_indices), ") | ",
                                           "↓ Down-regulated (", length(down_indices), ")")
      } else if(length(up_indices) > 0) {
        legend_title_with_groups <- paste0(legend_title, "\n↑ Up-regulated (", length(up_indices), ")")
      } else {
        legend_title_with_groups <- paste0(legend_title, "\n↓ Down-regulated (", length(down_indices), ")")
      }

      # Apply the colors and labels - single scale_color_manual call
      p <- p + scale_color_manual(values = setNames(ordered_colors, ordered_pathways),
                                  breaks = ordered_pathways,
                                  labels = ordered_labels,
                                  name = legend_title_with_groups) +
        guides(color = guide_legend(
          override.aes = list(size = 2, alpha = 1, linewidth = 1),
          ncol = legend_ncol,
          title.theme = element_text(size = 10, hjust = 0),
          label.theme = element_text(size = 8),
          keywidth = unit(1, "cm"),
          keyheight = unit(0.5, "cm")
        ))

    } else {
      # Traditional approach - single scale_color_manual call
      legend_labels <- names(pathway_colors)
      if(max(nchar(legend_labels)) > 50) {
        legend_labels <- ifelse(nchar(legend_labels) > 50,
                                paste0(substr(legend_labels, 1, 47), "..."),
                                legend_labels)
        names(pathway_colors) <- legend_labels
      }

      p <- p + scale_color_manual(values = pathway_colors, name = legend_title) +
        guides(color = guide_legend(override.aes = list(size = 2, alpha = 1),
                                    ncol = legend_ncol,
                                    title = legend_title))
    }
  } else {
    # No legend case - still need to set colors
    p <- p + scale_color_manual(values = pathway_colors, guide = "none")
  }

  # Add pathway names as annotations if requested
  if(show_pathway_names && length(sigpathway) <= 10) {
    annotation_y <- seq(max(tops) * 0.8, min(bottoms) * 0.8, length.out = length(sigpathway))
    p <- p + annotate("text", x = max(toPlot$x) * 0.7, y = annotation_y,
                      label = sigpathway, size = pathway_name_size/3, hjust = 0, vjust = 0.5)
  }

  # Create interactive plot if requested
  if(interactive) {
    if(!requireNamespace("plotly", quietly = TRUE)) {
      warning("plotly package not available, returning static plot")
    } else {
      p <- plotly::ggplotly(p)
    }
  }

  # Return data if requested
  if(return_data) {
    return(list(
      plot = p,
      data = toPlot,
      pathways = sigpathway,
      directions = if(show_direction) pathway_direction else NULL,
      nes_values = sig_nes,
      pathway_colors = pathway_colors,
      stats = list(tops = tops, bottoms = bottoms),
      pathway_data = pathway_segments
    ))
  } else {
    return(p)
  }
}

#' Get available pathways from GSEA results
#' @param gseaRes GSEAResult object
#' @param padj_cutoff padj cutoff for filtering (if NULL, no filtering applied)
#' @param pval_cutoff pvalue cutoff for filtering (used if padj_cutoff is NULL)
#' @param min_size minimum pathway size
#' @param max_size maximum pathway size
#' @export
#' @examples
#' \dontrun{
#' pathways <- getPathways(gseaRes, padj_cutoff = 0.05)
#' }
#' @author Kai Guo
getPathways <- function(gseaRes, padj_cutoff = 0.05, pval_cutoff = 0.05, min_size = NULL, max_size = NULL) {
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  result <- gseaRes@result

  # Apply significance filtering
  if(!is.null(padj_cutoff)) {
    result <- result[result$padj < padj_cutoff, ]
  } else if(!is.null(pval_cutoff)) {
    result <- result[result$pval < pval_cutoff, ]
  }

  # Apply size filtering
  if(!is.null(min_size)) {
    result <- result[result$size >= min_size, ]
  }
  if(!is.null(max_size)) {
    result <- result[result$size <= max_size, ]
  }

  return(result$pathway)
}

#' Search pathways by keyword
#' @param gseaRes GSEAResult object
#' @param keyword search term (can be regex pattern)
#' @param ignore_case whether to ignore case when searching
#' @export
#' @examples
#' \dontrun{
#' signaling_pathways <- searchPathways(gseaRes, "signaling")
#' mapk_pathways <- searchPathways(gseaRes, "MAPK")
#' }
#' @author Kai Guo
searchPathways <- function(gseaRes, keyword, ignore_case = TRUE) {
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  pathways <- gseaRes@result$pathway
  matched <- pathways[grepl(keyword, pathways, ignore.case = ignore_case)]

  if(length(matched) == 0) {
    message("No pathways found matching the keyword: ", keyword)
  }

  return(matched)
}

#' Save GSEA plot with optimal dimensions
#' @param plot_obj ggplot object from plotGSEA
#' @param filename output filename
#' @param width plot width in inches
#' @param height plot height in inches
#' @param dpi resolution in dots per inch
#' @param device output device (auto-detected from filename if NULL)
#' @importFrom ggplot2 ggsave
#' @export
#' @examples
#' \dontrun{
#' p <- plotGSEA(object = hsako, gseaRes = res)
#' saveGSEAplot(p, "gsea_plot.pdf")
#' saveGSEAplot(p, "gsea_plot.png", dpi = 600)
#' }
#' @author Kai Guo
saveGSEAplot <- function(plot_obj, filename, width = 12, height = 8, dpi = 300, device = NULL) {
  if(!inherits(plot_obj, c("gg", "ggplot"))) {
    stop("plot_obj must be a ggplot object")
  }

  ggsave(filename, plot_obj, width = width, height = height, dpi = dpi, device = device)
  message("Plot saved to: ", filename)
}

#' Create multiple GSEA plots for different comparisons
#' @param gsea_list named list of GSEAResult objects
#' @param object Annot object
#' @param output_dir directory to save plots
#' @param width plot width
#' @param height plot height
#' @param dpi resolution
#' @param format output format (pdf, png, etc.)
#' @param ... additional arguments passed to plotGSEA
#' @export
#' @examples
#' \dontrun{
#' gsea_results <- list(
#'   "D7_WT_vs_TFAM" = res1,
#'   "D14_WT_vs_D7_WT" = res2
#' )
#' plots <- batchGSEAplot(gsea_results, object = hsako)
#' }
#' @author Kai Guo
batchGSEAplot <- function(gsea_list, object, output_dir = "gsea_plots",
                          width = 12, height = 8, dpi = 300, format = "pdf", ...) {
  if(!is.list(gsea_list) || is.null(names(gsea_list))) {
    stop("gsea_list must be a named list")
  }

  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    message("Created output directory: ", output_dir)
  }

  plots <- lapply(names(gsea_list), function(name) {
    message("Processing: ", name)

    tryCatch({
      p <- plotGSEA(object = object, gseaRes = gsea_list[[name]],
                    plot_title = paste("GSEA:", name), ...)
      filename <- file.path(output_dir, paste0(name, "_gsea.", format))
      ggsave(filename, p, width = width, height = height, dpi = dpi)
      return(p)
    }, error = function(e) {
      warning("Failed to create plot for ", name, ": ", e$message)
      return(NULL)
    })
  })

  names(plots) <- names(gsea_list)
  plots <- plots[!sapply(plots, is.null)]  # Remove failed plots

  message("Successfully created ", length(plots), " plots in ", output_dir)
  return(plots)
}

#' Get pathway statistics summary
#' @param gseaRes GSEAResult object
#' @param top number of top pathways to summarize
#' @export
#' @examples
#' \dontrun{
#' summary <- getPathwayStats(gseaRes, top = 20)
#' print(summary)
#' }
#' @author Kai Guo
getPathwayStats <- function(gseaRes, top = 10) {
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  result <- gseaRes@result

  if(nrow(result) == 0) {
    return(data.frame())
  }

  # Get top pathways
  top_pathways <- head(result, top)

  # Create summary
  summary_stats <- data.frame(
    pathway = top_pathways$pathway,
    NES = round(top_pathways$NES, 3),
    pvalue = formatC(top_pathways$pval, format = "e", digits = 2),
    padj = formatC(top_pathways$padj, format = "e", digits = 2),
    size = top_pathways$size,
    direction = ifelse(top_pathways$NES > 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )

  return(summary_stats)
}

#' Filter pathways by multiple criteria
#' @param gseaRes GSEAResult object
#' @param nes_cutoff absolute NES cutoff
#' @param padj_cutoff padj cutoff
#' @param pval_cutoff pval cutoff (used if padj_cutoff is NULL)
#' @param min_size minimum pathway size
#' @param max_size maximum pathway size
#' @param direction filter by direction ("up", "down", or "both")
#' @export
#' @examples
#' \dontrun{
#' filtered <- filterPathways(gseaRes, nes_cutoff = 1.5, padj_cutoff = 0.01, direction = "up")
#' }
#' @author Kai Guo
filterPathways <- function(gseaRes, nes_cutoff = NULL, padj_cutoff = 0.05, pval_cutoff = 0.05,
                           min_size = NULL, max_size = NULL, direction = "both") {
  if(!inherits(gseaRes, "GSEAResult")) {
    stop("gseaRes must be a GSEAResult object")
  }

  result <- gseaRes@result

  # Apply NES filtering
  if(!is.null(nes_cutoff)) {
    result <- result[abs(result$NES) >= nes_cutoff, ]
  }

  # Apply significance filtering
  if(!is.null(padj_cutoff)) {
    result <- result[result$padj < padj_cutoff, ]
  } else if(!is.null(pval_cutoff)) {
    result <- result[result$pval < pval_cutoff, ]
  }

  # Apply size filtering
  if(!is.null(min_size)) {
    result <- result[result$size >= min_size, ]
  }
  if(!is.null(max_size)) {
    result <- result[result$size <= max_size, ]
  }

  # Apply direction filtering
  if(direction == "up") {
    result <- result[result$NES > 0, ]
  } else if(direction == "down") {
    result <- result[result$NES < 0, ]
  }

  return(result)
}

#' GSEA Enrichment analysis function for data frame with log2FC and name specifc
#' @param x data.frame include the gene name and log2FC
#' @param log2FC the column name for log2FoldChange
#' @param gene.col the gene name column if rownames won't be used as gene name
#' @param object annotation file for all genes
#' @param simplify boolean value to indicate if the simple table shoule be returned
#' @param pvalue pvalue cutoff value
#' @param padj adjust p value cut off method
#' @param padj.method p value adjust method
#' @param minSize Minimal size of a gene set to test. All pathways below the threshold are excluded.
#' @param maxSize Maximal size of a gene set to test. All pathways above the threshold are excluded.
#' @param table leadingEdge as vector
#' @param sep character string used to separate the genes when concatenating
#' @export
#' @examples
#' \dontrun{
#' hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#' hsako<-as.data.frame(hsako)
#' name=sample(unique(hsako$GeneID),1000)
#' gene<-rnorm(1000)
#' names(gene)<-name
#' d<-data.frame(Gene=name,log2FoldChange=gene)
#' rownames(d)<-d$Gene
#' res<-parGSEA(d,object = hsako)
#' }
#' @author Kai Guo
parGSEA<-function(x,object,log2FC="log2FoldChange",gene.col=NULL,
                  simplify=FALSE,keytype="",
                  pvalue=0.05,padj=NULL,KEGG=FALSE,minSize=15,ontology="",
                  maxSize=500,padj.method="BH",organism=NULL,table=TRUE,sep=","){
  fc<-x[,log2FC]
  if(is.null(gene.col)){
    names(fc)<-rownames(x)
  }else{
    names(fc)<-x[,gene.col]
  }
  tmp<-richGSEA(fc,object,keytype=keytype,pvalue=pvalue,padj=padj,KEGG=KEGG,minSize=minSize,ontology=ontology,
                maxSize=maxSize,padj.method=padj.method,organism=organism,table=table,sep=sep)
  if(isTRUE(simplify)){
    tmp<-result(tmp)
  }
  return(tmp)
}
