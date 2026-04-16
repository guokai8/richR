##' @importFrom ggplot2 ggplot aes geom_segment geom_point coord_flip labs
##' @importFrom ggplot2 theme_minimal geom_bar coord_polar scale_fill_gradient
##' @importFrom ggplot2 scale_fill_manual geom_text scale_color_gradient
##' @importFrom ggplot2 scale_color_gradient2 scale_color_manual geom_vline
##' @importFrom ggplot2 geom_tile scale_size_continuous geom_step
##' @importFrom ggplot2 element_text ggsave theme scale_fill_gradient2
NULL

# ============================================================
#  1. gglollipop — Lollipop chart of enrichment significance
# ============================================================

#' Lollipop plot of enrichment significance
#'
#' Draws a lollipop chart where each stem represents a pathway /
#' GO term ordered by significance.
#'
#' @param resultFis data.frame of enrichment results (columns: Term, Pvalue, Padj)
#' @param top number of top terms to display (default 20)
#' @param pvalue p-value cutoff when \code{padj} is NULL
#' @param padj adjusted-p cutoff (takes precedence when non-NULL)
#' @param usePadj logical; colour by adjusted p-value?
#' @param low low-end gradient colour (default \code{"lightpink"})
#' @param high high-end gradient colour (default \code{"red"})
#' @param point.size point size (default 4)
#' @param short shorten long term names (default FALSE)
#' @param filename save path (NULL = no save)
#' @param width figure width (inches)
#' @param height figure height (inches)
#' @return ggplot2 object (invisibly saved if \code{filename} is given)
#' @keywords internal
gglollipop_internal <- function(resultFis,
                                top      = 20,
                                pvalue   = 0.05,
                                padj     = NULL,
                                usePadj  = TRUE,
                                low      = "lightpink",
                                high     = "red",
                                point.size = 4,
                                short    = FALSE,
                                filename = NULL,
                                width    = 10,
                                height   = 8) {

  resultFis <- .filter_and_top(resultFis, top = top, pvalue = pvalue, padj = padj)
  if (nrow(resultFis) == 0) { message("No significant terms to plot."); return(invisible(NULL)) }
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  neg_log <- if (isTRUE(usePadj) && "Padj" %in% names(resultFis)) -log10(resultFis$Padj) else -log10(resultFis$Pvalue)
  resultFis$neg_log_p <- neg_log
  if (!"RichFactor" %in% names(resultFis) && "Significant" %in% names(resultFis) && "Annotated" %in% names(resultFis))
    resultFis$RichFactor <- resultFis$Significant / resultFis$Annotated

  color_label <- if (isTRUE(usePadj)) "-log10(Padj)" else "-log10(Pvalue)"

  p <- ggplot2::ggplot(resultFis, ggplot2::aes(x = stats::reorder(Term, RichFactor), y = RichFactor)) +
    ggplot2::geom_segment(ggplot2::aes(xend = Term, yend = 0), color = "grey50") +
    ggplot2::geom_point(ggplot2::aes(color = neg_log_p), size = point.size) +
    ggplot2::scale_color_gradient(low = low, high = high, name = color_label) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Lollipop Plot of Enrichment Significance", x = NULL,
                  y = "RichFactor") +
    ggplot2::theme_minimal()

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
#  2. ggvolcano — Volcano plot for GSEA results (NES vs -log10 padj)
# ============================================================

#' Volcano plot of GSEA gene sets
#'
#' Scatter NES vs significance, highlighting significant gene sets.
#'
#' @param resultFis data.frame with GSEA or ORA results (columns vary; NES/pval/padj for GSEA, FoldEnrichment/Pvalue/Padj for ORA)
#' @param top max number of top terms to display (default 50)
#' @param pvalue p-value cutoff
#' @param padj adjusted-p cutoff (overrides \code{pvalue} when supplied)
#' @param usePadj logical, use adjusted p-value column for significance (default TRUE)
#' @param low colour for down-regulated (negative-effect) points (default \code{"#2166ac"})
#' @param mid colour at zero effect (default \code{"white"})
#' @param high colour for up-regulated (positive-effect) points (default \code{"#b2182b"})
#' @param size.range numeric length-2, min/max point sizes mapped from -log10(p) (default \code{c(2, 8)})
#' @param label.size font size for point labels (default 3)
#' @param label.top number of top terms to label (default 5)
#' @param short logical, shorten term names via \code{.paste.char()} (default FALSE)
#' @param filename if non-NULL, save the plot to this file
#' @param width figure width
#' @param height figure height
#' @return a ggplot2 object
#' @keywords internal
ggvolcano_internal <- function(resultFis,
                               top        = 50,
                               pvalue     = 0.05,
                               padj       = NULL,
                               usePadj    = TRUE,
                               low        = "#2166ac",
                               mid        = "white",
                               high       = "#b2182b",
                               size.range = c(2, 8),
                               label.size = 3,
                               label.top  = 5,
                               short      = FALSE,
                               filename   = NULL,
                               width      = 10,
                               height     = 8) {

  resultFis <- .standardise_gsea_cols(resultFis)
  if (nrow(resultFis) > top) resultFis <- resultFis[seq_len(top), ]
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  # x-axis: log2(FoldEnrichment) for ORA, NES for GSEA
  is_gsea <- "NES" %in% names(resultFis)
  if (is_gsea) {
    resultFis$effect <- resultFis$NES
    x_label <- "Normalized Enrichment Score (NES)"
  } else {
    if (!"FoldEnrichment" %in% names(resultFis)) {
      if ("Significant" %in% names(resultFis) && "Annotated" %in% names(resultFis)) {
        resultFis$FoldEnrichment <- resultFis$Significant / resultFis$Annotated
      } else {
        resultFis$FoldEnrichment <- 1
      }
    }
    resultFis$effect <- log2(pmax(resultFis$FoldEnrichment, 0.001))
    x_label <- "log2(Fold Enrichment)"
  }

  # y / size: -log10(p)
  sig_col <- if (isTRUE(usePadj) && "Padj" %in% names(resultFis)) "Padj" else "Pvalue"
  resultFis$neg_log_p <- -log10(resultFis[[sig_col]])
  cutoff <- if (!is.null(padj)) padj else pvalue

  # Label top terms by significance
  resultFis <- resultFis[order(resultFis$neg_log_p, decreasing = TRUE), ]
  resultFis$label <- ""
  n_label <- min(label.top, nrow(resultFis))
  resultFis$label[seq_len(n_label)] <- resultFis$Term[seq_len(n_label)]

  p <- ggplot2::ggplot(resultFis, ggplot2::aes(x = effect, y = neg_log_p,
                                                color = effect, size = neg_log_p)) +
    ggplot2::geom_point(alpha = 0.8) +
    ggplot2::scale_color_gradient2(low = low, mid = mid, high = high, midpoint = 0,
                                   name = if (is_gsea) "NES" else "log2(FE)") +
    ggplot2::scale_size_continuous(range = size.range,
                                   name = paste0("-log10(", sig_col, ")")) +
    ggplot2::geom_hline(yintercept = -log10(cutoff), linetype = "dashed", color = "grey50") +
    ggplot2::labs(title = "Enrichment Volcano Plot",
                  x = x_label,
                  y = paste0("-log10(", sig_col, ")")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  if (!is_gsea) p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")

  if (any(resultFis$label != "")) {
    p <- p + ggrepel::geom_text_repel(
      data = resultFis[resultFis$label != "", ],
      ggplot2::aes(label = label),
      size = label.size, max.overlaps = 20, show.legend = FALSE,
      segment.color = "grey50", segment.alpha = 0.6)
  }

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
#  3. ggcircbar — Circular barplot of enrichment results
# ============================================================

#' Circular barplot (polar bar) of enrichment results
#'
#' @param resultFis data.frame of enrichment results
#' @param top number of top terms to display (default 15)
#' @param pvalue p-value cutoff
#' @param padj adjusted-p cutoff
#' @param usePadj colour by adjusted p-value?
#' @param low low-end gradient colour
#' @param high high-end gradient colour
#' @param short shorten long term names
#' @param filename if non-NULL, save figure
#' @param width figure width
#' @param height figure height
#' @return a ggplot2 object
#' @keywords internal
ggcircbar_internal <- function(resultFis,
                               top    = 15,
                               pvalue = 0.05,
                               padj   = NULL,
                               usePadj = TRUE,
                               low    = "#fee0d2",
                               high   = "#b2182b",
                               short  = FALSE,
                               filename = NULL,
                               width  = 10,
                               height = 8) {

  resultFis <- .filter_and_top(resultFis, top = top, pvalue = pvalue, padj = padj)
  if (nrow(resultFis) == 0) { message("No significant terms to plot."); return(invisible(NULL)) }
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  fill_val <- if (isTRUE(usePadj) && "Padj" %in% names(resultFis)) -log10(resultFis$Padj) else -log10(resultFis$Pvalue)
  resultFis$neg_log_p <- fill_val

  p <- ggplot2::ggplot(resultFis, ggplot2::aes(x = stats::reorder(Term, Significant),
                                                y = Significant, fill = neg_log_p)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_polar(start = 0) +
    ggplot2::scale_fill_gradient(low = low, high = high,
                                 name = if (isTRUE(usePadj)) "-log10(Padj)" else "-log10(Pvalue)") +
    ggplot2::labs(title = "Circular Barplot of Enrichment", x = NULL, y = NULL) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 0))

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
#  4. ggNES — Bidirectional NES barplot for GSEA results
# ============================================================

#' Bidirectional NES barplot of GSEA results
#'
#' @param resultFis data.frame with GSEA results (must contain NES column)
#' @param top number of top terms to display (default 20)
#' @param pvalue p-value cutoff
#' @param padj adjusted-p cutoff
#' @param up.color colour for positive NES (default \code{"#b2182b"})
#' @param down.color colour for negative NES (default \code{"#2166ac"})
#' @param short shorten long term names
#' @param filename if non-NULL, save figure
#' @param width figure width
#' @param height figure height
#' @return a ggplot2 object
#' @keywords internal
ggNES_internal <- function(resultFis,
                           top        = 20,
                           pvalue     = 0.05,
                           padj       = NULL,
                           up.color   = "#b2182b",
                           down.color = "#2166ac",
                           short      = FALSE,
                           filename   = NULL,
                           width      = 10,
                           height     = 8) {

  resultFis <- .standardise_gsea_cols(resultFis)
  if (!is.null(padj))  resultFis <- resultFis[resultFis$Padj < padj, ]
  else                 resultFis <- resultFis[resultFis$Pvalue < pvalue, ]

  resultFis <- resultFis[order(resultFis$NES, decreasing = TRUE), ]
  if (nrow(resultFis) > top) resultFis <- resultFis[seq_len(top), ]
  if (nrow(resultFis) == 0) { message("No significant terms to plot."); return(invisible(NULL)) }
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  resultFis$direction <- ifelse(resultFis$NES > 0, "Up", "Down")
  resultFis$direction <- factor(resultFis$direction, levels = c("Down", "Up"))
  # Drop unused levels so legend only shows directions present in data
  resultFis$direction <- droplevels(resultFis$direction)

  dir_colors <- c(Down = down.color, Up = up.color)
  dir_colors <- dir_colors[levels(resultFis$direction)]

  p <- ggplot2::ggplot(resultFis, ggplot2::aes(x = stats::reorder(Term, NES),
                                                y = NES, fill = direction)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_manual(values = dir_colors, name = "Direction") +
    ggplot2::labs(title = "GSEA NES Barplot", x = NULL, y = "Normalized Enrichment Score") +
    ggplot2::theme_minimal()

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
#  5. ggscatter — Labeled scatter: gene count vs significance
# ============================================================

#' Labeled scatter plot of enrichment results
#'
#' @param resultFis data.frame of enrichment results
#' @param top number of top terms to display (default 20)
#' @param pvalue p-value cutoff
#' @param padj adjusted-p cutoff
#' @param usePadj colour by adjusted p-value?
#' @param low low-end gradient colour
#' @param high high-end gradient colour
#' @param point.size point size (default 5)
#' @param label.size text label size (default 3)
#' @param short shorten long term names
#' @param filename if non-NULL, save figure
#' @param width figure width
#' @param height figure height
#' @return a ggplot2 object
#' @keywords internal
ggscatter_internal <- function(resultFis,
                               top        = 50,
                               pvalue     = 0.05,
                               padj       = NULL,
                               usePadj    = TRUE,
                               low        = "#fee0d2",
                               high       = "#b2182b",
                               point.size = 3,
                               label.size = 3,
                               label.top  = 5,
                               short      = FALSE,
                               filename   = NULL,
                               width      = 10,
                               height     = 8) {

  resultFis <- .filter_and_top(resultFis, top = top, pvalue = pvalue, padj = padj)
  if (nrow(resultFis) == 0) { message("No significant terms to plot."); return(invisible(NULL)) }
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  # x-axis: log10(Annotated) — pathway size
  resultFis$log_annotated <- log10(resultFis$Annotated)

  # y-axis: log2(FoldEnrichment) — effect size
  if (!"FoldEnrichment" %in% names(resultFis)) {
    if ("Significant" %in% names(resultFis) && "Annotated" %in% names(resultFis)) {
      resultFis$FoldEnrichment <- resultFis$Significant / resultFis$Annotated
    } else {
      resultFis$FoldEnrichment <- 1
    }
  }
  resultFis$log2_fe <- log2(pmax(resultFis$FoldEnrichment, 0.001))

  # color: -log10(p)
  sig_col <- if (isTRUE(usePadj) && "Padj" %in% names(resultFis)) "Padj" else "Pvalue"
  resultFis$neg_log_p <- -log10(resultFis[[sig_col]])
  color_label <- paste0("-log10(", sig_col, ")")

  # Label top terms by significance
  resultFis <- resultFis[order(resultFis$neg_log_p, decreasing = TRUE), ]
  resultFis$label <- ""
  n_label <- min(label.top, nrow(resultFis))
  resultFis$label[seq_len(n_label)] <- resultFis$Term[seq_len(n_label)]

  p <- ggplot2::ggplot(resultFis, ggplot2::aes(x = log_annotated, y = log2_fe,
                                                color = neg_log_p)) +
    ggplot2::geom_point(size = point.size, alpha = 0.8) +
    ggplot2::scale_color_gradient(low = low, high = high, name = color_label) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::labs(x = "log10(Pathway Size)",
                  y = "log2(Fold Enrichment)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  if (any(resultFis$label != "")) {
    p <- p + ggrepel::geom_text_repel(
      data = resultFis[resultFis$label != "", ],
      ggplot2::aes(label = label),
      size = label.size, max.overlaps = 20, show.legend = FALSE,
      segment.color = "grey50", segment.alpha = 0.6)
  }

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
#  6. ggecdf — ECDF step plot for GSEA ranked gene sets
# ============================================================

#' ECDF step plot of GSEA ranked gene sets
#'
#' @param resultFis data.frame with GSEA results (must contain NES column)
#' @param top number of top terms to display (default 30)
#' @param pvalue p-value cutoff
#' @param padj adjusted-p cutoff
#' @param line.color line colour (default \code{"#1b7837"})
#' @param line.size line size (default 1.2)
#' @param point.size point size (default 3)
#' @param short shorten long term names
#' @param filename if non-NULL, save figure
#' @param width figure width
#' @param height figure height
#' @return a ggplot2 object
#' @keywords internal
ggecdf_internal <- function(resultFis,
                            top        = 30,
                            pvalue     = 0.05,
                            padj       = NULL,
                            line.color = "#1b7837",
                            line.size  = 1.2,
                            point.size = 3,
                            short      = FALSE,
                            filename   = NULL,
                            width      = 10,
                            height     = 8) {

  resultFis <- .standardise_gsea_cols(resultFis)
  if (!is.null(padj))  resultFis <- resultFis[resultFis$Padj < padj, ]
  else                 resultFis <- resultFis[resultFis$Pvalue < pvalue, ]
  resultFis <- resultFis[order(resultFis$NES), ]
  if (nrow(resultFis) > top) resultFis <- resultFis[seq_len(top), ]
  if (nrow(resultFis) == 0) { message("No significant terms to plot."); return(invisible(NULL)) }
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  resultFis$cum_rank <- cumsum(resultFis$NES) / sum(abs(resultFis$NES))

  p <- ggplot2::ggplot(resultFis, ggplot2::aes(x = stats::reorder(Term, NES), y = cum_rank, group = 1)) +
    ggplot2::geom_step(linewidth = line.size, color = line.color) +
    ggplot2::geom_point(size = point.size, color = line.color) +
    ggplot2::labs(title = "ECDF of GSEA Ranked Gene Sets", x = "Gene Sets", y = "Cumulative Proportion") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
#  7. ggtermsim — Term similarity scatter (MDS-like)
# ============================================================

#' Term similarity scatter plot (MDS)
#'
#' Computes a simple MDS on the gene-overlap Jaccard matrix
#' between enriched terms and plots each term as a point.
#'
#' @param resultFis data.frame of enrichment results (must have Term, Pvalue/Padj, Significant, GeneID)
#' @param top number of top terms to include (default 20)
#' @param pvalue p-value cutoff
#' @param padj adjusted-p cutoff
#' @param usePadj colour by adjusted p-value?
#' @param low low-end gradient colour
#' @param high high-end gradient colour
#' @param label.size text label size
#' @param short shorten long term names
#' @param sep separator for gene IDs
#' @param filename if non-NULL, save figure
#' @param width figure width
#' @param height figure height
#' @return a ggplot2 object
#' @keywords internal
ggtermsim_internal <- function(resultFis,
                               top        = 20,
                               pvalue     = 0.05,
                               padj       = NULL,
                               usePadj    = TRUE,
                               low        = "grey80",
                               high       = "#b2182b",
                               label.size = 3,
                               sim.cutoff = 0.2,
                               size.range = c(3, 10),
                               short      = FALSE,
                               sep        = ",",
                               filename   = NULL,
                               width      = 10,
                               height     = 8) {

  resultFis <- .filter_and_top(resultFis, top = top, pvalue = pvalue, padj = padj)
  if (nrow(resultFis) < 3) { message("Need at least 3 terms for MDS."); return(invisible(NULL)) }
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  # Jaccard similarity matrix
  gene_lists <- strsplit(as.character(resultFis$GeneID), split = sep)
  n <- length(gene_lists)
  sim <- matrix(0, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      inter <- length(intersect(gene_lists[[i]], gene_lists[[j]]))
      uni   <- length(union(gene_lists[[i]], gene_lists[[j]]))
      sim[i, j] <- if (uni > 0) inter / uni else 0
    }
  }
  # MDS on distance (1 - similarity)
  dist_mat <- 1 - sim
  mds <- stats::cmdscale(stats::as.dist(dist_mat), k = 2)

  mds_df <- data.frame(
    Term        = resultFis$Term,
    Dim1        = mds[, 1],
    Dim2        = mds[, 2],
    Significant = resultFis$Significant,
    neg_log_p   = if (isTRUE(usePadj) && "Padj" %in% names(resultFis)) -log10(resultFis$Padj) else -log10(resultFis$Pvalue)
  )

  # Build edges for term pairs above similarity cutoff
  edges <- data.frame(x = numeric(0), y = numeric(0),
                      xend = numeric(0), yend = numeric(0),
                      similarity = numeric(0))
  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      if (sim[i, j] >= sim.cutoff) {
        edges <- rbind(edges, data.frame(
          x = mds[i, 1], y = mds[i, 2],
          xend = mds[j, 1], yend = mds[j, 2],
          similarity = sim[i, j]))
      }
    }
  }

  p <- ggplot2::ggplot()

  # Draw edges first (behind nodes)
  if (nrow(edges) > 0) {
    p <- p + ggplot2::geom_segment(data = edges,
             ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = similarity),
             color = "grey70", alpha = 0.5, show.legend = FALSE) +
      ggplot2::scale_linewidth_continuous(range = c(0.3, 2))
  }

  # Nodes
  p <- p +
    ggplot2::geom_point(data = mds_df,
             ggplot2::aes(x = Dim1, y = Dim2, size = Significant, color = neg_log_p),
             alpha = 0.85) +
    ggplot2::scale_color_gradient(low = low, high = high,
                                  name = if (isTRUE(usePadj)) "-log10(Padj)" else "-log10(Pvalue)") +
    ggplot2::scale_size_continuous(range = size.range, name = "Gene Count") +
    ggplot2::guides(color = ggplot2::guide_colourbar(order = 1),
                    size  = ggplot2::guide_legend(order = 2)) +
    ggrepel::geom_text_repel(data = mds_df,
              ggplot2::aes(x = Dim1, y = Dim2, label = Term),
              size = label.size, max.overlaps = 20, show.legend = FALSE,
              segment.color = "grey50", segment.alpha = 0.5) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text  = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank())

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
#  8. gggenedot — Gene–pathway association dotplot
# ============================================================

#' Gene–pathway association dotplot
#'
#' For each enriched term, shows contributing genes as dots
#' sized by fold-change magnitude and coloured by direction.
#'
#' @param resultFis data.frame of enrichment results (columns: Term, GeneID)
#' @param fc named numeric vector of fold-changes (names = gene IDs)
#' @param top number of top terms (default 10)
#' @param max.genes max genes per term to display (default 30)
#' @param pvalue p-value cutoff
#' @param padj adjusted-p cutoff
#' @param low colour for negative FC
#' @param mid colour for zero FC
#' @param high colour for positive FC
#' @param short shorten long term names
#' @param sep separator for gene IDs
#' @param filename if non-NULL, save figure
#' @param width figure width
#' @param height figure height
#' @return a ggplot2 object
#' @keywords internal
gggenedot_internal <- function(resultFis,
                               fc         = NULL,
                               top        = 10,
                               max.genes  = 30,
                               pvalue     = 0.05,
                               padj       = NULL,
                               usePadj    = TRUE,
                               low        = "#2166ac",
                               mid        = "white",
                               high       = "#b2182b",
                               short      = FALSE,
                               sep        = ",",
                               filename   = NULL,
                               width      = 12,
                               height     = 8) {

  resultFis <- .filter_and_top(resultFis, top = top, pvalue = pvalue, padj = padj)
  if (nrow(resultFis) == 0) { message("No significant terms to plot."); return(invisible(NULL)) }
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  melted <- .melt_gene_term(resultFis, fc = fc, sep = sep, max.genes = max.genes)
  if (nrow(melted) == 0) { message("No gene-term pairs."); return(invisible(NULL)) }

  if (is.null(fc)) {
    # No fold-change: colour by -log10(Pvalue/Padj), no size mapping
    neg_log <- if (isTRUE(usePadj) && "Padj" %in% names(resultFis)) -log10(resultFis$Padj) else -log10(resultFis$Pvalue)
    pval_map <- stats::setNames(neg_log, resultFis$Term)
    melted$neg_log_p <- pval_map[as.character(melted$Term)]
    color_label <- if (isTRUE(usePadj)) "-log10(Padj)" else "-log10(Pvalue)"
    p <- ggplot2::ggplot(melted, ggplot2::aes(x = Gene, y = Term, color = neg_log_p)) +
      ggplot2::geom_point(size = 4, alpha = 0.85) +
      ggplot2::scale_color_gradient(low = low, high = high, name = color_label)
  } else {
    p <- ggplot2::ggplot(melted, ggplot2::aes(x = Gene, y = Term, size = abs_fc, color = fc_value)) +
      ggplot2::geom_point(alpha = 0.85) +
      ggplot2::scale_size_continuous(range = c(2, 7), name = "|Fold Change|") +
      ggplot2::scale_color_gradient2(low = low, mid = mid, high = high, name = "Fold Change")
  }
  p <- p +
    ggplot2::labs(title = "Gene-Pathway Association Dotplot", x = "Genes", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                   axis.text.y = ggplot2::element_text(size = 9),
                   plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold"))

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
#  9. gggeneheat — Gene–term heatmap
# ============================================================

#' Gene–term heatmap
#'
#' Heatmap of fold-changes for genes in enriched terms.
#'
#' @param resultFis data.frame of enrichment results (columns: Term, GeneID)
#' @param fc named numeric vector of fold-changes (names = gene IDs)
#' @param top number of top terms (default 10)
#' @param max.genes max genes per term to display (default 30)
#' @param pvalue p-value cutoff
#' @param padj adjusted-p cutoff
#' @param low colour for negative FC
#' @param mid colour for zero FC
#' @param high colour for positive FC
#' @param label.cutoff absolute FC threshold above which gene labels are shown
#' @param short shorten long term names
#' @param sep separator for gene IDs
#' @param filename if non-NULL, save figure
#' @param width figure width
#' @param height figure height
#' @return a ggplot2 object
#' @keywords internal
gggeneheat_internal <- function(resultFis,
                                fc           = NULL,
                                top          = 10,
                                max.genes    = 30,
                                pvalue       = 0.05,
                                padj         = NULL,
                                usePadj      = TRUE,
                                low          = "#2166ac",
                                mid          = "white",
                                high         = "#b2182b",
                                na.fill      = "grey90",
                                border.color = "grey40",
                                label.cutoff = 1.5,
                                short        = FALSE,
                                sep          = ",",
                                filename     = NULL,
                                width        = 12,
                                height       = 8) {

  resultFis <- .filter_and_top(resultFis, top = top, pvalue = pvalue, padj = padj)
  if (nrow(resultFis) == 0) { message("No significant terms to plot."); return(invisible(NULL)) }
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  melted <- .melt_gene_term(resultFis, fc = fc, sep = sep, max.genes = max.genes)
  if (nrow(melted) == 0) { message("No gene-term pairs."); return(invisible(NULL)) }

  # Build full gene x term grid so missing combinations appear as NA (grey)
  all_genes <- unique(melted$Gene)
  all_terms <- unique(melted$Term)
  full_grid <- expand.grid(Gene = all_genes, Term = all_terms, stringsAsFactors = FALSE)
  melted$present <- TRUE
  melted <- merge(full_grid, melted, by = c("Gene", "Term"), all.x = TRUE)

  if (is.null(fc)) {
    neg_log <- if (isTRUE(usePadj) && "Padj" %in% names(resultFis)) -log10(resultFis$Padj) else -log10(resultFis$Pvalue)
    pval_map <- stats::setNames(neg_log, resultFis$Term)
    melted$neg_log_p <- ifelse(is.na(melted$present), NA_real_, pval_map[as.character(melted$Term)])
    fill_label <- if (isTRUE(usePadj)) "-log10(Padj)" else "-log10(Pvalue)"
    p <- ggplot2::ggplot(melted, ggplot2::aes(x = Gene, y = Term, fill = neg_log_p)) +
      ggplot2::geom_tile(color = border.color, linewidth = 0.5) +
      ggplot2::scale_fill_gradient(low = low, high = high, name = fill_label,
                                   na.value = na.fill)
  } else {
    melted$fc_value[is.na(melted$present)] <- NA_real_
    melted$label <- ifelse(!is.na(melted$fc_value) & abs(melted$fc_value) > label.cutoff,
                           as.character(melted$Gene), "")
    p <- ggplot2::ggplot(melted, ggplot2::aes(x = Gene, y = Term, fill = fc_value)) +
      ggplot2::geom_tile(color = border.color, linewidth = 0.5) +
      ggplot2::geom_text(data = melted[!is.na(melted$label) & melted$label != "", ],
                         ggplot2::aes(label = label), size = 2.5, color = "black") +
      ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high,
                                    name = "Fold Change", na.value = na.fill)
  }
  p <- p +
    ggplot2::labs(title = "Gene-Term Heatmap", x = "Genes", y = NULL) +
    ggplot2::theme(panel.background = ggplot2::element_rect(fill = "white", color = NA),
                   panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
                   axis.text.y = ggplot2::element_text(size = 9),
                   axis.ticks = ggplot2::element_line(color = "grey50"),
                   plot.title  = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}

# ============================================================
# 10. gggenebar — Barplot with gene names embedded inside bars
# ============================================================

#' Barplot with gene names embedded inside bars
#'
#' @param resultFis data.frame of enrichment results
#' @param top number of top terms (default 10)
#' @param pvalue p-value cutoff
#' @param padj adjusted-p cutoff
#' @param max.genes max gene names shown inside each bar (default 5)
#' @param bar.color bar colour (default \code{"#b2182b"})
#' @param bar.alpha bar transparency (default 0.85)
#' @param text.size size of gene labels inside bars (default 3)
#' @param short shorten long term names
#' @param sep separator for gene IDs
#' @param filename if non-NULL, save figure
#' @param width figure width
#' @param height figure height
#' @return a ggplot2 object
#' @keywords internal
gggenebar_internal <- function(resultFis,
                               top       = 10,
                               pvalue    = 0.05,
                               padj      = NULL,
                               usePadj   = TRUE,
                               low       = "lightpink",
                               high      = "red",
                               max.genes = 5,
                               bar.color = "#b2182b",
                               bar.alpha = 0.85,
                               text.size = 3,
                               short     = FALSE,
                               sep       = ",",
                               filename  = NULL,
                               width     = 10,
                               height    = 8) {

  resultFis <- .filter_and_top(resultFis, top = top, pvalue = pvalue, padj = padj)
  if (nrow(resultFis) == 0) { message("No significant terms to plot."); return(invisible(NULL)) }
  if (isTRUE(short)) resultFis$Term <- vapply(resultFis$Term, function(x) .paste.char(x, n = 6), character(1))

  # build gene label per term (first max.genes genes)
  resultFis$gene_label <- vapply(resultFis$GeneID, function(gid) {
    genes <- unlist(strsplit(as.character(gid), split = sep))
    paste(utils::head(genes, max.genes), collapse = ", ")
  }, character(1))

  neg_log <- if (isTRUE(usePadj) && "Padj" %in% names(resultFis)) -log10(resultFis$Padj) else -log10(resultFis$Pvalue)
  resultFis$neg_log_p <- neg_log
  if (!"RichFactor" %in% names(resultFis) && "Significant" %in% names(resultFis) && "Annotated" %in% names(resultFis))
    resultFis$RichFactor <- resultFis$Significant / resultFis$Annotated
  color_label <- if (isTRUE(usePadj)) "-log10(Padj)" else "-log10(Pvalue)"

  p <- ggplot2::ggplot(resultFis, ggplot2::aes(x = stats::reorder(Term, RichFactor), y = RichFactor)) +
    ggplot2::geom_bar(stat = "identity", ggplot2::aes(fill = neg_log_p), alpha = bar.alpha) +
    ggplot2::scale_fill_gradient(low = low, high = high, name = color_label) +
    ggplot2::geom_text(ggplot2::aes(label = gene_label, y = 0),
                       hjust = 0, nudge_y = 0.005, size = text.size, color = "black", fontface = "italic") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Enriched Terms with Key Genes", x = NULL, y = "RichFactor") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  if (!is.null(filename)) ggplot2::ggsave(filename, p, width = width, height = height)
  p
}


# ============================================================
#  Internal helpers
# ============================================================

#' Filter and subset enrichment results
#' @keywords internal
#' @noRd
.filter_and_top <- function(df, top, pvalue, padj) {
  if (!is.null(padj))  df <- df[df$Padj < padj, , drop = FALSE]
  else                 df <- df[df$Pvalue < pvalue, , drop = FALSE]
  df <- df[!is.na(df$Pvalue), , drop = FALSE]
  if (nrow(df) > top) df <- df[seq_len(top), , drop = FALSE]
  df
}

#' Standardise GSEA column names (pval->Pvalue, padj->Padj, pathway->Term)
#' @keywords internal
#' @noRd
.standardise_gsea_cols <- function(df) {
  if ("pval"    %in% names(df) && !"Pvalue" %in% names(df)) names(df)[names(df) == "pval"]    <- "Pvalue"
  if ("padj"    %in% names(df) && !"Padj"   %in% names(df)) names(df)[names(df) == "padj"]    <- "Padj"
  if ("pathway" %in% names(df) && !"Term"   %in% names(df)) names(df)[names(df) == "pathway"] <- "Term"
  if ("size"    %in% names(df) && !"Significant" %in% names(df)) names(df)[names(df) == "size"] <- "Significant"
  df
}

#' Melt gene-term pairs for dotplot / heatmap
#' @keywords internal
#' @noRd
.melt_gene_term <- function(resultFis, fc = NULL, sep = ",", max.genes = 30) {
  rows <- lapply(seq_len(nrow(resultFis)), function(i) {
    genes <- unlist(strsplit(as.character(resultFis$GeneID[i]), split = sep))
    genes <- utils::head(genes, max.genes)
    data.frame(Term = resultFis$Term[i], Gene = genes, stringsAsFactors = FALSE)
  })
  melted <- do.call(rbind, rows)
  if (is.null(fc)) {
    melted$fc_value <- 0
  } else {
    melted$fc_value <- fc[melted$Gene]
    melted$fc_value[is.na(melted$fc_value)] <- 0
  }
  melted$abs_fc <- abs(melted$fc_value)
  melted
}
