#' Plot KEGG Cluster Visualization
#'
#' This function merges two cluster-plot approaches:
#' \itemize{
#'   \item \strong{Enrichment mode}: requires columns \code{Term}, \code{Padj}, \code{Significant}, \code{Annotated}, etc.
#'   \item \strong{GSEA mode}: requires columns \code{pathway}, \code{padj}, \code{NES}, etc.
#' }
#'
#' If \code{method = "auto"}, the function guesses which mode to use based on the columns found.
#'
#' @param data A data frame with the appropriate columns for either mode.
#' @param method One of \code{"enrich"} or \code{"gsea"}. If \code{"auto"}, guessed from columns.
#' @param data A data frame containing KEGG cluster data. Must include columns: `Term`, `Level2`, `group`, `Padj`, `Significant`, and `Annotated`.
#' @param color_low Color for the lowest value in the gradient (default: "pink").
#' @param color_high Color for the highest value in the gradient (default: "red").
#' @param color_mid Color for the middle value in the gradient (default: "white").
#' @param size_range A numeric vector of length 2 to control the size range of points (default: c(0.8, 4)).
#' @param curve_color Color for the connecting curves between Level2 and Terms (default: "grey70").
#' @param curve_size Line width for the connecting curves (default: 0.5).
#' @param vertical_line_color Color for the vertical lines between Level2 labels and Terms (default: "darkcyan").
#' @param vertical_line_size Width of vertical lines between Level2 and Terms (default: 1.5).
#' @param dot_line_color Color for dotted lines connecting clusters to pathways (default: "grey70").
#' @param dot_line_size Line width for the dotted lines (default: 0.3).
#' @param dot_line_type Line type for the dotted lines (default: "dotted").
#' @param vline_color Color for vertical lines separating clusters (default: "grey80").
#' @param vline_type Line type for vertical lines separating clusters (default: "dashed").
#' @param label_font_size Font size for the Level2 labels (default: 3).
#' @param label_font_face Font face for the Level2 labels (default: "bold").
#' @param pathway_font_size Font size for pathway labels (default: 2.5).
#' @param pathway_font_face Font face for pathway labels (default: "italic").
#' @param legend_position Position of the legend (default: "right").
#' @param x_pathway_offset Offset to control the x-axis positioning of pathways (default: 0.3).
#' @param label_x_size Font size for the x labels (default: 6).
#' @param plot_margins Numeric vector of length 4 to control plot margins (default: c(5, 200, 10, 250)).
#'
#' @return A ggplot2 object (or a \code{cowplot} combined object).
#'
#' @export
ggcluster <- function(data, method = c("auto", "enrich", "gsea"),color_low = "pink",
    color_high = "red",
    color_mid = "white",
    size_range = c(0.8, 4),
    curve_color = "grey70",
    curve_size = 0.5,
    vertical_line_color = "darkcyan",
    vertical_line_size = 1.5,
    dot_line_color = "grey70",
    dot_line_size = 0.3,
    dot_line_type = "dotted",
    vline_color = "grey80",
    vline_type = "dashed",
    label_font_size = 3,
    label_font_face = "bold",
    pathway_font_size = 2.5,
    pathway_font_face = "italic",
    legend_position = "right",
    x_pathway_offset = 0.3,
    label_x_size = 6,
    plot_margins = c(5, 200, 10, 250),...) {
  method <- match.arg(method)
  if (method == "auto") {
    if ("Term" %in% colnames(data)) {
      method <- "enrich"
    } else if ("pathway" %in% colnames(data)) {
      method <- "gsea"
    } else {
      stop("Cannot autodetect method, please specify method='enrich' or 'gsea'.")
    }
  }

  if (method == "enrich") {
    .ggcluster_enrich(data, color_low = color_low,
    color_high = color_high,
    size_range = size_range,
    curve_color = curve_color,
    curve_size = curve_size,
    vertical_line_color = vertical_line_color,
    vertical_line_size = vertical_line_size,
    dot_line_color = dot_line_color,
    dot_line_size = dot_line_size,
    dot_line_type = dot_line_type,
    vline_color =vline_color,
    vline_type = vline_type,
    label_font_size = label_font_size,
    label_font_face = label_font_face,
    pathway_font_size = pathway_font_size,
    pathway_font_face = pathway_font_face,
    legend_position = legend_position,
    x_pathway_offset = x_pathway_offset,
    label_x_size = label_x_size,
    plot_margins = plot_margins, ...)
  } else {
    .ggcluster_gsea(data,color_low = color_low,
    color_high = color_high,
    color_mid = color_mid,
    size_range = size_range,
    curve_color = curve_color,
    curve_size = curve_size,
    vertical_line_color = vertical_line_color,
    vertical_line_size = vertical_line_size,
    dot_line_color = dot_line_color,
    dot_line_size = dot_line_size,
    dot_line_type = dot_line_type,
    vline_color =vline_color,
    vline_type = vline_type,
    label_font_size = label_font_size,
    label_font_face = label_font_face,
    pathway_font_size = pathway_font_size,
    pathway_font_face = pathway_font_face,
    legend_position = legend_position,
    x_pathway_offset = x_pathway_offset,
    #label_x_size = label_x_size,
    plot_margins = plot_margins, ...)
  }
}



#' Plot KEGG Cluster Visualization
#'
#' This function generates a KEGG cluster plot, showing pathways and their relationships
#' across different clusters with customizable aesthetics.
#'
#' @param data A data frame containing KEGG cluster data. Must include columns: `Term`, `Level2`, `group`, `Padj`, `Significant`, and `Annotated`.
#' @param color_low Color for the lowest value in the gradient (default: "pink").
#' @param color_high Color for the highest value in the gradient (default: "red").
#' @param size_range A numeric vector of length 2 to control the size range of points (default: c(0.8, 4)).
#' @param curve_color Color for the connecting curves between Level2 and Terms (default: "grey70").
#' @param curve_size Line width for the connecting curves (default: 0.5).
#' @param vertical_line_color Color for the vertical lines between Level2 labels and Terms (default: "darkcyan").
#' @param vertical_line_size Width of vertical lines between Level2 and Terms (default: 1.5).
#' @param dot_line_color Color for dotted lines connecting clusters to pathways (default: "grey70").
#' @param dot_line_size Line width for the dotted lines (default: 0.3).
#' @param dot_line_type Line type for the dotted lines (default: "dotted").
#' @param vline_color Color for vertical lines separating clusters (default: "grey80").
#' @param vline_type Line type for vertical lines separating clusters (default: "dashed").
#' @param label_font_size Font size for the Level2 labels (default: 3).
#' @param label_font_face Font face for the Level2 labels (default: "bold").
#' @param pathway_font_size Font size for pathway labels (default: 2.5).
#' @param pathway_font_face Font face for pathway labels (default: "italic").
#' @param legend_position Position of the legend (default: "right").
#' @param x_pathway_offset Offset to control the x-axis positioning of pathways (default: 0.3).
#' @param label_x_size Font size for the x labels (default: 6).
#' @param plot_margins Numeric vector of length 4 to control plot margins (default: c(5, 200, 10, 250)).
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import cowplot
#'
#' @examples
#' data <- data.frame(
#'     Term = c("Peroxisome", "PPAR signaling pathway", "Fatty acid elongation"),
#'     Level2 = c("Transport and catabolism", "Endocrine system", "Lipid metabolism"),
#'     group = c("Cluster1", "Cluster1", "Cluster1"),
#'     Padj = c(0.015, 0.24, 0.23),
#'     Significant = c(7, 5, 3),
#'     Annotated = c(87, 89, 29)
#' )
#' .ggcluster_enrich(data)
#'
.ggcluster_enrich <- function(
    data,
    color_low = "pink",
    color_high = "red",
    size_range = c(0.8, 4),
    curve_color = "grey70",
    curve_size = 0.5,
    vertical_line_color = "darkcyan",
    vertical_line_size = 1.5,
    dot_line_color = "grey70",
    dot_line_size = 0.3,
    dot_line_type = "dotted",
    vline_color = "grey80",
    vline_type = "dashed",
    label_font_size = 3,
    label_font_face = "bold",
    pathway_font_size = 2.5,
    pathway_font_face = "italic",
    legend_position = "right",
    x_pathway_offset = 0.3,
    label_x_size = 6,
    plot_margins = c(5, 200, 10, 250)
) {
  # Check that required columns are present
  required_cols <- c("Term", "Level2", "group", "Padj", "Significant", "Annotated")
  if (!all(required_cols %in% colnames(data))) {
    stop("Input data must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Preprocess data
  data <- data %>%
    mutate(
      Padj = as.numeric(Padj),
      neg_log10_Padj = -log10(Padj)
    )

  # Define terms and positions for y-axis
  terms <- data %>%
    select(Term, Level2) %>%
    distinct() %>%
    arrange(Level2, Term) %>%
    mutate(y = row_number())

  level2_labels <- terms %>%
    group_by(Level2) %>%
    summarise(y_level2 = mean(y)) %>%
    ungroup()

  data <- data %>%
    left_join(terms %>% select(Term, y), by = "Term") %>%
    left_join(level2_labels, by = "Level2")

  groups <- sort(unique(data$group))
  n_groups <- length(groups)
  group_x_positions <- 2:(1 + n_groups)

  group_mapping <- data.frame(
    group = groups,
    x_group = group_x_positions
  )

  data <- data %>%
    left_join(group_mapping, by = "group")

  x_pathway <- 1 + n_groups + x_pathway_offset

  pathway_labels <- terms %>%
    mutate(x_pathway = x_pathway) %>%
    mutate(x = x_pathway, label = Term)

  # Prepare data for curves and segments
  lines_data <- data %>%
    select(Level2, Term, y) %>%
    distinct() %>%
    left_join(level2_labels, by = "Level2") %>%
    group_by(Level2) %>%
    arrange(y) %>%
    mutate(
      n_terms = n(),
      index = row_number(),
      curvature = case_when(
        n_terms == 1 ~ 0,
        index <= n_terms / 2 ~ +0.1,
        TRUE ~ -0.1
      )
    ) %>%
    ungroup() %>%
    mutate(
      x_start = 1,
      x_end = 1.5,
      y_start = y_level2,
      y_end = y
    )

  curve_lines_data <- lines_data %>% filter(n_terms > 1)
  segment_lines_data <- lines_data %>% filter(n_terms == 1)

  terms_lines <- terms %>%
    group_by(Level2) %>%
    summarise(
      y_min = min(y),
      y_max = max(y),
      n_terms = n()
    ) %>%
    ungroup() %>%
    mutate(
      delta = ifelse(n_terms == 1, 0.1, 0),
      y_start = y_min - delta,
      y_end = y_max + delta,
      x = 1.5
    )

  lines_horizontal <- data %>%
    select(group, x_group, y, Term) %>%
    distinct() %>%
    mutate(
      x_start = x_group,
      x_end = x_pathway,
      y_end = y
    )

  vlines_data <- data.frame(
    x = group_x_positions,
    y_start = min(data$y) - 0.5,
    y_end = max(data$y) + 0.5
  )

  p_main <- ggplot() +
    lapply(1:nrow(curve_lines_data), function(i) {
      geom_curve(
        data = curve_lines_data[i, ],
        aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
        curvature = curve_lines_data$curvature[i],
        color = curve_color,
        size = curve_size
      )
    }) +
    geom_segment(data = segment_lines_data,
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 color = curve_color, size = curve_size) +
    geom_segment(data = terms_lines,
                 aes(x = x, y = y_start, xend = x, yend = y_end),
                 color = vertical_line_color, size = vertical_line_size) +
    geom_segment(data = lines_horizontal,
                 aes(x = x_start, y = y, xend = x_end, yend = y_end),
                 color = dot_line_color, size = dot_line_size, linetype = dot_line_type) +
    geom_segment(data = vlines_data,
                 aes(x = x, y = y_start, xend = x, yend = y_end),
                 color = vline_color, linetype = vline_type) +
    geom_text(data = level2_labels,
              aes(x = 1, y = y_level2, label = Level2),
              hjust = 1, size = label_font_size, fontface = label_font_face) +
    geom_point(data = data,
               aes(x = x_group, y = y, color = neg_log10_Padj, size = Significant / Annotated)) +
    geom_text(data = pathway_labels,
              aes(x = x, y = y, label = label),
              hjust = 0, size = pathway_font_size, fontface = pathway_font_face) +
    scale_x_continuous(
      limits = c(0.5, x_pathway + 1.2),
      breaks = c(1, 1.5, group_x_positions, x_pathway),
      labels = c("", "", groups, ""),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(min(data$y) - 0.5, max(data$y) + 0.5),
      expand = c(0, 0)
    ) +
    scale_color_gradient(low = color_low, high = color_high, name = "-log10(Padj)") +
    scale_size_continuous(name = "Significant/Annotated", range = size_range) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle=90,vjust=0.5,hjust=1,size=label_x_size),
      axis.ticks = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = plot_margins[1], r = plot_margins[2], b = plot_margins[3], l = plot_margins[4], unit = "pt")
    ) +
    coord_cartesian(clip = "off")

  p_legend <- ggplot(data, aes(x = x_group, y = y, color = neg_log10_Padj, size = Significant / Annotated)) +
    geom_point() +
    scale_color_gradient(low = color_low, high = color_high, name = "-log10(Padj)") +
    scale_size_continuous(name = "RichFactor", range = size_range) +
    theme_minimal() +
    theme(
      legend.position = legend_position,
      legend.box.margin = margin(0, 0, 0, 0)
    )

  legend <- get_legend(p_legend)

  final_plot <- plot_grid(
    p_main,
    legend,
    ncol = 2,
    rel_widths = c(1, 0.2)
  )

  return(final_plot)
}
