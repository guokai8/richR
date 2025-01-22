#' Plot KEGG Cluster Visualization
#'
#' This function generates a KEGG GESA cluster plot, showing pathways and their relationships
#' across different clusters with customizable aesthetics.
#'
#' @param data A data frame containing KEGG cluster data. Must include columns: `pathway`, `Level2`, `group`, `Padj`, `Significant`, and `Annotated`.
#' @param color_low Color for the lowest value in the gradient (default: "cyan4").
#' @param color_high Color for the highest value in the gradient (default: "white").
#' @param color_mid Color for the middle value in the gradient (default: "red").
#' @param size_range A numeric vector of length 2 to control the size range of points (default: c(0.8, 4)).
#' @param curve_color Color for the connecting curves between Level2 and pathways (default: "grey70").
#' @param curve_size Line width for the connecting curves (default: 0.5).
#' @param vertical_line_color Color for the vertical lines between Level2 labels and pathways (default: "darkcyan").
#' @param vertical_line_size Width of vertical lines between Level2 and pathways (default: 1.5).
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
#' @param plot_margins Numeric vector of length 4 to control plot margins (default: c(5, 200, 10, 250)).
#'
#' @import ggplot2
#' @import dplyr
#' @import tidyr
#' @import cowplot
#'
.ggcluster_gsea <- function(
    data,
    color_low = "cyan4",
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
    plot_margins = c(5, 200, 10, 250)
) {
  # Check that required columns are present
  required_cols <- c("pathway", "Level2", "group", "padj","NES")
  if (!all(required_cols %in% colnames(data))) {
    stop("Input data must contain columns: ", paste(required_cols, collapse = ", "))
  }

  # Preprocess data
  data <- data %>%
    mutate(
      Padj = as.numeric(padj),
      neg_log10_Padj = -log10(padj)
    )

  # Define pathways and positions for y-axis
  pathways <- data %>%
    select(pathway, Level2) %>%
    distinct() %>%
    arrange(Level2, pathway) %>%
    mutate(y = row_number())

  level2_labels <- pathways %>%
    group_by(Level2) %>%
    summarise(y_level2 = mean(y)) %>%
    ungroup()

  data <- data %>%
    left_join(pathways %>% select(pathway, y), by = "pathway") %>%
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

  pathway_labels <- pathways %>%
    mutate(x_pathway = x_pathway) %>%
    mutate(x = x_pathway, label = pathway)

  # Prepare data for curves and segments
  lines_data <- data %>%
    select(Level2, pathway, y) %>%
    distinct() %>%
    left_join(level2_labels, by = "Level2") %>%
    group_by(Level2) %>%
    arrange(y) %>%
    mutate(
      n_pathways = n(),
      index = row_number(),
      curvature = case_when(
        n_pathways == 1 ~ 0,
        index <= n_pathways / 2 ~ +0.1,
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

  curve_lines_data <- lines_data %>% filter(n_pathways > 1)
  segment_lines_data <- lines_data %>% filter(n_pathways == 1)

  pathways_lines <- pathways %>%
    group_by(Level2) %>%
    summarise(
      y_min = min(y),
      y_max = max(y),
      n_pathways = n()
    ) %>%
    ungroup() %>%
    mutate(
      delta = ifelse(n_pathways == 1, 0.1, 0),
      y_start = y_min - delta,
      y_end = y_max + delta,
      x = 1.5
    )

  lines_horizontal <- data %>%
    select(group, x_group, y, pathway) %>%
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
    geom_segment(data = pathways_lines,
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
               aes(x = x_group, y = y, color = NES, size = -log10(padj))) +
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
    scale_color_gradient2(low = color_low, high = color_high, mid=color_mid,midpoint = 0,name = "NES") +
    scale_size_continuous(name = "-log10(Padj)", range = size_range) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.text.x = element_text(angle=90,hjust=1,vjust = 0.5),
      axis.ticks = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = plot_margins[1], r = plot_margins[2], b = plot_margins[3], l = plot_margins[4], unit = "pt")
    ) +
    coord_cartesian(clip = "off")

    p_legend <- ggplot(data, aes(x = x_group, y = y, color = NES, size = neg_log10_Padj)) +
    geom_point() +
    scale_color_gradient(low = color_low, high = color_high, name = "NES") +
    scale_size_continuous(name = "-log10(Padj)", range = size_range) +
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
