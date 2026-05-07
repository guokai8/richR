#' UpSet plot for visualizing set intersections with color support
#'
#' Creates an UpSet-style plot from named gene lists, built entirely with ggplot2.
#' Supports per-set coloring of bars, matrix dots, and connecting lines, without
#' requiring or modifying the UpSetR package.
#'
#' @param x A named list of character vectors (gene sets), or a list of richResult objects
#' @param sets Character vector of set names to include (default: all sets)
#' @param order.by How to order intersections: "freq" (default) or "degree"
#' @param nsets Maximum number of sets to display (default: all)
#' @param nintersects Maximum number of intersections to display (default: 40)
#' @param mycol Character vector of colors for each set. If NULL, uses a default palette.
#' @param show.numbers Logical, show intersection size labels on bars (default: TRUE)
#' @param number.size Font size for intersection count labels (default: 3)
#' @param point.size Size of matrix dots (default: 3)
#' @param line.size Line width for connecting segments in the matrix (default: 0.8)
#' @param line.color Color for connecting lines between dots (default: "gray23")
#' @param color.by.set Logical, color matrix dots by set (default: TRUE). If FALSE,
#'   active dots use a single \code{matrix.color}.
#' @param matrix.color Single color for active dots when \code{color.by.set = FALSE}
#'   (default: "black")
#' @param inactive.color Color for inactive dots in matrix (default: "grey80")
#' @param main.bar.color Color for intersection bar chart. If NULL (default),
#'   single-set bars use the set color and multi-set intersection bars use an
#'   RGB-blend of the participating set colors.
#' @param text.scale Numeric vector of length 6 for scaling text elements:
#'   c(intersection size title, intersection size tick labels,
#'   set size title, set size tick labels, set names, numbers above bars).
#'   Or a single number to scale all uniformly.
#' @param set.size.show Logical, show set size bar chart (default: TRUE)
#' @param sep Character used to separate genes in GeneID columns (default: ",")
#' @param filename Optional output filename for saving the plot
#' @param width Plot width in inches (default: 10)
#' @param height Plot height in inches (default: 6)
#'
#' @return A combined ggplot object (via \code{cowplot::plot_grid})
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_point geom_segment geom_text
#'   coord_flip element_blank element_text labs margin scale_color_identity
#'   scale_fill_identity scale_x_continuous scale_x_discrete scale_y_continuous
#'   theme theme_minimal theme_void ggsave
#' @importFrom cowplot plot_grid align_plots ggdraw draw_plot
#' @importFrom grDevices col2rgb rgb
#'
#' @examples
#' \dontrun{
#' # From gene lists
#' sets <- list(
#'   "Treatment A" = c("TP53", "BRCA1", "EGFR", "MYC", "PTEN"),
#'   "Treatment B" = c("TP53", "EGFR", "KRAS", "AKT1", "RB1"),
#'   "Control"     = c("TP53", "MYC", "KRAS", "CDK2", "CCND1")
#' )
#' richUpset(sets)
#'
#' # Custom colors
#' richUpset(sets, mycol = c("dodgerblue", "goldenrod1", "seagreen3"))
#'
#' # From richResult objects
#' res1 <- richKEGG(gene1, kodata = hsako)
#' res2 <- richKEGG(gene2, kodata = hsako)
#' richUpset(list("Sample1" = res1, "Sample2" = res2))
#' }
#' @export
#' @author Junguk Hur
richUpset <- function(x,
                    sets = NULL,
                    order.by = c("freq", "degree"),
                    nsets = NULL,
                    nintersects = 40,
                    mycol = NULL,
                    show.numbers = TRUE,
                    number.size = 3,
                    point.size = 3,
                    line.size = 0.8,
                    line.color = "gray23",
                    color.by.set = TRUE,
                    matrix.color = "black",
                    inactive.color = "grey80",
                    main.bar.color = NULL,
                    text.scale = 1,
                    set.size.show = TRUE,
                    sep = ",",
                    filename = NULL,
                    width = 10,
                    height = 6) {
  order.by <- match.arg(order.by)

  # --- Extract gene lists ---
  gene_lists <- .extract_gene_lists(x, sep = sep)
  if (is.null(names(gene_lists))) {
    names(gene_lists) <- paste0("Set", seq_along(gene_lists))
  }

  # Subset and limit sets
  if (!is.null(sets)) {
    gene_lists <- gene_lists[sets]
  }
  if (!is.null(nsets) && length(gene_lists) > nsets) {
    set_sizes <- sapply(gene_lists, length)
    gene_lists <- gene_lists[order(set_sizes, decreasing = TRUE)[1:nsets]]
  }
  n_sets <- length(gene_lists)
  set_names <- names(gene_lists)

  # --- Set colors ---
  default_pal <- c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3",
                    "orchid3", "#E31A1C", "#1F78B4", "#33A02C", "#FF7F00",
                    "#6A3D9A", "#B15928", "#A6CEE3")
  if (is.null(mycol)) {
    mycol <- rep(default_pal, length.out = n_sets)
  } else {
    mycol <- rep(mycol, length.out = n_sets)
  }
  names(mycol) <- set_names

  # --- Compute all intersections ---
  all_genes <- unique(unlist(gene_lists))
  membership <- sapply(gene_lists, function(s) all_genes %in% s)
  rownames(membership) <- all_genes

  # Generate all unique intersection patterns
  patterns <- unique(membership)
  pattern_ids <- apply(membership, 1, paste, collapse = "")
  pattern_key <- apply(patterns, 1, paste, collapse = "")

  intersect_data <- data.frame(
    pattern = pattern_key,
    count = as.integer(table(pattern_ids)[pattern_key]),
    degree = rowSums(patterns),
    stringsAsFactors = FALSE
  )

  # Add set membership columns
  for (i in seq_len(n_sets)) {
    intersect_data[[set_names[i]]] <- patterns[, i]
  }

  # Remove empty intersection (no sets)
  intersect_data <- intersect_data[intersect_data$degree > 0, ]

  # Order intersections
  if (order.by == "freq") {
    intersect_data <- intersect_data[order(intersect_data$count, decreasing = TRUE), ]
  } else {
    intersect_data <- intersect_data[order(intersect_data$degree, intersect_data$count,
                                            decreasing = c(TRUE, TRUE)), ]
  }

  # Limit number of intersections
  if (nrow(intersect_data) > nintersects) {
    intersect_data <- intersect_data[1:nintersects, ]
  }
  n_intersects <- nrow(intersect_data)
  intersect_data$x <- seq_len(n_intersects)

  # --- Text scaling ---
  if (length(text.scale) == 1) {
    text.scale <- rep(text.scale, 6)
  }
  ts <- text.scale * c(10, 8, 10, 8, 9, 3)

  # --- Determine bar colors (single-set: set color; multi-set: blended) ---
  bar_colors <- sapply(seq_len(n_intersects), function(i) {
    if (!is.null(main.bar.color)) return(main.bar.color)
    active_sets <- set_names[as.logical(intersect_data[i, set_names])]
    if (length(active_sets) == 1) return(mycol[active_sets])
    .blend_colors(mycol[active_sets])
  })

  # --- 1. Intersection size bar chart (top) ---
  bar_df <- data.frame(
    x = intersect_data$x,
    count = intersect_data$count,
    fill_color = bar_colors,
    stringsAsFactors = FALSE
  )

  p_bars <- ggplot(bar_df, aes(x = x, y = count)) +
    geom_bar(stat = "identity", aes(fill = fill_color), width = 0.7) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(0.5, n_intersects + 0.5), expand = c(0, 0)) +
    labs(y = "Intersection Size") +
    theme_minimal() +
    theme(
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(size = ts[1]),
      axis.text.y = element_text(size = ts[2]),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(5, 5, 0, 5)
    )

  if (show.numbers) {
    p_bars <- p_bars + geom_text(
      aes(label = count, color = fill_color),
      vjust = -0.3, size = ts[6],
      show.legend = FALSE
    ) +
      scale_color_identity()
  }

  # --- 2. Matrix dot plot (middle) ---
  matrix_rows <- list()
  for (i in seq_len(n_intersects)) {
    for (j in seq_len(n_sets)) {
      active <- as.logical(intersect_data[i, set_names[j]])
      if (color.by.set) {
        dot_color <- if (active) mycol[set_names[j]] else inactive.color
      } else {
        dot_color <- if (active) matrix.color else inactive.color
      }
      matrix_rows[[length(matrix_rows) + 1]] <- data.frame(
        x = intersect_data$x[i],
        y = j,
        active = active,
        color = dot_color,
        set = set_names[j],
        stringsAsFactors = FALSE
      )
    }
  }
  matrix_df <- do.call(rbind, matrix_rows)

  # Connecting lines for active dots in each intersection
  line_rows <- list()
  for (i in seq_len(n_intersects)) {
    active_y <- which(as.logical(intersect_data[i, set_names]))
    if (length(active_y) >= 2) {
      line_rows[[length(line_rows) + 1]] <- data.frame(
        x = intersect_data$x[i],
        ymin = min(active_y),
        ymax = max(active_y),
        stringsAsFactors = FALSE
      )
    }
  }

  p_matrix <- ggplot(matrix_df, aes(x = x, y = y)) +
    geom_point(aes(color = color), size = point.size) +
    scale_color_identity() +
    scale_x_continuous(limits = c(0.5, n_intersects + 0.5), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = seq_len(n_sets),
      labels = set_names,
      limits = c(0.5, n_sets + 0.5),
      expand = c(0, 0)
    ) +
    theme_minimal() +
    theme(
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = ts[5], color = mycol[set_names]),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin = margin(0, 5, 5, 5)
    )

  if (length(line_rows) > 0) {
    line_df <- do.call(rbind, line_rows)
    p_matrix <- p_matrix +
      geom_segment(
        data = line_df,
        aes(x = x, xend = x, y = ymin, yend = ymax),
        color = line.color,
        linewidth = line.size,
        inherit.aes = FALSE
      )
    # Redraw active dots on top of lines
    active_df <- matrix_df[matrix_df$active, ]
    p_matrix <- p_matrix +
      geom_point(data = active_df, aes(x = x, y = y, color = color),
                 size = point.size, inherit.aes = FALSE)
  }

  # --- 3. Set size bar chart (left) ---
  if (set.size.show) {
    set_size_df <- data.frame(
      set = factor(set_names, levels = set_names),
      y = seq_len(n_sets),
      size = sapply(gene_lists, length),
      fill_color = mycol[set_names],
      stringsAsFactors = FALSE
    )

    p_setsize <- ggplot(set_size_df, aes(x = y, y = size)) +
      geom_bar(stat = "identity", aes(fill = fill_color), width = 0.7) +
      scale_fill_identity() +
      scale_x_continuous(
        breaks = seq_len(n_sets),
        labels = NULL,
        limits = c(0.5, n_sets + 0.5),
        expand = c(0, 0)
      ) +
      coord_flip() +
      labs(y = "Set Size") +
      theme_minimal() +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = ts[3]),
        axis.text.x = element_text(size = ts[4]),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(0, 0, 5, 5)
      )

    # Align the bar chart and matrix on their shared intersection x-axis so
    # their panels have identical widths, then place them with ggdraw for
    # pixel-accurate positioning regardless of y-axis label width.
    aligned <- align_plots(p_bars, p_matrix, align = "v", axis = "lr")
    set_w  <- 0.15
    main_w <- 1 - set_w
    top_h  <- 0.6
    bot_h  <- 1 - top_h

    final_plot <- ggdraw() +
      draw_plot(aligned[[1]], x = set_w, y = bot_h, width = main_w, height = top_h) +
      draw_plot(p_setsize,    x = 0,     y = 0,     width = set_w,  height = bot_h) +
      draw_plot(aligned[[2]], x = set_w, y = 0,     width = main_w, height = bot_h)
  } else {
    final_plot <- plot_grid(p_bars, p_matrix, ncol = 1, rel_heights = c(0.6, 0.4),
                            align = "v", axis = "lr")
  }

  if (!is.null(filename)) {
    ggsave(filename, final_plot, width = width, height = height)
    message("Plot saved to: ", filename)
  }

  return(final_plot)
}

#' Blend multiple colors by averaging RGB channels
#'
#' Used by \code{richUpset} to color multi-set intersection bars. Given a vector
#' of color specifications, returns a single hex color whose R/G/B channels are
#' the mean of the inputs. Length-1 input is returned unchanged.
#'
#' @param colors Character vector of color specifications (any form accepted by
#'   \code{grDevices::col2rgb}).
#' @return A single hex color string.
#' @keywords internal
.blend_colors <- function(colors) {
  if (length(colors) == 1) return(colors)
  rgbs <- col2rgb(colors)
  blended <- round(rowMeans(rgbs))
  rgb(blended[1], blended[2], blended[3], maxColorValue = 255)
}

#' Extract gene lists from various input types
#' @param x Named list of character vectors, richResult objects, or GSEAResult objects
#' @param sep Separator for gene IDs
#' @return Named list of character vectors
#' @keywords internal
.extract_gene_lists <- function(x, sep = ",") {
  if (!is.list(x)) {
    stop("x must be a named list of gene sets or enrichment results")
  }
  lapply(x, function(item) {
    if (is.character(item)) {
      return(unique(item))
    }
    if (inherits(item, "richResult")) {
      return(unique(item@gene))
    }
    if (inherits(item, "GSEAResult")) {
      return(unique(item@gene))
    }
    if (is.data.frame(item)) {
      if ("GeneID" %in% colnames(item)) {
        return(unique(unlist(strsplit(as.character(item$GeneID), sep))))
      }
      if ("leadingEdge" %in% colnames(item)) {
        return(unique(unlist(strsplit(as.character(item$leadingEdge), sep))))
      }
    }
    stop("Unsupported input type. Provide character vectors, richResult, GSEAResult, or data.frames with GeneID column.")
  })
}
