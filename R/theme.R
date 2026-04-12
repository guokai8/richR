#' richR publication-ready ggplot2 theme
#'
#' A clean, consistent theme for all richR visualizations. Designed for
#' publication-quality figures with legible fonts and minimal clutter.
#'
#' @param base_size base font size (default: 12)
#' @param base_family base font family (default: "")
#' @return A ggplot2 theme object
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' richDot(res) + theme_richR()
#'
#' # Larger fonts for presentations
#' richBar(res) + theme_richR(base_size = 16)
#' }
#' @importFrom ggplot2 theme theme_minimal element_text element_line element_blank
#'   element_rect margin %+replace%
#' @export
#' @author Junguk Hur
theme_richR <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Text
      plot.title = element_text(size = base_size * 1.2, face = "bold",
                                 hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(size = base_size, face = "bold"),
      axis.text = element_text(size = base_size * 0.85, color = "gray20"),
      legend.title = element_text(size = base_size * 0.9, face = "bold"),
      legend.text = element_text(size = base_size * 0.8),
      # Grid
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      # Background
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      # Legend
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      # Margins
      plot.margin = margin(10, 10, 10, 10)
    )
}
