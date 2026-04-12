## Suppress R CMD check NOTEs for non-standard evaluation variables
## used in dplyr pipelines and ggplot2 aes() calls
utils::globalVariables(c(
  # ggplot2 aes variables
  "x", "y", "xend", "yend", "label", "Group", "NES",
  "x_start", "y_start", "x_end", "y_end", "x_group", "y_level2",
  "neg_log10_Padj", "annotateText", "diff", "rich", "variable",
  "hjustvar", "vjustvar", "xpos", "ypos",
  # dplyr pipeline variables
  "Padj", "Pvalue", "Annot", "Term", "GeneID", "pathway", "Level2", "pval", "leadingEdge",
  "group", "padj", "Annotated", "Significant", "RichFactor",
  "FoldEnrichment", "n_terms", "n_pathways", "delta", "y_min", "y_max",
  # ggnetwork / ggrich
  "weight", "name", "value", "from", "to",
  # ggbar
  "neg_log_p",
  # makeAnnot / Reactome
  "gs_cat", "gs_subcat", "reactomePATHNAME2ID", "reactomePATHID2EXTID",
  # ggupset
  "fill_color", "color", "ymin", "ymax", "size", "active", "count",
  # data variables
  "kegg.db", "module", "path",
  # GO.db / AnnotationDbi
  "ONTOLOGY", "mapIds",
  # compareGSEA
  "join",
  # richDAVID (RDAVIDWebService)
  "DAVIDWebService", "getIdTypes", "addList",
  "setAnnotationCategories", "getFunctionalAnnotationChart", "getSpecieNames",
  # new ggplots
  "is_sig", "direction", "cum_rank", "Dim1", "Dim2",
  "Gene", "abs_fc", "fc_value", "gene_label", "present", "effect"
))

#' Internal richR functions and methods
#'
#' Internal functions, Rcpp helpers, and dplyr S3 method dispatchers.
#' These are not intended to be called directly by the user.
#'
#' @name richR-internal
#' @useDynLib richR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @aliases fast_factor hyper_bench_vector name_table reverseList sf
#'   arrange.Annot arrange.GSEAResult arrange.richResult
#'   as.data.frame.Annot as.data.frame.GSEAResult as.data.frame.richResult
#'   filter.Annot filter.GSEAResult filter.richResult
#'   group_by.Annot group_by.GSEAResult group_by.richResult
#'   mutate.Annot mutate.GSEAResult mutate.richResult
#'   rename.Annot rename.GSEAResult rename.richResult
#'   select.Annot select.GSEAResult select.richResult
#'   slice.Annot slice.GSEAResult slice.richResult
#'   summarise.Annot summarise.GSEAResult summarise.richResult
#' @keywords internal
NULL

.onLoad <- function(libname, pkgname) {
  invisible(NULL)
}
