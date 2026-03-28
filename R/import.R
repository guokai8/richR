#' Read gene sets from a GMT file
#'
#' Reads a Gene Matrix Transposed (GMT) file and returns an Annot object
#' for use with \code{enrich()} or \code{richGSEA()}.
#'
#' @param file path to a GMT file
#' @param species species name (default: "custom")
#' @param keytype gene ID type in the file (default: "SYMBOL")
#' @return An \code{Annot} object
#'
#' @details
#' GMT format: each line is tab-separated with:
#' \enumerate{
#'   \item Gene set name
#'   \item Description (or URL)
#'   \item Gene 1, Gene 2, ... (tab-separated)
#' }
#'
#' @examples
#' \dontrun{
#' # Read MSigDB hallmark gene sets
#' annot <- readGMT("h.all.v2023.2.Hs.symbols.gmt")
#' res <- enrich(my_genes, annot)
#'
#' # Read with species info
#' annot <- readGMT("my_pathways.gmt", species = "human", keytype = "SYMBOL")
#' }
#' @export
#' @author Junguk Hur
readGMT <- function(file, species = "custom", keytype = "SYMBOL") {
  if (!file.exists(file)) {
    stop("GMT file not found: ", file, call. = FALSE)
  }

  lines <- readLines(file, warn = FALSE)
  lines <- lines[nchar(trimws(lines)) > 0]

  if (length(lines) == 0) {
    stop("GMT file is empty.", call. = FALSE)
  }

  rows <- list()
  for (line in lines) {
    fields <- strsplit(line, "\t")[[1]]
    if (length(fields) < 3) next
    term_name <- fields[1]
    # fields[2] is description, skip
    genes <- fields[3:length(fields)]
    genes <- genes[nchar(trimws(genes)) > 0]
    if (length(genes) > 0) {
      rows[[length(rows) + 1]] <- data.frame(
        GeneID = genes,
        Term = term_name,
        Annot = term_name,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(rows) == 0) {
    stop("No valid gene sets found in GMT file.", call. = FALSE)
  }

  annot_df <- do.call(rbind, rows)

  n_sets <- length(unique(annot_df$Term))
  n_genes <- length(unique(annot_df$GeneID))
  message("readGMT: loaded ", n_sets, " gene sets with ", n_genes, " unique genes.")

  new("Annot",
      species = species,
      anntype = "custom",
      keytype = keytype,
      annot = annot_df)
}


#' Build annotation from a named list of gene sets
#'
#' Converts a named list of character vectors (gene sets) into an Annot object
#' for use with \code{enrich()} or \code{richGSEA()}.
#'
#' @param gene_sets named list of character vectors, where names are gene set
#'   names and values are gene ID vectors
#' @param species species name (default: "custom")
#' @param keytype gene ID type (default: "SYMBOL")
#' @return An \code{Annot} object
#'
#' @examples
#' \dontrun{
#' my_sets <- list(
#'   "Apoptosis" = c("TP53", "BAX", "BCL2", "CASP3"),
#'   "Cell Cycle" = c("CDK1", "CDK2", "CCND1", "RB1"),
#'   "DNA Repair" = c("BRCA1", "BRCA2", "ATM", "ATR")
#' )
#' annot <- buildAnnotFromList(my_sets)
#' res <- enrich(my_genes, annot)
#' }
#' @export
#' @author Junguk Hur
buildAnnotFromList <- function(gene_sets, species = "custom", keytype = "SYMBOL") {
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list of character vectors.", call. = FALSE)
  }

  rows <- lapply(names(gene_sets), function(term) {
    genes <- unique(as.character(gene_sets[[term]]))
    genes <- genes[!is.na(genes) & nchar(genes) > 0]
    if (length(genes) == 0) return(NULL)
    data.frame(GeneID = genes, Term = term, Annot = term, stringsAsFactors = FALSE)
  })

  annot_df <- do.call(rbind, rows[!sapply(rows, is.null)])

  if (is.null(annot_df) || nrow(annot_df) == 0) {
    stop("No valid gene sets provided.", call. = FALSE)
  }

  n_sets <- length(unique(annot_df$Term))
  n_genes <- length(unique(annot_df$GeneID))
  message("buildAnnotFromList: created annotation with ", n_sets,
          " gene sets and ", n_genes, " unique genes.")

  new("Annot",
      species = species,
      anntype = "custom",
      keytype = keytype,
      annot = annot_df)
}
