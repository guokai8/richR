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
  line_numbers <- which(nzchar(trimws(lines)))
  lines <- lines[line_numbers]
  if (length(lines) == 0) {
    stop("GMT file is empty.", call. = FALSE)
  }
  rows <- list()
  parsed_terms <- character()
  for (i in seq_along(lines)) {
    fields <- strsplit(paste0(lines[i], "\t"), "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) next
    term_name <- fields[1]
    if (is.na(term_name) || !nzchar(trimws(term_name))) {
      stop("GMT term identifier is blank on line ", line_numbers[i], ".", call. = FALSE)
    }
    parsed_terms <- c(parsed_terms, term_name)
    # fields[2] is description, skip
    genes <- fields[3:length(fields)]
    rows[[length(rows) + 1]] <- data.frame(
      GeneID = genes,
      Term = term_name,
      Annot = term_name,
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) {
    stop("No valid gene sets found in GMT file.", call. = FALSE)
  }
  duplicated_terms <- duplicated(parsed_terms) | duplicated(parsed_terms, fromLast = TRUE)
  if (any(duplicated_terms)) {
    warning("readGMT: duplicate term identifier(s) merged: ",
            paste(unique(parsed_terms[duplicated_terms]), collapse = ", "),
            call. = FALSE)
  }
  annot_df <- do.call(rbind, rows)
  normalization <- .normalize_annotation(annot_df, func_name = "readGMT")
  annot_df <- normalization$data
  n_sets <- length(unique(annot_df$Term))
  n_genes <- length(unique(annot_df$GeneID))
  message("readGMT: loaded ", n_sets, " gene sets with ", n_genes, " unique genes.")
  provenance <- .annotation_provenance(
    constructor = "readGMT", species = species, keytype = keytype,
    anntype = "custom", source_type = "custom",
    source_name = basename(file), source_version = NA_character_,
    constructor_normalization = normalization$counts
  )
  .attach_annotation_metadata(new("Annot",
      species = species,
      anntype = "custom",
      keytype = keytype,
      annot = annot_df), normalization, provenance,
      options = list(file = basename(file)))
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
  .validate_named_list(gene_sets, "gene_sets")
  rows <- lapply(seq_along(gene_sets), function(i) {
    term <- names(gene_sets)[i]
    genes <- as.character(gene_sets[[i]])
    if (length(genes) == 0L) return(NULL)
    data.frame(GeneID = genes, Term = term, Annot = term, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    stop("No valid gene sets provided.", call. = FALSE)
  }
  normalization <- .normalize_annotation(do.call(rbind, rows),
                                          func_name = "buildAnnotFromList")
  annot_df <- normalization$data
  if (is.null(annot_df) || nrow(annot_df) == 0) {
    stop("No valid gene sets provided.", call. = FALSE)
  }
  n_sets <- length(unique(annot_df$Term))
  n_genes <- length(unique(annot_df$GeneID))
  message("buildAnnotFromList: created annotation with ", n_sets,
          " gene sets and ", n_genes, " unique genes.")
  provenance <- .annotation_provenance(
    constructor = "buildAnnotFromList", species = species, keytype = keytype,
    anntype = "custom", source_type = "custom", source_name = NA_character_,
    source_version = NA_character_,
    constructor_normalization = normalization$counts
  )
  .attach_annotation_metadata(new("Annot",
      species = species,
      anntype = "custom",
      keytype = keytype,
      annot = annot_df), normalization, provenance)
}
