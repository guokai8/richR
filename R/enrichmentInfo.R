#' Read enrichment diagnostics and provenance
#'
#' Metadata is recorded on newly created Annot, richResult, and GSEAResult
#' objects. Legacy objects return an explicit unavailable record.
#' @param x Annot, richResult, or GSEAResult object
#' @return A named list containing diagnostics and provenance
#' @export
enrichmentInfo <- function(x) {
  supported <- methods::is(x, "Annot") || methods::is(x, "richResult") ||
    methods::is(x, "GSEAResult")
  if (!supported) {
    stop("enrichmentInfo() supports Annot, richResult, or GSEAResult objects.",
         call. = FALSE)
  }
  metadata <- attr(x, "richR_metadata", exact = TRUE)
  if (is.null(metadata)) {
    return(list(schema_version = 1L, available = FALSE,
                reason = "metadata not recorded", diagnostics = NULL,
                provenance = NULL))
  }
  metadata
}

.richr_version <- function() {
  if ("richR" %in% loadedNamespaces()) {
    return(as.character(utils::packageVersion("richR")))
  }
  NA_character_
}

.package_version_or_na <- function(package) {
  if (length(package) == 1L && !is.na(package) && nzchar(package) &&
      requireNamespace(package, quietly = TRUE)) {
    as.character(utils::packageVersion(package))
  } else {
    NA_character_
  }
}

.annotation_provenance <- function(object = NULL, constructor = NULL,
                                   species = NULL, keytype = NULL,
                                   anntype = NULL, source_type = NULL,
                                   source_name = NULL, source_version = NA_character_,
                                   provider_package = NULL,
                                   provider_version = NULL,
                                   constructor_normalization = NULL) {
  if (methods::is(object, "Annot")) {
    previous <- attr(object, "richR_metadata", exact = TRUE)
    if (is.list(previous) && is.list(previous$provenance$annotation)) {
      annotation <- previous$provenance$annotation
      annotation$species <- object@species
      annotation$keytype <- object@keytype
      annotation$anntype <- object@anntype
      return(annotation)
    }
    species <- object@species
    keytype <- object@keytype
    anntype <- object@anntype
    source_type <- source_type %||% "unknown"
    source_name <- source_name %||% NA_character_
  }
  list(
    constructor = constructor,
    species = species,
    keytype = keytype,
    anntype = anntype,
    source_type = source_type %||% "custom",
    source_name = source_name %||% NA_character_,
    source_version = source_version,
    provider_package = provider_package,
    provider_version = provider_version,
    constructor_normalization = constructor_normalization
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

.annotation_metadata <- function(normalization, provenance, options = list()) {
  list(schema_version = 1L, available = TRUE, analysis = "annotation",
       annotation = normalization$counts,
       provenance = list(richR_version = .richr_version(),
                         R_version = as.character(getRversion()),
                         fgsea_version = NULL,
                         annotation = provenance, rng = NULL),
       options = options)
}

.attach_annotation_metadata <- function(object, normalization, provenance,
                                        options = list()) {
  attr(object, "richR_metadata") <- .annotation_metadata(
    normalization, provenance, options
  )
  object
}

.analysis_metadata <- function(analysis, input, annotation, background,
                               testing, options, annotation_provenance,
                               fgsea_version = NULL, rng = NULL) {
  list(schema_version = 1L, available = TRUE, analysis = analysis,
       input = input, annotation = annotation, background = background,
       testing = testing, options = options,
       provenance = list(richR_version = .richr_version(),
                         R_version = as.character(getRversion()),
                         fgsea_version = fgsea_version,
                         annotation = annotation_provenance,
                         rng = rng))
}
