#' GO annotation data for human
#'
#' A data.frame containing GO annotation mappings for human genes.
#'
#' @format A data.frame with columns:
#' \describe{
#'   \item{GeneID}{Gene identifier}
#'   \item{GO}{GO term identifier}
#'   \item{Ontology}{GO ontology (BP, CC, MF)}
#' }
#' @source Generated from Bioconductor org.Hs.eg.db and GO.db
"godata"

#' KEGG pathway database
#'
#' A data.frame containing KEGG pathway and species mapping information.
#'
#' @format A data.frame with KEGG pathway identifiers and descriptions
#' @source KEGG database via KEGGREST
"kegg.db"

#' KEGG module data
#'
#' A data.frame containing KEGG module information.
#'
#' @format A data.frame with KEGG module identifiers and descriptions
#' @source KEGG database via KEGGREST
"module"

#' KEGG pathway data
#'
#' A data.frame containing KEGG pathway information.
#'
#' @format A data.frame with KEGG pathway identifiers and descriptions
#' @source KEGG database via KEGGREST
"path"

#' KEGG pathway hierarchy
#'
#' A data.frame containing KEGG pathway hierarchy with Level1 and Level2 categories.
#'
#' @format A data.frame with columns for pathway identifiers and hierarchy levels
#' @source KEGG database
"pathway"
