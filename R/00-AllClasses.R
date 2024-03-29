##' Class "richResult"
##' This class represents the result of enrichment analysis.
##'
##'
##' @name richResult-class
##' @aliases richResult-class
##'   show,richResult-method plot,richResult-method
##'   summary,richResult-method
##'
##' @docType class
##' @slot result enrichment analysis results
##' @slot detail genes included in significant terms and original information
##' @slot pvalueCutoff cutoff pvalue
##' @slot pAdjustMethod pvalue adjust method
##' @slot padjCutoff pvalue adjust cutoff value
##' @slot genenumber number of input genes
##' @slot organism organism used
##' @slot ontology biological ontology
##' @slot gene Gene IDs
##' @slot keytype Gene ID type
##' @slot sep character string used to separate the genes when concatenating
##' @exportClass richResult
##' @author Kai Guo
##' @keywords classes
setClass("richResult",
         representation=representation(
           result         = "data.frame",
           detail         = "data.frame",
           pvalueCutoff   = "numeric",
           pAdjustMethod  = "character",
           padjCutoff   = "numeric",
           genenumber    = "numeric",
           organism       = "character",
           ontology       = "character",
           gene           = "character",
           keytype        = "character",
           sep = "character"
         )
)
##' Class "GSEAResult"
##' This class represents the result from GSEA analysis
##' @name GSEAResult-class
##' @aliases GSEAResult-class
##'   show,GSEAResult-method summary,GSEAResult-method
##'   plot,GSEAResult-method
##'
##' @docType class
##' @slot result enrichment analysis results
##' @slot pvalueCutoff cutoff pvalue
##' @slot pAdjustMethod pvalue adjust method
##' @slot padjCutoff pvalue adjust cutoff value
##' @slot genenumber number of input genes
##' @slot organism organism used
##' @slot ontology biological ontology
##' @slot gene Gene IDs
##' @slot keytype Gene ID type
##' @slot sep character string used to separate the genes when concatenating
##' @exportClass GSEAResult
##' @author Kai Guo
##' @keywords classes
setClass("GSEAResult",
         representation   = representation(
           result          = "data.frame",
           pvalueCutoff   = "numeric",
           pAdjustMethod  = "character",
           padjCutoff   = "numeric",
           genenumber    = "numeric",
           organism       = "character",
           gene           = "character",
           input = 'numeric',
           ontology = "character",
           keytype        = "character",
           sep = "character"
         )
)
##' Class "Annot"
##' This class represents the Annotation information
##' @name Annot-class
##' @aliases Annot-class
##'   summary, Annot-method
##' @docType class
##' @slot species the species of the annotation file
##' @slot anntype the type of the annotation file
##' @slot keytype Gene ID type
##' @slot annot Annotation information data.frame
##' @exportClass Annot
##' @author Kai Guo
##' @keywords classes
setClass("Annot",
         representation = representation(
           species="character",
           anntype="character",
           keytype="character",
           annot="data.frame"
         ))

