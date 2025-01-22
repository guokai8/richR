##' @importFrom dplyr rename
##' @importFrom magrittr %<>%
##' @method rename richResult
##' @export
rename.richResult <- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% rename(!!!dots,)
    return(.data)
}

##' @method rename GSEAResult
##' @export
rename.GSEAResult <- rename.richResult
##' @importFrom dplyr rename
##' @importFrom magrittr %<>%
##' @method rename Annot
##' @export
rename.Annot<- function(.data, ...) {
  dots <- quos(...)
  .data@annot %<>% rename(!!!dots,)
  return(.data)
}
