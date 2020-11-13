#' @importFrom dplyr group_by
##' @importFrom magrittr %<>%
##' @method group_by richResult
##' @export
group_by.richResult<- function(.data, ..., add = FALSE, .drop = FALSE) {
    dots <- quos(...)
    .data@result %<>% group_by(!!!dots, add = add, .drop = .drop)
    return(.data)
}


##' @method group_by GSEAResult
##' @export
group_by.GSEAResult<- group_by.richResult

#' @importFrom dplyr group_by
##' @importFrom magrittr %<>%
##' @method group_by Annot
group_by.Annot<- function(.data, ..., add = FALSE, .drop = FALSE) {
  dots <- quos(...)
  .data@annot %<>% group_by(!!!dots, add = add, .drop = .drop)
  return(.data)
}
