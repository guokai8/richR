##' @method select richResult
##' @importFrom magrittr %<>%
##' @importFrom dplyr select
##' @export
select.richResult<- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% select(!!!dots,)
    return(.data)
}

##' @method select GSEAResult
##' @export
select.GSEAResult <- select.richResult


##' @method select richResult
##' @importFrom magrittr %<>%
##' @importFrom dplyr select
##' @export
select.Annot<- function(.data, ...) {
  dots <- quos(...)
  .data@annot %<>% select(!!!dots,)
  return(.data)
}
