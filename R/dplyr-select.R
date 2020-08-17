##' @method select richResult
##' @importFrom magrittr %<>%
##' @export
select.richResult<- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% dplyr::select(!!!dots,)
    return(.data)
}

##' @method select GSEAResult
##' @export
select.GSEAResult <- select.richResult



