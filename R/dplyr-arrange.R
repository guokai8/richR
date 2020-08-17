##' @importFrom dplyr arrange
##' @importFrom magrittr %<>%
##' @method arrange richResult
##' @export
arrange.richResult<- function(.data, ...) {
    dots <- quos(...)
    .data@result%<>%arrange(!!!dots,)
    return(.data)
}

##' @method arrange GSEAResult
##' @export
arrange.GSEAResult <- arrange.richResult


