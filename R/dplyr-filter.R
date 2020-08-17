##' @importFrom dplyr filter
##' @importFrom magrittr %<>%
##' @method filter richResult
##' @export
filter.richResult <- function(.data, ..., .preserve = FALSE) {
    dots <- quos(...)
    .data@result %<>% filter(!!!dots, .preserve = .preserve)
    return(.data)
}


##' @method filter GSEAResult
##' @export
filter.GSEAResult <- filter.richResult



