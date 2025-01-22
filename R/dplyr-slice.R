##' @importFrom dplyr slice
##' @importFrom magrittr %<>%
##' @method slice richResult
##' @export
slice.richResult <- function(.data, ..., .preserve = FALSE) {
    dots <- quos(...)
    .data@result %<>% slice(!!!dots, .preserve = .preserve)
    return(.data)
}

##' @method slice GSEAResult
##' @export
slice.GSEAResult <- slice.richResult

##' @importFrom dplyr slice
##' @importFrom magrittr %<>%
##' @method slice richResult
##' @export
slice.Annot <- function(.data, ..., .preserve = FALSE) {
  dots <- quos(...)
  .data@annot %<>% slice(!!!dots, .preserve = .preserve)
  return(.data)
}
