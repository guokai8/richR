##' @importFrom dplyr summarise
##' @importFrom magrittr %<>%
##' @method summarise richResult
##' @export
summarise.richResult <- function(.data, ...) {
    dots <- quos(...)
    .data@result %>% summarise(!!!dots)
}


##' @method summarise GSEAResult
##' @export
summarise.GSEAResult <- summarise.richResult

##' @importFrom dplyr summarise
##' @importFrom magrittr %<>%
##' @method summarise richResult
##' @export
summarise.Annot <- function(.data, ...) {
  dots <- quos(...)
  .data@annot %>% summarise(!!!dots)
}
