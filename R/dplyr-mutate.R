##' @method mutate richResult
##' @importFrom magrittr %<>%
##' @importFrom dplyr mutate
##' @importFrom rlang quos
##' @export
##' @author Kai Guo
mutate.richResult <- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% mutate(!!!dots)
    return(.data)
}

##' @method mutate GSEAResult
##' @export
mutate.GSEAResult <- mutate.richResult

##' @importFrom magrittr %<>%
##' @importFrom dplyr mutate
##' @importFrom rlang quos
##' @export
mutate.Annot<- function(.data, ...) {
  dots <- quos(...)
  .data@annot %<>% mutate(!!!dots)
  return(.data)
}
