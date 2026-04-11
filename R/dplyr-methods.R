## ----------------------------------------------------------------
##  Re-exported dplyr / magrittr functions
## ----------------------------------------------------------------

#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom magrittr %<>%
#' @export
magrittr::`%<>%`

#' @importFrom dplyr filter
#' @export
dplyr::filter

#' @importFrom dplyr select
#' @export
dplyr::select

#' @importFrom dplyr mutate
#' @export
dplyr::mutate

#' @importFrom dplyr arrange
#' @export
dplyr::arrange

#' @importFrom dplyr rename
#' @export
dplyr::rename

#' @importFrom dplyr slice
#' @export
dplyr::slice

#' @importFrom dplyr summarise
#' @export
dplyr::summarise

#' @importFrom dplyr group_by
#' @export
dplyr::group_by

#' @importFrom dplyr n
#' @export
dplyr::n

## ----------------------------------------------------------------
##  S3 methods: filter
## ----------------------------------------------------------------

#' @importFrom magrittr %<>%
#' @method filter richResult
#' @export
filter.richResult <- function(.data, ..., .preserve = FALSE) {
  dots <- quos(...)
  .data@result %<>% filter(!!!dots, .preserve = .preserve)
  return(.data)
}

#' @method filter GSEAResult
#' @export
filter.GSEAResult <- filter.richResult

#' @importFrom magrittr %<>%
#' @method filter Annot
#' @export
filter.Annot <- function(.data, ..., .preserve = FALSE) {
  dots <- quos(...)
  .data@annot %<>% filter(!!!dots, .preserve = .preserve)
  return(.data)
}

## ----------------------------------------------------------------
##  S3 methods: select
## ----------------------------------------------------------------

#' @importFrom magrittr %<>%
#' @importFrom dplyr select
#' @method select richResult
#' @export
select.richResult <- function(.data, ...) {
  dots <- quos(...)
  .data@result %<>% select(!!!dots, )
  return(.data)
}

#' @method select GSEAResult
#' @export
select.GSEAResult <- select.richResult

#' @importFrom magrittr %<>%
#' @importFrom dplyr select
#' @method select Annot
#' @export
select.Annot <- function(.data, ...) {
  dots <- quos(...)
  .data@annot %<>% select(!!!dots, )
  return(.data)
}

## ----------------------------------------------------------------
##  S3 methods: mutate
## ----------------------------------------------------------------

#' @importFrom magrittr %<>%
#' @method mutate richResult
#' @export
mutate.richResult <- function(.data, ...) {
  dots <- quos(...)
  .data@result %<>% mutate(!!!dots)
  return(.data)
}

#' @method mutate GSEAResult
#' @export
mutate.GSEAResult <- mutate.richResult

#' @importFrom magrittr %<>%
#' @method mutate Annot
#' @export
mutate.Annot <- function(.data, ...) {
  dots <- quos(...)
  .data@annot %<>% mutate(!!!dots)
  return(.data)
}

## ----------------------------------------------------------------
##  S3 methods: arrange
## ----------------------------------------------------------------

#' @importFrom magrittr %<>%
#' @method arrange richResult
#' @export
arrange.richResult <- function(.data, ...) {
  dots <- quos(...)
  .data@result %<>% arrange(!!!dots, )
  return(.data)
}

#' @method arrange GSEAResult
#' @export
arrange.GSEAResult <- arrange.richResult

#' @importFrom magrittr %<>%
#' @method arrange Annot
#' @export
arrange.Annot <- function(.data, ...) {
  dots <- quos(...)
  .data@annot %<>% arrange(!!!dots, )
  return(.data)
}

## ----------------------------------------------------------------
##  S3 methods: rename
## ----------------------------------------------------------------

#' @importFrom magrittr %<>%
#' @method rename richResult
#' @export
rename.richResult <- function(.data, ...) {
  dots <- quos(...)
  .data@result %<>% rename(!!!dots, )
  return(.data)
}

#' @method rename GSEAResult
#' @export
rename.GSEAResult <- rename.richResult

#' @importFrom magrittr %<>%
#' @method rename Annot
#' @export
rename.Annot <- function(.data, ...) {
  dots <- quos(...)
  .data@annot %<>% rename(!!!dots, )
  return(.data)
}

## ----------------------------------------------------------------
##  S3 methods: slice
## ----------------------------------------------------------------

#' @importFrom magrittr %<>%
#' @method slice richResult
#' @export
slice.richResult <- function(.data, ..., .preserve = FALSE) {
  dots <- quos(...)
  .data@result %<>% slice(!!!dots, .preserve = .preserve)
  return(.data)
}

#' @method slice GSEAResult
#' @export
slice.GSEAResult <- slice.richResult

#' @importFrom magrittr %<>%
#' @method slice Annot
#' @export
slice.Annot <- function(.data, ..., .preserve = FALSE) {
  dots <- quos(...)
  .data@annot %<>% slice(!!!dots, .preserve = .preserve)
  return(.data)
}

## ----------------------------------------------------------------
##  S3 methods: summarise
## ----------------------------------------------------------------

#' @importFrom magrittr %<>%
#' @method summarise richResult
#' @export
summarise.richResult <- function(.data, ...) {
  dots <- quos(...)
  .data@result %>% summarise(!!!dots)
}

#' @method summarise GSEAResult
#' @export
summarise.GSEAResult <- summarise.richResult

#' @importFrom magrittr %<>%
#' @method summarise Annot
#' @export
summarise.Annot <- function(.data, ...) {
  dots <- quos(...)
  .data@annot %>% summarise(!!!dots)
}

## ----------------------------------------------------------------
##  S3 methods: group_by
## ----------------------------------------------------------------

#' @importFrom magrittr %<>%
#' @method group_by richResult
#' @export
group_by.richResult <- function(.data, ..., add = FALSE, .drop = FALSE) {
  dots <- quos(...)
  .data@result %<>% group_by(!!!dots, add = add, .drop = .drop)
  return(.data)
}

#' @method group_by GSEAResult
#' @export
group_by.GSEAResult <- group_by.richResult

#' @importFrom dplyr group_by
#' @importFrom magrittr %<>%
#' @method group_by Annot
#' @export
group_by.Annot <- function(.data, ..., add = FALSE, .drop = FALSE) {
  dots <- quos(...)
  .data@annot %<>% group_by(!!!dots, add = add, .drop = .drop)
  return(.data)
}
