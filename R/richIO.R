#' Write an enrichment result to a folder of tables
#'
#' Stores a \code{richResult} or \code{GSEAResult} as a folder holding
#' \code{manifest.json} and one UTF-8 tab-separated table per data frame of the
#' object, so that \code{\link{readRich}} can rebuild it exactly.
#'
#' The folder holds \code{result.tsv} (\code{result(x)}), \code{tested.tsv}
#' (the tested table, when one is recorded; see \code{\link{testedTerms}}),
#' \code{detail.tsv} (\code{richResult}) or \code{input.tsv} (\code{GSEAResult}:
#' the ranked input as gene and score, in the order of the \code{input} slot).
#' Each table has a header row and one line per row. In character, factor and
#' list cells a backslash is written \code{\\\\}, TAB \code{\\t}, newline
#' \code{\\n} and carriage return \code{\\r}; a missing value of any class is
#' written \code{\\N}; doubles are written with 17 significant digits
#' (\code{NaN}, \code{Inf}, \code{-Inf} spelled out), which read back
#' bit-identical; a list
#' cell (a list of character vectors, such as a GSEA
#' \code{leadingEdge} with \code{table = FALSE}) is written as a JSON array. A
#' table with zero columns is recorded in the manifest only, without a file.
#'
#' \code{manifest.json} records the format name \code{richR-result}, the format
#' version \code{2}, the class, every other slot, the
#' \code{\link{enrichmentInfo}} metadata (or its absence), whether a tested
#' table is recorded, and for each table its file name, column names, column
#' classes, row names, row count and SHA-256 checksum of each table file.
#'
#' Supported column classes are character, integer, double, logical, factor and
#' lists of character vectors. Every check runs before any file is written:
#' \code{writeRich()} stops on an unsupported column, on a table that is not a
#' plain \code{data.frame} (for example the \code{grouped_df} that
#' \code{group_by()} leaves in \code{@result}; convert it with
#' \code{x@result <- as.data.frame(x@result)}), on an \code{@result} without
#' its identifier column (\code{Annot} or \code{pathway}), or an ORA result
#' without \code{Term}, and
#' on an existing non-empty \code{path}. With \code{overwrite = TRUE} only a
#' folder previously written by \code{writeRich()} is replaced; its table files
#' and manifest are removed, nothing else.
#'
#' @param x a \code{richResult} or \code{GSEAResult} object
#' @param path folder to write; created when it does not exist
#' @param overwrite replace a folder previously written by \code{writeRich()}
#' @return \code{path}, invisibly
#' @seealso \code{\link{readRich}}, \code{\link{testedTerms}}
#' @importFrom jsonlite toJSON fromJSON serializeJSON unserializeJSON
#' @export
writeRich <- function(x, path, overwrite = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("writeRich(): `path` must be a single folder path.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("writeRich(): `overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  supplied_path <- path
  path <- path.expand(path)
  ## Everything is checked and encoded in memory before the first file
  ## operation, so a refused call leaves the disk as it was.
  files <- .rich_io_encode(x)
  replace <- .rich_io_check_target(path, overwrite)
  if (replace) {
    old <- file.path(path, .rich_io_known_files)
    if (file.access(path, 2L) != 0L) {
      stop("writeRich(): folder '", path, "' is not writable.", call. = FALSE)
    }
    if (unlink(old[file.exists(old)], expand = FALSE) != 0L) {
      stop("writeRich(): could not remove the old files in folder '", path, "'.",
           call. = FALSE)
    }
  }
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  ## The manifest goes last: a folder without one is never read as a result.
  for (f in names(files)) .rich_io_write_lines(files[[f]], file.path(path, f))
  invisible(supplied_path)
}

## The file names writeRich() writes; the only files overwrite = TRUE removes.
.rich_io_known_files <- c("manifest.json", "result.tsv", "tested.tsv",
                          "detail.tsv", "input.tsv")

.rich_io_format <- "richR-result"
.rich_io_format_version <- 2L

## Non-table slots of a result class (stored typed in the manifest).
.rich_io_table_slots <- c("result", "detail", "tested", "input")
.rich_io_scalar_slots <- function(cls) {
  setdiff(methods::slotNames(cls), .rich_io_table_slots)
}

.rich_io_id_column <- function(cls) {
  if (identical(cls, "GSEAResult")) "pathway" else "Annot"
}

## Validate x and return the folder content as a named list of text lines,
## manifest.json last.
.rich_io_encode <- function(x) {
  cls <- class(x)[1L]
  if (!cls %in% c("richResult", "GSEAResult")) {
    stop("writeRich() supports richResult or GSEAResult objects, not ",
         cls, ".", call. = FALSE)
  }
  tables <- list(result = x@result)
  if (.has_tested(x)) tables$tested <- x@tested
  if (identical(cls, "richResult")) {
    tables$detail <- x@detail
  } else {
    tables$input <- data.frame(
      gene = if (is.null(names(x@input))) rep(NA_character_, length(x@input))
             else names(x@input),
      score = unname(x@input), stringsAsFactors = FALSE)
  }
  for (s in names(tables)) {
    if (!identical(class(tables[[s]]), "data.frame")) {
      stop(sprintf(paste("writeRich(): @%s is a %s, not a plain data.frame;",
                         "convert it first, e.g. x@%s <- as.data.frame(x@%s)."),
                   s, paste(class(tables[[s]]), collapse = "/"), s, s),
           call. = FALSE)
    }
  }
  needed <- c(.rich_io_id_column(cls), if (identical(cls, "richResult")) "Term")
  lacking <- setdiff(needed, names(x@result))
  if (length(lacking)) {
    stop(sprintf("writeRich(): @result lacks the column(s) %s, which a %s needs.",
                 paste0("'", lacking, "'", collapse = ", "), cls), call. = FALSE)
  }
  ## Column classes of every table, checked before anything is encoded.
  classes <- lapply(names(tables), function(s) {
    tab <- tables[[s]]
    vapply(seq_along(tab), function(j) {
      k <- .rich_io_column_class(tab[[j]])
      if (is.na(k)) {
        stop(sprintf(paste("writeRich(): column '%s' of @%s has class %s;",
                           "supported: character, integer, double, logical,",
                           "factor and lists of character vectors."),
                     names(tab)[j], s, paste(class(tab[[j]]), collapse = "/")),
             call. = FALSE)
      }
      k
    }, character(1))
  })
  names(classes) <- names(tables)

  files <- list()
  entries <- list()
  for (s in names(tables)) {
    tab <- tables[[s]]
    file <- if (ncol(tab)) paste0(s, ".tsv") else NULL
    entry <- list(slot = s, file = file, nrow = nrow(tab), ncol = ncol(tab),
                  columns = .rich_io_column_entries(tab, classes[[s]]),
                  row_names = .rich_io_row_names(tab))
    if (identical(s, "input")) entry$names <- !is.null(names(x@input))
    if (!is.null(file)) {
      lines <- .rich_io_table_lines(tab, classes[[s]])
      files[[file]] <- lines
      bytes <- charToRaw(enc2utf8(paste0(paste(lines, collapse = "\n"), "\n")))
      entry$sha256 <- digest::digest(bytes, algo = "sha256", serialize = FALSE)
    }
    entries[[length(entries) + 1L]] <- entry
  }
  metadata <- attr(x, "richR_metadata", exact = TRUE)
  slots <- lapply(.rich_io_scalar_slots(cls), function(s) methods::slot(x, s))
  names(slots) <- .rich_io_scalar_slots(cls)
  manifest <- list(
    format = .rich_io_format,
    format_version = .rich_io_format_version,
    class = cls,
    tested_recorded = .has_tested(x),
    metadata_recorded = !is.null(metadata),
    metadata = if (is.null(metadata)) NULL else .rich_io_typed_json(metadata),
    slots = .rich_io_typed_json(slots),
    tables = entries
  )
  files[["manifest.json"]] <- as.character(jsonlite::toJSON(
    manifest, auto_unbox = TRUE, null = "null", digits = I(17), pretty = TRUE,
    json_verbatim = TRUE))
  files
}

## The class name the manifest records for a column, or NA when unsupported.
.rich_io_column_class <- function(v) {
  cls <- class(v)
  if (identical(cls, "factor")) {
    return(if (anyNA(levels(v))) NA_character_ else "factor")
  }
  if (identical(cls, "list")) {
    ok <- all(vapply(v, function(e) is.character(e) && is.null(attributes(e)),
                     logical(1)))
    return(if (ok) "list" else NA_character_)
  }
  switch(paste(cls, collapse = "/"),
         character = "character", integer = "integer", numeric = "double",
         logical = "logical", NA_character_)
}

.rich_io_column_entries <- function(tab, classes) {
  lapply(seq_along(tab), function(j) {
    entry <- list(name = names(tab)[j], class = classes[[j]])
    if (identical(classes[[j]], "factor")) entry$levels <- levels(tab[[j]])
    entry
  })
}

## The internal row.names attribute: automatic (compact), integer or character.
.rich_io_row_names <- function(tab) {
  rn <- .row_names_info(tab, type = 0L)
  if (is.character(rn)) return(list(type = "character", values = rn))
  if (length(rn) == 2L && is.na(rn[1L])) {
    return(list(type = "automatic", n = rn[2L]))
  }
  list(type = "integer", values = rn)
}

## Typed JSON (jsonlite::serializeJSON, 17 significant digits): types, NULLs,
## empty vectors and names survive exactly (plan.md B7). Embedded verbatim in
## the manifest.
.rich_io_typed_json <- function(value) {
  structure(jsonlite::serializeJSON(value, digits = I(17)), class = "json")
}

## Escape a character vector for a TSV cell: backslash, TAB, newline and CR
## become \\, \t, \n, \r; a missing value becomes the reserved token \N.
.rich_io_escape <- function(v) {
  out <- enc2utf8(as.character(v))
  hit <- which(!is.na(out) & grepl("[\\\\\t\n\r]", out, perl = TRUE))
  if (length(hit)) {
    y <- out[hit]
    y <- gsub("\\", "\\\\", y, fixed = TRUE)
    y <- gsub("\t", "\\t", y, fixed = TRUE)
    y <- gsub("\n", "\\n", y, fixed = TRUE)
    out[hit] <- gsub("\r", "\\r", y, fixed = TRUE)
  }
  out[is.na(v)] <- "\\N"
  out
}

.rich_io_cells <- function(v, class) {
  switch(class,
    character = .rich_io_escape(v),
    factor = .rich_io_escape(as.character(v)),
    integer = ifelse(is.na(v), "\\N", as.character(v)),
    logical = ifelse(is.na(v), "\\N", ifelse(v, "TRUE", "FALSE")),
    double = {
      out <- sprintf("%.17g", v)
      out[is.na(v) & !is.nan(v)] <- "\\N"
      out
    },
    list = .rich_io_escape(vapply(v, function(e) {
      as.character(jsonlite::toJSON(e, na = "null"))
    }, character(1)))
  )
}

.rich_io_table_lines <- function(tab, classes) {
  header <- paste(.rich_io_escape(names(tab)), collapse = "\t")
  if (!nrow(tab)) return(header)
  cells <- lapply(seq_along(tab), function(j) .rich_io_cells(tab[[j]], classes[[j]]))
  c(header, do.call(paste, c(unname(cells), sep = "\t")))
}

.rich_io_write_lines <- function(lines, file) {
  con <- file(file, open = "wb")
  on.exit(close(con))
  writeLines(enc2utf8(lines), con, sep = "\n", useBytes = TRUE)
}

## Decide what to do with `path`: FALSE = write into a new or empty folder,
## TRUE = replace a previous writeRich() folder; anything else stops.
.rich_io_check_target <- function(path, overwrite) {
  if (!file.exists(path)) return(FALSE)
  if (!dir.exists(path)) {
    stop(sprintf("writeRich(): '%s' exists and is not a folder.", path),
         call. = FALSE)
  }
  if (!length(list.files(path, all.files = TRUE, no.. = TRUE))) return(FALSE)
  if (!overwrite) {
    stop(sprintf(paste("writeRich(): folder '%s' is not empty; use",
                       "overwrite = TRUE to replace a folder written by",
                       "writeRich()."), path), call. = FALSE)
  }
  if (!.rich_io_is_folder(path)) {
    stop(sprintf(paste("writeRich(): folder '%s' was not written by writeRich()",
                       "(no manifest.json naming format '%s'); it is left",
                       "unchanged."), path, .rich_io_format), call. = FALSE)
  }
  TRUE
}

.rich_io_is_folder <- function(path) {
  f <- file.path(path, "manifest.json")
  if (!file.exists(f) || dir.exists(f)) return(FALSE)
  m <- tryCatch(.rich_io_read_json(f), error = function(e) NULL)
  is.list(m) && identical(m$format, .rich_io_format)
}

.rich_io_read_json <- function(file) {
  txt <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  ## fromJSON() follows a path/URL when its text is not valid JSON.
  if (!isTRUE(jsonlite::validate(txt))) {
    stop("the file is not a JSON document", call. = FALSE)
  }
  jsonlite::fromJSON(txt, simplifyVector = FALSE)
}

#' Read an enrichment result written by writeRich()
#'
#' Format 2 verifies SHA-256 table checksums before decoding. Legacy format 1
#' is readable with a warning that checksums are unavailable.
#'
#' Rebuilds the \code{richResult} or \code{GSEAResult} that
#' \code{\link{writeRich}} stored in \code{path}. Every value comes back exact,
#' doubles included; the tested table is unavailable exactly when it was
#' unavailable in the stored object; and \code{\link{enrichmentInfo}} returns
#' the stored record unchanged, so reading never alters, upgrades or re-dates
#' the analysis record.
#'
#' Version policy: the folder format version must be \code{1} or \code{2}, and when the
#' folder records metadata, its \code{schema_version} must be \code{1}
#' (results saved without integrated evidence support) or \code{2}; a missing, non-integer,
#' older or newer version is an error naming the found and the supported
#' versions. A folder written from an object without metadata (such as
#' \code{as(<data.frame>, "richResult")}) is read back without metadata, as
#' it was stored. \code{readRich()} also stops,
#' naming the folder, file, column and first offending row, when
#' \code{manifest.json} is missing or names another format or a class other
#' than \code{richResult} or \code{GSEAResult}, when a table file the manifest
#' lists is missing or a table file it does not list is present, when a
#' table's header or row count differs from the manifest, on a line with more
#' or fewer fields than the manifest's columns, on an unknown escape sequence,
#' and on a cell that does not parse to its column's class (a double only in
#' the form \code{writeRich()} writes it). It never returns a partially read or
#' silently coerced object.
#'
#' The metadata and slots stored in \code{manifest.json} are read as plain
#' data only: \code{NULL}, logical, integer, double, character and lists, with
#' names. A manifest holding a function, call, environment, namespace, raw or
#' serialized bytes, an S4 object or a classed value is refused with an error
#' naming the field, so reading a folder never rebuilds code or loads a
#' package.
#'
#' @param path folder written by \code{writeRich()}
#' @return the stored \code{richResult} or \code{GSEAResult}
#' @seealso \code{\link{writeRich}}
#' @export
readRich <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !dir.exists(path)) {
    stop("readRich(): `path` must be a folder written by writeRich().",
         call. = FALSE)
  }
  mf <- file.path(path, "manifest.json")
  if (!file.exists(mf)) {
    stop(sprintf("readRich(): no manifest.json in folder '%s'.", path),
         call. = FALSE)
  }
  m <- tryCatch(.rich_io_read_json(mf), error = function(e) {
    stop(sprintf("readRich(): manifest.json of '%s' is not valid JSON: %s",
                 path, conditionMessage(e)), call. = FALSE)
  })
  where <- sprintf("manifest.json of '%s'", path)
  if (!is.list(m) || !identical(m$format, .rich_io_format)) {
    stop(sprintf("readRich(): %s records format %s, not '%s'.", where,
                 .rich_io_found(if (is.list(m)) m$format), .rich_io_format),
         call. = FALSE)
  }
  .rich_io_check_version(m$format_version, 1:.rich_io_format_version,
                         "format version", where)
  cls <- m$class
  if (!(is.character(cls) && length(cls) == 1L &&
        cls %in% c("richResult", "GSEAResult"))) {
    stop(sprintf("readRich(): class %s in %s is not supported (supported: richResult, GSEAResult).",
                 .rich_io_found(cls), where), call. = FALSE)
  }
  bad <- function(what) {
    stop(sprintf("readRich(): %s is inconsistent: %s.", where, what), call. = FALSE)
  }
  for (flag in c("tested_recorded", "metadata_recorded")) {
    if (!(is.logical(m[[flag]]) && length(m[[flag]]) == 1L && !is.na(m[[flag]]))) {
      bad(paste(flag, "is not TRUE or FALSE"))
    }
  }

  metadata <- NULL
  if (m$metadata_recorded) {
    metadata <- .rich_io_untyped(m$metadata, "metadata", where,
                                 function() bad("the metadata cannot be decoded"))
    if (!is.list(metadata)) bad("the metadata is not a list")
    .rich_io_check_version(metadata$schema_version, 1:2,
                           "metadata schema_version", where)
  } else if (!is.null(m$metadata)) {
    bad("metadata is present but metadata_recorded is FALSE")
  }
  slots <- .rich_io_untyped(m$slots, "slots", where,
                            function() bad("the slots cannot be decoded"))
  expected_slots <- .rich_io_scalar_slots(cls)
  if (!is.list(slots) || !setequal(names(slots), expected_slots) ||
      anyDuplicated(names(slots))) {
    bad(paste("the slots do not match those of", cls))
  }

  ## Table entries: the tables of the class, tested only when recorded.
  entries <- m$tables
  wanted <- c("result", if (m$tested_recorded) "tested",
              if (identical(cls, "richResult")) "detail" else "input")
  listed <- vapply(entries, function(e) {
    if (is.list(e) && is.character(e$slot) && length(e$slot) == 1L) e$slot else NA_character_
  }, character(1))
  if (!is.list(entries) || anyNA(listed) || anyDuplicated(listed) ||
      !setequal(listed, wanted)) {
    bad(sprintf("it lists the tables %s, expected %s",
                paste(listed, collapse = ", "), paste(wanted, collapse = ", ")))
  }
  names(entries) <- listed
  files <- unlist(lapply(entries, function(e) e$file))
  present <- .rich_io_known_files[-1L][file.exists(file.path(path, .rich_io_known_files[-1L]))]
  unlisted <- setdiff(present, files)
  if (length(unlisted)) {
    stop(sprintf("readRich(): %s in folder '%s' is not listed in manifest.json.",
                 unlisted[1L], path), call. = FALSE)
  }
  if (m$format_version == 1L) {
    warning("readRich(): checksums unavailable in legacy format 1; typed validation still applies.",
            call. = FALSE)
  }
  tables <- lapply(entries, .rich_io_read_table, path = path, bad = bad,
                   verify_checksum = m$format_version >= 2L)
  required <- c(.rich_io_id_column(cls), if (identical(cls, "richResult")) "Term")
  for (name in intersect(c("result", "tested"), names(tables))) {
    absent <- setdiff(required, names(tables[[name]]))
    if (length(absent)) bad(paste(name, "lacks required column(s)", paste(absent, collapse = ", ")))
  }

  obj <- methods::new(cls)
  tryCatch({
    for (s in expected_slots) methods::slot(obj, s) <- slots[[s]]
    obj@result <- tables$result
    if (m$tested_recorded) obj@tested <- tables$tested
    if (identical(cls, "richResult")) {
      obj@detail <- tables$detail
    } else {
      input <- tables$input
      if (!identical(names(input), c("gene", "score"))) bad("input has other columns than gene and score")
      v <- input$score
      if (isTRUE(entries$input$names)) names(v) <- input$gene
      obj@input <- v
    }
  }, error = function(e) {
    stop(sprintf("readRich(): cannot rebuild the %s stored in '%s': %s", cls,
                 path, conditionMessage(e)), call. = FALSE)
  })
  if (!is.null(metadata)) attr(obj, "richR_metadata") <- metadata
  obj
}

## A value from a parsed manifest as text for an error message.
.rich_io_found <- function(v) {
  if (is.null(v)) return("(none)")
  v <- unlist(v)
  if (is.character(v)) paste0('"', v, '"', collapse = ", ")
  else paste(format(v), collapse = ", ")
}

## A version must be one whole number within `supported`.
.rich_io_check_version <- function(v, supported, what, where) {
  sup <- paste(supported, collapse = ", ")
  if (is.null(v)) {
    stop(sprintf("readRich(): %s records no %s (supported: %s).", where, what, sup),
         call. = FALSE)
  }
  whole <- is.numeric(v) && length(v) == 1L && is.finite(v) && v == round(v)
  reason <- if (!whole) "is not a whole number"
            else if (v < min(supported)) "is older than the oldest supported"
            else if (v > max(supported)) "is newer than the newest supported"
            else NULL
  if (!is.null(reason)) {
    stop(sprintf("readRich(): %s %s in %s %s (supported: %s).", what,
                 .rich_io_found(v), where, reason, sup), call. = FALSE)
  }
}

## Decode a typed JSON node of the manifest (see .rich_io_typed_json()). The
## node is checked to hold plain data before jsonlite::unserializeJSON() sees
## it, which would otherwise rebuild functions, calls, environments, S4 objects
## and serialized bytes from a folder the user did not write, and load the
## namespaces such a folder names.
.rich_io_untyped <- function(node, field, where, fail) {
  if (is.null(node)) fail()
  .rich_io_check_plain(node, field, where)
  tryCatch(
    jsonlite::unserializeJSON(jsonlite::toJSON(node, auto_unbox = TRUE,
                                               null = "null", digits = I(17))),
    error = function(e) fail())
}

## What a typed JSON node of the manifest may hold: the node types and
## attributes the metadata and slots of richR results are built from (lists of
## NULL, logical, integer, double and character values; names on them; and
## the package attribute of the class name a GSEA run records for BPPARAM).
.rich_io_plain_types <- c("NULL", "logical", "integer", "double", "character",
                          "list")
.rich_io_plain_attributes <- c("names", "package")

## Stop, naming the manifest field, on the first node that is not plain data.
.rich_io_check_plain <- function(node, at, where) {
  refuse <- function(what) {
    stop(sprintf(paste("readRich(): %s holds %s at %s; readRich() decodes only",
                       "plain data (NULL, logical, integer, double, character",
                       "and lists, with names)."), where, what, at),
         call. = FALSE)
  }
  type <- if (is.list(node)) node$type
  if (!(is.character(type) && length(type) == 1L)) refuse("a value without a type")
  if (!type %in% .rich_io_plain_types) refuse(sprintf("a value of type '%s'", type))
  if (identical(type, "integer")) {
    values <- node$value
    valid <- vapply(values, function(v) {
      identical(v, "NA") ||
        (is.numeric(v) && length(v) == 1L && is.finite(v) &&
         v == trunc(v) && abs(v) <= .Machine$integer.max)
    }, logical(1))
    if (!all(valid)) refuse("an integer value that cannot be represented exactly")
  }
  attrs <- node$attributes
  if (!is.null(attrs) && !is.list(attrs)) refuse("malformed attributes")
  if (length(attrs)) {
    keys <- names(attrs)
    if (is.null(keys) || anyNA(keys) || !all(nzchar(keys)) || anyDuplicated(keys)) {
      refuse("attributes without unique names")
    }
  }
  for (i in seq_along(attrs)) {
    a <- names(attrs)[i]
    if (!a %in% .rich_io_plain_attributes) refuse(sprintf("attribute '%s'", a))
    an <- attrs[[i]]
    if (!(is.list(an) && identical(an$type, "character") && !length(an$attributes))) {
      refuse(sprintf("attribute '%s' that is not text", a))
    }
  }
  if (identical(type, "list")) {
    if (!is.null(node$value) && !is.list(node$value)) refuse("a list without elements")
    labels <- vapply(attrs$names$value, function(v) {
      if (is.character(v) && length(v) == 1L) v else ""
    }, character(1))
    for (i in seq_along(node$value)) {
      label <- if (i <= length(labels) && nzchar(labels[i])) {
        paste0(at, "$", labels[i])
      } else {
        sprintf("%s[[%d]]", at, i)
      }
      .rich_io_check_plain(node$value[[i]], label, where)
    }
  }
  invisible(TRUE)
}

## Read one table as the manifest describes it; every check runs before the
## data.frame is built.
.rich_io_read_table <- function(entry, path, bad, verify_checksum = FALSE) {
  slot <- entry$slot
  columns <- entry$columns
  n <- entry$nrow
  if (!(is.numeric(n) && length(n) == 1L && is.finite(n) &&
        n >= 0 && n <= .Machine$integer.max && n == trunc(n))) {
    bad(sprintf("the row count of %s is not an exactly representable whole number", slot))
  }
  nc <- entry$ncol
  if (!(is.numeric(nc) && length(nc) == 1L && is.finite(nc) &&
        nc >= 0 && nc <= .Machine$integer.max && nc == trunc(nc))) {
    bad(sprintf("the column count of %s is not an exactly representable whole number", slot))
  }
  if (!is.list(columns) || !identical(length(columns), as.integer(nc))) {
    bad(sprintf("the columns of %s do not match its column count", slot))
  }
  col_names <- vapply(columns, function(cl) as.character(cl$name)[1L], character(1))
  col_classes <- vapply(columns, function(cl) as.character(cl$class)[1L], character(1))
  known <- c("character", "integer", "double", "logical", "factor", "list")
  if (anyNA(col_names) || anyNA(col_classes) || !all(col_classes %in% known)) {
    bad(sprintf("the column records of %s are incomplete or name an unknown class", slot))
  }
  row_names <- .rich_io_restore_row_names(entry$row_names, n, slot, bad)
  if (!length(columns)) {
    if (!is.null(entry$file)) bad(sprintf("the zero-column table %s names a file", slot))
    return(structure(list(), names = character(0), row.names = row_names,
                     class = "data.frame"))
  }
  file <- entry$file
  if (!identical(file, paste0(slot, ".tsv"))) {
    bad(sprintf("the file of %s is not %s.tsv", slot, slot))
  }
  f <- file.path(path, file)
  if (!file.exists(f)) {
    stop(sprintf("readRich(): table file %s listed in manifest.json of '%s' is missing.",
                 file, path), call. = FALSE)
  }
  if (isTRUE(verify_checksum)) {
    hash <- entry$sha256
    if (!(is.character(hash) && length(hash) == 1L && !is.na(hash) &&
          grepl("^[0-9a-f]{64}$", hash)) ||
        !identical(digest::digest(f, algo = "sha256", file = TRUE), hash)) {
      bad(paste("checksum mismatch or missing SHA-256 for", file))
    }
  }
  cannot_read <- function(e) {
    stop(sprintf("readRich(): cannot read %s in '%s': %s", file, path,
                 conditionMessage(e)), call. = FALSE)
  }
  fail <- function(column, row, what) {
    stop(sprintf("readRich(): %s in '%s', column '%s', row %d: %s.", file, path,
                 column, row, what), call. = FALSE)
  }
  ## The header and the field count of every line are checked before the table
  ## is parsed, so a malformed line is reported by its own row. Cells never
  ## hold a raw TAB (it is escaped), so splitting on TAB is exact; an empty
  ## line is one empty field.
  first <- tryCatch(readLines(f, n = 1L, warn = FALSE, encoding = "UTF-8"),
                    error = cannot_read)
  header <- .rich_io_unescape(strsplit(paste0(first, "\t"), "\t", fixed = TRUE)[[1L]],
                              function(row, what) fail("(header)", row, what))
  if (!identical(header, col_names)) {
    common <- seq_len(min(length(header), length(col_names)))
    k <- which(header[common] != col_names[common])[1L]
    if (is.na(k)) k <- length(common) + 1L
    stop(sprintf(paste("readRich(): the header of %s in '%s' does not match",
                       "manifest.json: %d columns found, %d recorded; column %d is",
                       "%s, expected %s."),
                 file, path, length(header), length(col_names), k,
                 if (k <= length(header)) paste0("'", header[k], "'") else "absent",
                 if (k <= length(col_names)) paste0("'", col_names[k], "'") else "none"),
         call. = FALSE)
  }
  counts <- tryCatch(
    utils::count.fields(f, sep = "\t", quote = "", comment.char = "",
                        blank.lines.skip = FALSE),
    error = cannot_read)
  counts[counts == 0L] <- 1L
  wrong <- which(counts[-1L] != length(col_names))[1L]
  if (!is.na(wrong)) {
    stop(sprintf(paste("readRich(): %s in '%s', row %d: %d fields found, %d",
                       "expected (the columns manifest.json records)."),
                 file, path, wrong, counts[wrong + 1L], length(col_names)),
         call. = FALSE)
  }
  raw <- tryCatch(
    utils::read.delim(f, header = FALSE, quote = "", na.strings = character(0),
                      blank.lines.skip = FALSE, comment.char = "", fill = FALSE,
                      colClasses = "character", encoding = "UTF-8",
                      check.names = FALSE),
    error = cannot_read)
  body <- raw[-1L, , drop = FALSE]
  if (nrow(body) != n) {
    stop(sprintf("readRich(): %s in '%s' has %d data rows; manifest.json records %d.",
                 file, path, nrow(body), as.integer(n)), call. = FALSE)
  }
  values <- lapply(seq_along(col_names), function(j) {
    cells <- body[[j]]
    missing <- cells == "\\N"
    v <- .rich_io_unescape(cells, function(row, what) fail(col_names[j], row, what),
                           skip = missing)
    .rich_io_parse(v, missing, col_classes[j], columns[[j]],
                   function(row) {
                     fail(col_names[j], row, sprintf("'%s' does not parse as %s",
                                                    cells[row], col_classes[j]))
                   })
  })
  structure(values, names = col_names, row.names = row_names, class = "data.frame")
}

## Undo .rich_io_escape(): \\, \t, \n, \r; any other backslash sequence is an
## error. Cells flagged in `skip` (the missing token) are left untouched.
.rich_io_unescape <- function(cells, fail, skip = rep(FALSE, length(cells))) {
  idx <- which(!skip & grepl("\\", cells, fixed = TRUE))
  if (!length(idx)) return(cells)
  map <- c("\\\\" = "\\", "\\t" = "\t", "\\n" = "\n", "\\r" = "\r")
  y <- cells[idx]
  m <- gregexpr("\\\\.?", y, perl = TRUE)
  tokens <- regmatches(y, m)
  for (i in seq_along(tokens)) {
    unknown <- setdiff(tokens[[i]], names(map))
    if (length(unknown)) {
      fail(idx[i], sprintf("unknown escape sequence '%s'", unknown[1L]))
    }
  }
  regmatches(y, m) <- lapply(tokens, function(t) unname(map[t]))
  cells[idx] <- y
  cells
}

## Convert decoded cells to the recorded class; `fail(row)` on the first cell
## that does not parse.
.rich_io_parse <- function(v, missing, class, column, fail) {
  check <- function(ok) {
    wrong <- which(!missing & !ok)
    if (length(wrong)) fail(wrong[1L])
  }
  out <- switch(class,
    character = v,
    factor = {
      lev <- as.character(unlist(column$levels))
      check(v %in% lev)
      factor(ifelse(missing, NA_character_, v), levels = lev)
    },
    integer = {
      x <- suppressWarnings(as.integer(v))
      check(grepl("^-?[0-9]+$", v) & !is.na(x))
      x
    },
    double = {
      ## Only what writeRich() writes (sprintf("%.17g"), NaN, Inf, -Inf):
      ## as.numeric() alone would also take spaces, hex and other spellings.
      check(grepl("^-?[0-9]+(\\.[0-9]+)?(e[-+][0-9]+)?$", v, perl = TRUE) |
              v %in% c("NaN", "Inf", "-Inf"))
      special <- v %in% c("NaN", "Inf", "-Inf")
      ## R's decimal conversion loses low bits on platforms without extended
      ## long double (notably Apple Silicon). jsonlite uses the system strtod.
      ## The grammar above is checked first; only a JSON array is ever parsed.
      tokens <- v
      tokens[missing | special] <- "0.0"
      tokens <- sub("^(-?)0+([0-9])", "\\1\\2", tokens)
      whole <- !grepl("[.e]", tokens)
      tokens[whole] <- paste0(tokens[whole], ".0")  # retain double -0, not integer 0
      x <- as.double(jsonlite::fromJSON(paste0("[", paste(tokens, collapse = ","), "]")))
      x[special] <- c(NaN, Inf, -Inf)[match(v[special], c("NaN", "Inf", "-Inf"))]
      nonzero <- grepl("[1-9]", sub("e.*$", "", v))
      check(special | (is.finite(x) & !(x == 0 & nonzero)))
      x
    },
    logical = {
      check(v %in% c("TRUE", "FALSE"))
      v == "TRUE"
    },
    list = {
      check(rep(TRUE, length(v)))
      lapply(seq_along(v), function(i) {
        ## Only a JSON array is parsed (fromJSON() would read other text as a
        ## file path or URL).
        p <- if (missing[i] || !startsWith(v[i], "[") ||
                 !isTRUE(jsonlite::validate(v[i]))) NULL else
          tryCatch(jsonlite::fromJSON(v[i]), error = function(e) NULL)
        if (is.list(p) && !length(p)) return(character(0))
        if (is.character(p) && is.null(dim(p))) return(as.character(p))
        if (is.logical(p) && is.null(dim(p)) && all(is.na(p))) return(as.character(p))
        fail(i)
      })
    })
  if (!identical(class, "list")) out[missing] <- NA
  out
}

## Rebuild the internal row.names attribute recorded by .rich_io_row_names().
.rich_io_restore_row_names <- function(rn, n, slot, bad) {
  type <- if (is.list(rn)) rn$type else NULL
  values <- unlist(rn$values)
  out <- switch(if (is.character(type) && length(type) == 1L) type else "",
    automatic = {
      k <- rn$n
      if (!(is.numeric(k) && length(k) == 1L && is.finite(k) &&
            abs(k) <= .Machine$integer.max && k == trunc(k) && abs(k) == n)) NULL
      else c(NA_integer_, as.integer(k))
    },
    integer = if ((n == 0 && !length(values)) ||
                  (is.numeric(values) && all(is.finite(values)) &&
                   all(abs(values) <= .Machine$integer.max) &&
                   all(values == trunc(values)))) as.integer(values) else NULL,
    character = if (n == 0 || is.character(values)) as.character(values) else NULL,
    NULL)
  if (is.null(out) || (!identical(type, "automatic") && length(out) != n)) {
    bad(sprintf("the row names of %s do not match its row count", slot))
  }
  out
}
