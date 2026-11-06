#' Write a study report folder from enrichment results
#'
#' Writes a folder a collaborator can open offline: a self-contained HTML
#' report, the saved results, every table behind the report as a
#' tab-separated file, an analysis manifest and a two-line rerun recipe.
#' \code{richReport()} only assembles what the results already record; it
#' performs no new statistical test.
#'
#' @section Inputs:
#' \code{x} is a \code{richResult} or \code{GSEAResult} (reported as
#' \code{Group 1}), a list of them under the group rules of
#' \code{\link{bindRich}} (\code{NULL} entries are failed groups; ORA and GSEA
#' results are not mixed), a character vector of folders written by
#' \code{\link{writeRich}} (group names from the vector's names, else the
#' folders' base names, which must be unique), or one report folder written
#' by \code{richReport()} (see Rerun). A data frame, such as a
#' \code{bindRich()} table, is refused: it carries no analysis record.
#'
#' \code{themes} is \code{NULL} or the \code{\link{richThemes}} object of
#' exactly these groups; only its labels (\code{themes$themes$label}) may have
#' been edited, so it must equal \code{richThemes(x, cutoff)} in every table
#' apart from the labels. \code{stats} (the statistics of
#' \code{\link{themeGenes}}) needs \code{themes}. \code{sensitivity} is
#' \code{NULL} or a \code{\link{richSensitivity}} object whose design holds
#' plain data only. Every check, including the checks of \code{writeRich()}
#' for each group and the column checks of every table, runs before
#' \code{path} is created or changed.
#'
#' @section Folder:
#' \describe{
#'   \item{\code{report.html}}{the report;}
#'   \item{\code{manifest.json}}{format \code{"richR-report"}, version 1;}
#'   \item{\code{rerun.R}}{two lines that call
#'     \code{richR::richReport(".", <new folder>)};}
#'   \item{\code{results/group_<i>/}}{the \code{writeRich()} folder of group
#'     \code{i} (input order, from 1); none for a \code{NULL} group;}
#'   \item{\code{tables/comparison.tsv}}{\code{bindRich(x)};}
#'   \item{\code{tables/audit_group_<i>.tsv}}{\code{mappingAudit()} of group
#'     \code{i}, where it returns a table;}
#'   \item{\code{tables/themes.tsv}, \code{theme_members.tsv},
#'     \code{theme_genes.tsv}, \code{theme_similarity.tsv},
#'     \code{theme_ranking.tsv}}{the five tables of \code{themes};}
#'   \item{\code{tables/theme_gene_view.tsv}}{\code{themeGenes(themes, stats =
#'     stats)};}
#'   \item{\code{tables/theme_stats_group_<i>.tsv}}{the \code{stats} entry of
#'     group \code{i};}
#'   \item{\code{tables/sensitivity_cells.tsv}, \code{sensitivity_summary.tsv},
#'     \code{sensitivity_states.tsv}}{the tables of \code{sensitivity}.}
#' }
#' Tables use the encoding of \code{writeRich()}; \code{manifest.json} records,
#' for each table, its file, column names and classes, factor levels, row
#' names and row count, so the tables read back exactly. It also records
#' \code{format}, \code{format_version}, \code{generated} (the UTC time of the
#' call, \code{YYYY-MM-DDTHH:MM:SSZ}), \code{richR_version},
#' \code{R_version}, \code{title}, \code{figures}, \code{figure_top},
#' \code{groups} (name, folder and class of each group, \code{null} for a
#' failed group), \code{theme_cutoff}, \code{stats_groups},
#' \code{sensitivity_design} (the level names and sizes of the backgrounds,
#' \code{maxGSSize}, \code{engine}, \code{arguments}, \code{criterion},
#' \code{definitions}, \code{richR_version} and \code{interpretation}; not the
#' background identifiers), \code{figures_drawn} (per figure the group, the
#' plot function, the drawn count \code{k}, the reported count \code{n} and
#' the rule), \code{tables} and \code{checksums} (the SHA-256 of every file
#' under \code{tables/} and \code{results/} and of \code{rerun.R}). The
#' generation time is the only value that depends on when the call ran.
#'
#' @section Report:
#' \code{report.html} is one file without scripts, external resources or
#' network access: styles are inline and figures are embedded PNG images.
#' Every text is HTML-escaped; integers and whole numbers are shown in full,
#' other numbers with 4 significant digits, missing values as \code{NA}.
#' Its sections are, in order:
#' \enumerate{
#'   \item a header with the title, the generation time, the richR and R
#'     versions and the statement that the report displays stored results
#'     and performs no new statistical test;
#'   \item per group: the analysis record of \code{\link{enrichmentInfo}} as
#'     key / value rows (nested names joined by \code{$}; vectors of up to 10
#'     values in full, longer vectors, lists without names and tables
#'     summarised by their size), every reported row with every column (each
#'     row has an anchor), a figure, and the mapping audit of
#'     \code{\link{mappingAudit}} as row counts per role and decision, linked
#'     to its full table;
#'   \item a comparison of the identifiers reported in at least one group,
#'     when there are two or more groups, giving each group's
#'     \code{bindRich()} state and adjusted p-value;
#'   \item the themes, when \code{themes} is given: settings, each theme with
#'     its label, representative, member links and gene view;
#'   \item the sensitivity grid, when \code{sensitivity} is given: its design,
#'     its interpretation, and its cells and summary as supplied;
#'   \item how to reproduce the report, with the list of files.
#' }
#' Figures use the internal ORA dot or GSEA NES renderer through a report-only
#' path. Every drawable reported row, including P=1, is eligible before
#' \code{figure.top} is applied, under a local \code{LC_COLLATE = "C"}.
#' Captions and manifest records distinguish drawing omissions from eligible
#' rows truncated by the top limit, and record the actual selected identifiers. A group without a reported row that has a p-value gets no
#' figure, and the report says why. A group missing required plotting
#' columns (\code{Significant} or \code{Annotated} for ORA, \code{NES}
#' for GSEA) also gets no figure; the caption names those columns. Its
#' result table is kept unchanged, with no counts or scores filled in.
#'
#' @section Interpretation:
#' A term significant in one group and not in another is not evidence of a difference between the groups.
#'
#' Theme grouping summarizes the listed terms; it is not a statistical test, and no theme-level p-value exists.
#'
#' Evidence genes are ORA overlap genes or GSEA leading-edge genes; they are not shown to be causal drivers.
#'
#' The mapping audit and the sensitivity grid are shown as recorded or
#' supplied; \code{richReport()} recomputes neither.
#'
#' @section Rerun:
#' rerun.R regenerates this report from the saved results. It does not re-execute the enrichment.
#'
#' Run \code{Rscript rerun.R <new folder>} inside the report folder. Given a
#' report folder, \code{richReport()} accepts only \code{path} and
#' \code{overwrite}; it checks the SHA-256 of every file the manifest lists
#' and stops naming the first file that is missing or changed, reads the
#' groups back with \code{readRich()} (failed groups at their positions),
#' rebuilds the themes with \code{richThemes()} at the recorded cutoff,
#' re-applies the saved labels and stops naming the first theme table that
#' differs from the saved one, and writes the report with the recorded title
#' and figure settings. On the same installation every file under
#' \code{tables/} and \code{results/} and \code{rerun.R} is byte-identical to
#' the original, and \code{report.html} and \code{manifest.json} differ only by
#' the generation time. The gene inputs and the universe of the enrichment
#' are not stored, so a rerun cannot re-execute it. Edits to
#' \code{manifest.json} itself are not detected.
#'
#' @section Overwrite:
#' An existing non-empty \code{path} is left unchanged and the call stops,
#' unless \code{overwrite = TRUE} and \code{path} holds a report folder; then
#' exactly \code{report.html}, \code{manifest.json}, \code{rerun.R},
#' \code{results/} and \code{tables/} are replaced and every other file in
#' \code{path} is kept. Replacement is not transactional: if removal fails,
#' some old files may already be gone. A retry needs the original results
#' or an intact report folder; a partly removed report cannot serve as its
#' own rerun source.
#'
#' @param x results, \code{writeRich()} folders or a report folder; see
#'   Inputs.
#' @param path the report folder to write; created when it does not exist.
#' @param themes \code{NULL} or a \code{richThemes} object of \code{x}.
#' @param stats \code{NULL} or the \code{stats} of \code{themeGenes()}
#'   (needs \code{themes}).
#' @param sensitivity \code{NULL} or a \code{richSensitivity} object.
#' @param title the report title, one string.
#' @param figures \code{TRUE} to draw figures, \code{FALSE} for none.
#' @param figure.top the most rows a figure draws: one positive whole number,
#'   or \code{Inf}.
#' @param overwrite replace a report folder at \code{path}.
#' @return \code{path}, invisibly.
#' @seealso \code{\link{writeRich}}, \code{\link{bindRich}},
#'   \code{\link{mappingAudit}}, \code{\link{richThemes}},
#'   \code{\link{themeGenes}}, \code{\link{richSensitivity}}
#' @examples
#' annot <- data.frame(GeneID = paste0("G", 1:40),
#'                     Term = rep(c("T1", "T2", "T3", "T4"), each = 10),
#'                     Annot = rep(c("set one", "set two", "set three", "set four"),
#'                                 each = 10))
#' a <- enrich(paste0("G", c(1:8, 11:13)), annot, pvalue = 1, padj = NULL,
#'             minSize = 1, minGSSize = 3)
#' b <- enrich(paste0("G", c(2:4, 21:28)), annot, pvalue = 1, padj = NULL,
#'             minSize = 1, minGSSize = 3)
#' out <- file.path(tempdir(), "richReport-example")
#' richReport(list(a = a, b = b), out, figures = FALSE)
#' list.files(out, recursive = TRUE)
#' unlink(out, recursive = TRUE)
#' @export
richReport <- function(x, path, themes = NULL, stats = NULL, sensitivity = NULL,
                       title = "richR enrichment report", figures = TRUE,
                       figure.top = 30, overwrite = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("richReport(): `path` must be a single folder path.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("richReport(): `overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  if (.report_is_folder(x)) {
    given <- c(themes = !missing(themes), stats = !missing(stats),
               sensitivity = !missing(sensitivity), title = !missing(title),
               figures = !missing(figures), figure.top = !missing(figure.top))
    if (any(given)) {
      stop(sprintf(paste("richReport(): `%s` cannot be given with a report folder;",
                         "a rerun uses the values recorded in its manifest.json",
                         "(only `path` and `overwrite` may be given)."),
                   names(given)[given][1L]), call. = FALSE)
    }
    plan <- .report_rerun_plan(x)
  } else {
    .report_check_options(title, figures, figure.top)
    x <- .report_groups(x)
    groups <- names(x)
    if (!is.null(themes)) .report_check_themes(themes, x, groups)
    if (!is.null(stats) && is.null(themes)) {
      stop(paste("richReport(): `stats` needs `themes`; the statistics are shown",
                 "in the theme gene view of themeGenes()."), call. = FALSE)
    }
    view <- if (is.null(sensitivity)) NULL else .report_sensitivity_view(sensitivity)
    plan <- .report_plan(x, themes, stats, view, title, figures, figure.top)
  }
  .report_write(plan, path.expand(path), overwrite)
  invisible(path)
}

.report_format <- "richR-report"
.report_format_version <- 1L

## The entries richReport() writes; the only ones overwrite = TRUE replaces.
.report_known <- c("results", "tables", "report.html", "rerun.R", "manifest.json")

## Stems of the theme tables, by richThemes element.
.report_theme_tables <- c(themes = "themes", members = "theme_members",
                          genes = "theme_genes", similarity = "theme_similarity",
                          ranking = "theme_ranking")

.report_is_folder <- function(path) {
  if (!(is.character(path) && length(path) == 1L && !is.na(path) && dir.exists(path))) {
    return(FALSE)
  }
  f <- file.path(path, "manifest.json")
  if (!file.exists(f) || dir.exists(f)) return(FALSE)
  m <- tryCatch(.rich_io_read_json(f), error = function(e) NULL)
  is.list(m) && identical(m$format, .report_format)
}

.report_check_options <- function(title, figures, figure.top) {
  if (!is.character(title) || length(title) != 1L || is.na(title)) {
    stop("richReport(): `title` must be one non-missing string.", call. = FALSE)
  }
  if (!is.logical(figures) || length(figures) != 1L || is.na(figures)) {
    stop("richReport(): `figures` must be TRUE or FALSE.", call. = FALSE)
  }
  ok <- is.numeric(figure.top) && length(figure.top) == 1L && !is.na(figure.top) &&
    figure.top > 0 && (is.infinite(figure.top) || figure.top == round(figure.top))
  if (!ok) {
    stop("richReport(): `figure.top` must be one positive whole number, or Inf.",
         call. = FALSE)
  }
  invisible(TRUE)
}

## The input as a named list of groups (bindRich() group rules).
.report_groups <- function(x) {
  if (is.data.frame(x)) {
    stop(paste("richReport(): `x` is a data frame (such as a bindRich() table);",
               "a data frame carries no analysis record. Give the result objects",
               "or their writeRich() folders."), call. = FALSE)
  }
  if (is.character(x)) {
    if (!length(x) || anyNA(x)) {
      stop("richReport(): `x` must name one or more folders written by writeRich().",
           call. = FALSE)
    }
    nm <- if (is.null(names(x))) basename(x) else names(x)
    dup <- unique(nm[duplicated(nm)])
    if (length(dup)) {
      stop(sprintf(paste("richReport(): the group name '%s' occurs more than once in",
                         "`x`; name the folders, e.g. c(a = <folder>, b = <folder>)."),
                   dup[1L]), call. = FALSE)
    }
    x <- stats::setNames(lapply(unname(x), readRich), nm)
  }
  if (methods::is(x, "richResult") || methods::is(x, "GSEAResult")) x <- list(x)
  if (!is.list(x)) {
    stop(paste("richReport(): `x` must be a richResult or GSEAResult, a list of them,",
               "folders written by writeRich(), or a report folder."), call. = FALSE)
  }
  bindRich(x)  # stops with bindRich()'s own error for a list it refuses
  if (is.null(names(x))) names(x) <- paste("Group", seq_along(x))
  x
}

.report_check_themes <- function(themes, x, groups) {
  if (!inherits(themes, "richThemes")) {
    stop("richReport(): `themes` must be NULL or a richThemes object, as returned by richThemes().",
         call. = FALSE)
  }
  if (!identical(themes$settings$groups, groups)) {
    stop(sprintf(paste("richReport(): `themes` was built from the groups %s, not from",
                       "the groups of `x` (%s)."),
                 .bind_quote(as.character(themes$settings$groups)), .bind_quote(groups)),
         call. = FALSE)
  }
  cutoff <- themes$settings$cutoff
  rebuilt <- tryCatch(richThemes(x, cutoff = cutoff), error = function(e) {
    stop("richReport(): `themes` cannot be rebuilt from `x` at its cutoff: ",
         conditionMessage(e), call. = FALSE)
  })
  if (identical(rebuilt$themes$theme_id, themes$themes$theme_id)) {
    rebuilt$themes$label <- themes$themes$label
  }
  for (e in c("settings", names(.report_theme_tables))) {
    if (!identical(rebuilt[[e]], themes[[e]])) {
      stop(sprintf(paste("richReport(): `themes` differs from richThemes(x, cutoff = %s)",
                         "in its %s table; only the labels in themes$themes$label may",
                         "be edited."), format(cutoff), e), call. = FALSE)
    }
  }
  invisible(TRUE)
}

## What the report keeps of a richSensitivity object: three tables and the
## plain fields of the design (level names, not the background identifiers).
.report_sensitivity_view <- function(s) {
  if (!inherits(s, "richSensitivity")) {
    stop(paste("richReport(): `sensitivity` must be NULL or a richSensitivity object,",
               "as returned by richSensitivity()."), call. = FALSE)
  }
  d <- s$design
  design <- list(background_names = names(d$backgrounds),
                 background_sizes = d$background_sizes, maxGSSize = d$maxGSSize,
                 engine = d$engine, arguments = d$arguments, criterion = d$criterion,
                 definitions = d$definitions, richR_version = d$richR_version,
                 interpretation = d$interpretation)
  .report_check_plain(design, "design")
  list(cells = s$cells, summary = s$summary, states = s$states, design = design)
}

## Stop, naming the element, on the first value that is not plain data.
.report_check_plain <- function(v, at) {
  extra <- setdiff(names(attributes(v)), "names")
  plain <- is.null(v) ||
    ((is.list(v) || typeof(v) %in% c("logical", "integer", "double", "character")) &&
       !is.object(v) && !length(extra))
  if (!plain) {
    stop(sprintf(paste("richReport(): `sensitivity` %s is not plain data (NULL,",
                       "logical, integer, double, character or lists of these); it",
                       "is a %s."), at, paste(class(v), collapse = "/")), call. = FALSE)
  }
  if (is.list(v)) {
    nm <- names(v)
    for (i in seq_along(v)) {
      label <- if (!is.null(nm) && !is.na(nm[i]) && nzchar(nm[i])) paste0(at, "$", nm[i])
               else sprintf("%s[[%d]]", at, i)
      .report_check_plain(v[[i]], label)
    }
  }
  invisible(TRUE)
}

## One table: its storability checked, its lines and its manifest record.
.report_table <- function(stem, tab, what) {
  if (!identical(class(tab), "data.frame")) {
    stop(sprintf("richReport(): %s is a %s, not a plain data.frame.", what,
                 paste(class(tab), collapse = "/")), call. = FALSE)
  }
  classes <- vapply(seq_along(tab), function(j) {
    k <- .rich_io_column_class(tab[[j]])
    if (is.na(k)) {
      stop(sprintf(paste("richReport(): column '%s' of %s has class %s; supported:",
                         "character, integer, double, logical, factor without an NA",
                         "level and lists of character vectors."),
                   names(tab)[j], what, paste(class(tab[[j]]), collapse = "/")),
           call. = FALSE)
    }
    k
  }, character(1))
  file <- if (ncol(tab)) paste0("tables/", stem, ".tsv") else NULL
  list(entry = list(table = stem, file = file, nrow = nrow(tab), ncol = ncol(tab),
                    columns = .rich_io_column_entries(tab, classes),
                    row_names = .rich_io_row_names(tab)),
       lines = if (ncol(tab)) .rich_io_table_lines(tab, classes))
}

## Decode a table of a report folder with the record its manifest holds.
.report_decode <- function(entry, path) {
  stem <- entry$table
  bad <- function(what) {
    stop(sprintf("richReport(): manifest.json of '%s' is inconsistent: %s.", path, what),
         call. = FALSE)
  }
  if (!(is.character(stem) && length(stem) == 1L)) bad("a table record has no name")
  e <- entry
  e$slot <- stem
  e$file <- NULL
  if (!is.null(entry$file)) {
    if (!identical(entry$file, paste0("tables/", stem, ".tsv"))) {
      bad(sprintf("the file of table %s is not tables/%s.tsv", stem, stem))
    }
    e$file <- paste0(stem, ".tsv")
  }
  .rich_io_read_table(e, file.path(path, "tables"), bad)
}

## The reader the report uses, by table name (used by the tests).
.report_read_table <- function(path, stem) {
  m <- .rich_io_read_json(file.path(path, "manifest.json"))
  entry <- Find(function(e) identical(e$table, stem), m$tables)
  if (is.null(entry)) {
    stop(sprintf("richReport(): manifest.json of '%s' records no table %s.", path, stem),
         call. = FALSE)
  }
  .report_decode(entry, path)
}

## mappingAudit() of one group, with its warnings captured, never re-signalled.
.report_audit <- function(obj) {
  notes <- character(0)
  value <- withCallingHandlers(mappingAudit(obj), warning = function(w) {
    notes <<- c(notes, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  if (is.data.frame(value)) return(list(table = value, notes = notes))
  prefix <- "mappingAudit(): mapping audit unavailable: "
  msg <- notes[startsWith(notes, prefix)][1L]
  reason <- if (is.na(msg)) "no reason given" else sub("\\.$", "", substring(msg, nchar(prefix) + 1L))
  list(table = NULL, reason = reason)
}

## Everything the report folder holds, built and checked in memory.
.report_plan <- function(x, themes, stats, view, title, figures, figure.top) {
  groups <- names(x)
  gsea <- any(vapply(x, methods::is, logical(1), "GSEAResult"))
  id <- if (gsea) "pathway" else "Annot"
  files <- list()
  for (i in seq_along(x)) {
    if (is.null(x[[i]])) next
    enc <- tryCatch(.rich_io_encode(x[[i]]), error = function(e) {
      stop(sprintf("richReport(): group '%s' of `x` cannot be saved: %s", groups[i],
                   conditionMessage(e)), call. = FALSE)
    })
    names(enc) <- paste0("results/group_", i, "/", names(enc))
    files <- c(files, enc)
  }
  bound <- bindRich(x)
  tabs <- list(.report_table("comparison", bound, "the comparison table (bindRich(x))"))
  audits <- lapply(x, function(o) if (is.null(o)) NULL else .report_audit(o))
  for (i in seq_along(x)) {
    a <- audits[[i]]
    if (!is.null(a$table)) {
      tabs[[length(tabs) + 1L]] <- .report_table(
        paste0("audit_group_", i), a$table,
        sprintf("the mapping audit of group '%s'", groups[i]))
    }
  }
  gene_view <- NULL
  if (!is.null(themes)) {
    for (e in names(.report_theme_tables)) {
      tabs[[length(tabs) + 1L]] <- .report_table(
        .report_theme_tables[[e]], themes[[e]], sprintf("`themes` element %s", e))
    }
    gene_view <- if (is.null(stats)) themeGenes(themes) else {
      tryCatch(themeGenes(themes, stats = stats), error = function(e) {
        stop("richReport(): `stats` is refused by themeGenes(): ", conditionMessage(e),
             call. = FALSE)
      })
    }
    tabs[[length(tabs) + 1L]] <- .report_table("theme_gene_view", gene_view,
                                                "the theme gene view (themeGenes())")
    for (g in names(stats)) {
      tabs[[length(tabs) + 1L]] <- .report_table(
        paste0("theme_stats_group_", match(g, groups)), stats[[g]],
        sprintf("`stats` entry '%s'", g))
    }
  }
  if (!is.null(view)) {
    for (e in c("cells", "summary", "states")) {
      tabs[[length(tabs) + 1L]] <- .report_table(
        paste0("sensitivity_", e), view[[e]], sprintf("`sensitivity` element %s", e))
    }
  }
  for (t in tabs) if (!is.null(t$lines)) files[[t$entry$file]] <- t$lines
  files[["rerun.R"]] <- c("args <- commandArgs(trailingOnly = TRUE)",
                          "richR::richReport(\".\", args[[1L]])")
  files <- files[order(names(files), method = "radix")]

  figs <- lapply(seq_along(x), function(i) {
    .report_figure(x[[i]], groups[i], gsea, figures, figure.top)
  })
  generated <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")
  listing <- sort(c(names(files), "report.html", "manifest.json"), method = "radix")
  html <- .report_html(x, groups, gsea, id, bound, audits, figs, themes, gene_view, view,
                       title, generated, listing)
  manifest <- list(
    format = .report_format,
    format_version = .report_format_version,
    generated = generated,
    richR_version = .richr_version(),
    R_version = as.character(getRversion()),
    title = title,
    figures = figures,
    figure_top = if (is.infinite(figure.top)) "Inf" else figure.top,
    groups = lapply(seq_along(x), function(i) {
      list(name = groups[i],
           folder = if (is.null(x[[i]])) NULL else paste0("results/group_", i),
           class = if (is.null(x[[i]])) NULL else class(x[[i]])[1L])
    }),
    theme_cutoff = if (is.null(themes)) NULL else themes$settings$cutoff,
    stats_groups = if (is.null(stats)) NULL else I(names(stats)),
    sensitivity_design = if (is.null(view)) NULL else .rich_io_typed_json(view$design),
    figures_drawn = Filter(Negate(is.null), lapply(figs, `[[`, "record")),
    tables = lapply(tabs, `[[`, "entry")
  )
  list(files = files, html = html, manifest = manifest)
}

## Write a planned folder; the manifest goes last.
.report_write <- function(plan, path, overwrite) {
  if (file.exists(path)) {
    if (!dir.exists(path)) {
      stop(sprintf("richReport(): '%s' exists and is not a folder.", path), call. = FALSE)
    }
    if (length(list.files(path, all.files = TRUE, no.. = TRUE))) {
      if (!overwrite) {
        stop(sprintf(paste("richReport(): folder '%s' is not empty; use overwrite = TRUE",
                           "to replace a report folder written by richReport()."), path),
             call. = FALSE)
      }
      if (!.report_is_folder(path)) {
        stop(sprintf(paste("richReport(): folder '%s' was not written by richReport()",
                           "(no manifest.json naming format '%s'); it is left unchanged."),
                     path, .report_format), call. = FALSE)
      }
      # Keep the identifying manifest until the other entries are removed,
      # so a failed removal leaves the identifying manifest in place.
      # A retry still needs intact inputs; the old report may be incomplete.
      for (entry in .report_known) {
        if (file.access(path, 2L) != 0L ||
            unlink(file.path(path, entry), recursive = TRUE, expand = FALSE) != 0L) {
          stop(sprintf("richReport(): could not remove the old report entries in '%s'.", path),
               call. = FALSE)
        }
      }
    }
  }
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  for (f in names(plan$files)) {
    target <- file.path(path, f)
    dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
    .rich_io_write_lines(plan$files[[f]], target)
  }
  .rich_io_write_lines(plan$html, file.path(path, "report.html"))
  m <- plan$manifest
  m$checksums <- lapply(stats::setNames(nm = names(plan$files)), function(f) {
    digest::digest(file = file.path(path, f), algo = "sha256")
  })
  .rich_io_write_lines(as.character(jsonlite::toJSON(
    m, auto_unbox = TRUE, null = "null", digits = I(17), pretty = TRUE,
    json_verbatim = TRUE)), file.path(path, "manifest.json"))
  invisible(path)
}

## A report folder read back as a plan: checksums verified, groups restored,
## themes rebuilt and checked, nothing written.
.report_rerun_plan <- function(src) {
  m <- .rich_io_read_json(file.path(src, "manifest.json"))
  where <- sprintf("manifest.json of '%s'", src)
  if (!identical(m$format_version, .report_format_version)) {
    stop(sprintf("richReport(): %s records format version %s; supported: %d.", where,
                 .rich_io_found(m$format_version), .report_format_version), call. = FALSE)
  }
  for (f in names(m$checksums)) {
    file <- file.path(src, f)
    if (!file.exists(file) || dir.exists(file)) {
      stop(sprintf("richReport(): %s, listed in %s, is missing.", f, where), call. = FALSE)
    }
    if (!identical(digest::digest(file = file, algo = "sha256"), m$checksums[[f]])) {
      stop(sprintf(paste("richReport(): %s does not match its checksum in %s; the",
                         "report folder was changed after it was written."), f, where),
           call. = FALSE)
    }
  }
  groups <- vapply(m$groups, function(g) as.character(g$name)[1L], character(1))
  x <- vector("list", length(groups))
  names(x) <- groups
  for (i in seq_along(groups)) {
    folder <- m$groups[[i]]$folder
    if (!is.null(folder)) {
      expected <- paste0("results/group_", i)
      if (!identical(folder, expected)) {
        stop(sprintf("richReport(): group %d in %s must use folder %s.",
                     i, where, expected), call. = FALSE)
      }
      x[[i]] <- readRich(file.path(src, folder))
    }
  }
  entries <- m$tables
  names(entries) <- vapply(entries, function(e) as.character(e$table)[1L], character(1))
  decode <- function(stem) {
    if (is.null(entries[[stem]])) {
      stop(sprintf("richReport(): %s records no table %s.", where, stem), call. = FALSE)
    }
    .report_decode(entries[[stem]], src)
  }

  themes <- NULL
  stats <- NULL
  if (!is.null(m$theme_cutoff)) {
    themes <- tryCatch(richThemes(x, cutoff = m$theme_cutoff), error = function(e) {
      stop("richReport(): the themes cannot be rebuilt from the saved results: ",
           conditionMessage(e), call. = FALSE)
    })
    saved <- decode("themes")
    themes$themes$label <- saved$label[match(themes$themes$theme_id, saved$theme_id)]
    for (e in names(.report_theme_tables)) {
      stem <- .report_theme_tables[[e]]
      if (!identical(themes[[e]], decode(stem))) {
        stop(sprintf(paste("richReport(): the themes rebuilt from the saved results",
                           "differ from tables/%s.tsv; the saved labels are not attached",
                           "to different themes."), stem), call. = FALSE)
      }
    }
    if (!is.null(m$stats_groups)) {
      sg <- as.character(unlist(m$stats_groups))
      stats <- stats::setNames(lapply(sg, function(g) {
        decode(paste0("theme_stats_group_", match(g, groups)))
      }), sg)
    }
  }
  view <- NULL
  if (!is.null(m$sensitivity_design)) {
    design <- .rich_io_untyped(m$sensitivity_design, "sensitivity_design", where,
                               function() stop(sprintf("richReport(): the sensitivity_design of %s cannot be decoded.", where), call. = FALSE))
    view <- list(cells = decode("sensitivity_cells"), summary = decode("sensitivity_summary"),
                 states = decode("sensitivity_states"), design = design)
  }
  top <- m$figure_top
  if (identical(top, "Inf")) top <- Inf
  .report_check_options(m$title, m$figures, top)
  .report_plan(x, themes, stats, view, m$title, m$figures, top)
}

## ---- figures ---------------------------------------------------------------

## The figure of one group (an embedded PNG and its manifest record), or the
## caption stating why there is none.
.report_figure <- function(obj, group, gsea, figures, top) {
  none <- function(why) list(caption = paste("No figure:", why))
  if (is.null(obj)) return(none("the group failed."))
  if (!figures) return(none("figures = FALSE."))
  if (!isTRUE(unname(capabilities("png")))) return(none("the png device is unavailable."))
  pcol <- if (gsea) "pval" else "Pvalue"
  r <- obj@result
  n <- nrow(r)
  if (!n || !pcol %in% names(r) || all(is.na(r[[pcol]]))) {
    return(none("no reported term has a p-value."))
  }
  required <- if (gsea) "NES" else c("Significant", "Annotated")
  absent <- setdiff(required, names(r))
  if (length(absent)) {
    return(none(paste0("required plotting columns missing: ", paste(absent, collapse = ", "), ".")))
  }
  fun <- if (gsea) "richNES" else "richDot"
  rule <- if (gsea) sprintf("stored rows with finite NES; first %s by NES, descending", .report_num(top))
          else sprintf("stored rows with valid drawing values; first %s in result order", .report_num(top))
  drawn <- .report_draw(obj, gsea, top)
  k <- if (is.null(drawn$k)) n else drawn$k
  list(caption = sprintf(paste0("%s of %s reported terms; %s; ",
                              "%s omitted for missing, blank or non-finite drawing values; ",
                              "%s truncated by top"),
                         .report_num(k), .report_num(n), rule,
                         .report_num(drawn$omitted), .report_num(drawn$truncated)),
       png = drawn$png, alt = sprintf("%s of group %s", fun, group),
       record = list(group = group, `function` = fun, k = k, n = n, rule = rule, ids = drawn$ids,
                     eligible = drawn$eligible, omitted = drawn$omitted,
                     truncated = drawn$truncated,
                     omission_reason = if (drawn$omitted > 0L)
                       "missing, blank or non-finite drawing values" else NULL))
}

## Draw with every reported row eligible (pvalue = 1, padj = NULL) under a
## local LC_COLLATE = "C", so the term axis does not follow the session.
.report_draw <- function(obj, gsea, top) {
  old <- Sys.getlocale("LC_COLLATE")
  on.exit(Sys.setlocale("LC_COLLATE", old), add = TRUE)
  Sys.setlocale("LC_COLLATE", "C")
  k <- NULL
  selected <- NULL
  p <- withCallingHandlers(
    if (gsea) ggNES_internal(obj@result, pvalue = 1, padj = NULL, top = top, .report_all = TRUE)
    else ggdot_internal(obj@result, pvalue = 1, padj = NULL, top = top,
                 usePadj = "Padj" %in% names(obj@result), .report_all = TRUE),
    richR_plot_selection = function(cond) {
      k <<- cond$k
      selected <<- cond
      invokeRestart("muffleMessage")
    })
  # Each layer gets its own defaults; update_geom_defaults() in the caller
  # must neither change a saved report nor be changed by drawing one.
  for (i in seq_along(p$layers)) {
    geom <- p$layers[[i]]$geom
    defaults <- if (inherits(geom, "GeomPoint")) {
      ggplot2::aes(shape = 19, colour = "black", fill = NA, size = 1.5,
                   alpha = NA, stroke = 0.5)
    } else if (inherits(geom, "GeomBar")) {
      ggplot2::aes(colour = NA, fill = "grey35", linewidth = 0.5,
                   linetype = 1, alpha = NA, width = 0.9)
    } else NULL
    if (!is.null(defaults)) {
      p$layers[[i]]$geom <- ggplot2::ggproto(NULL, geom, default_aes = defaults)
    }
  }
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)
  type <- if (isTRUE(unname(capabilities("cairo")))) "cairo" else {
    if (.Platform$OS.type == "windows") "windows" else {
      if (isTRUE(unname(capabilities("aqua")))) "quartz" else "Xlib"
    }
  }
  grDevices::png(file, width = 800, height = 600, res = 96, type = type)
  dev <- grDevices::dev.cur()
  tryCatch(print(p), finally = grDevices::dev.off(dev))
  bytes <- readBin(file, "raw", file.info(file)$size)
  list(k = k, ids = selected$ids, eligible = selected$eligible,
       omitted = selected$omitted, truncated = selected$truncated, png = gsub("[\r\n]", "", jsonlite::base64_enc(bytes)))
}

## ---- HTML --------------------------------------------------------------------

# @MX:NOTE: [AUTO] the one escape helper every text of report.html passes through (REQ-011)
.report_escape <- function(x) {
  x <- enc2utf8(as.character(x))
  x[is.na(x)] <- "NA"
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  gsub("'", "&#39;", x, fixed = TRUE)
}

## A start tag: lower-case name, every attribute value escaped and
## double-quoted.
.report_open <- function(name, attrs = NULL) {
  a <- if (length(attrs)) {
    paste0(" ", names(attrs), "=\"", .report_escape(attrs), "\"", collapse = "")
  } else ""
  paste0("<", name, a, ">")
}

## An element holding escaped text (or `html`, already built from escaped parts).
.report_el <- function(name, text = NULL, attrs = NULL, html = NULL) {
  inner <- if (is.null(html)) .report_escape(text) else html
  paste0(.report_open(name, attrs), inner, "</", name, ">")
}

## Integers and whole numbers below 2^53 in full; other doubles with 4
## significant digits; NA as NA; NaN, Inf, -Inf spelled out.
.report_num <- function(v) {
  if (is.integer(v)) {
    out <- as.character(v)
    out[is.na(v)] <- "NA"
    return(out)
  }
  out <- character(length(v))
  na <- is.na(v) & !is.nan(v)
  whole <- is.finite(v) & v == trunc(v) & abs(v) < 2^53
  out[whole] <- sprintf("%.0f", v[whole])
  out[whole & v == 0] <- "0"
  other <- !whole & !na
  out[other] <- sprintf("%.4g", v[other])
  out[na] <- "NA"
  out
}

## Display text of each element of an atomic vector.
.report_values <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (is.numeric(v)) return(.report_num(v))
  if (is.logical(v)) return(ifelse(is.na(v), "NA", ifelse(v, "TRUE", "FALSE")))
  out <- as.character(v)
  out[is.na(out)] <- "NA"
  out
}

## Display text of each cell of a table column (list cells joined by ", ").
.report_cells <- function(v) {
  if (is.list(v)) {
    return(vapply(v, function(e) {
      if (is.null(e)) "NULL" else paste(.report_values(unlist(e)), collapse = ", ")
    }, character(1)))
  }
  .report_values(v)
}

## Summary of one record value (REQ-006).
.report_summary <- function(v) {
  if (is.null(v)) return("NULL")
  if (is.data.frame(v)) return(sprintf("table of %d rows and %d columns", nrow(v), ncol(v)))
  if (is.list(v)) return(sprintf("list of %d values", length(v)))
  if (is.atomic(v)) {
    cls <- class(v)[1L]
    if (!length(v)) return(sprintf("%s vector of 0 values", cls))
    if (length(v) > 10L) return(sprintf("%s vector of %d values", cls, length(v)))
    return(paste(.report_values(v), collapse = ", "))
  }
  sprintf("%s value", class(v)[1L])
}

## Key / value rows of a nested record; key paths join names with `$`.
.report_record <- function(x, prefix = "") {
  out <- list(key = character(0), value = character(0))
  nm <- names(x)
  for (i in seq_along(x)) {
    key <- paste0(prefix, nm[i])
    v <- x[[i]]
    named <- is.list(v) && !is.data.frame(v) && length(v) && !is.null(names(v)) &&
      !anyNA(names(v)) && all(nzchar(names(v)))
    if (named) {
      sub <- .report_record(v, paste0(key, "$"))
      out$key <- c(out$key, sub$key)
      out$value <- c(out$value, sub$value)
    } else {
      out$key <- c(out$key, key)
      out$value <- c(out$value, .report_summary(v))
    }
  }
  out
}

## An HTML table from a header and columns of already-escaped cell HTML.
.report_html_rows <- function(header, cols, ids = NULL) {
  head <- paste0("<tr>", paste0("<th>", .report_escape(header), "</th>", collapse = ""), "</tr>")
  n <- if (length(cols)) length(cols[[1L]]) else 0L
  body <- character(0)
  if (n) {
    cells <- do.call(paste0, lapply(cols, function(v) paste0("<td>", v, "</td>")))
    open <- if (is.null(ids)) rep("<tr>", n) else
      vapply(ids, function(i) .report_open("tr", c(id = i)), character(1))
    body <- paste0(open, cells, "</tr>")
  }
  c("<table>", head, body, "</table>")
}

.report_html_table <- function(df, ids = NULL) {
  .report_html_rows(names(df), lapply(df, function(v) .report_escape(.report_cells(v))), ids)
}

.report_kv_table <- function(rec) {
  .report_html_rows(c("key", "value"), list(.report_escape(rec$key), .report_escape(rec$value)))
}

## Row anchor: g<i>- and the lowercase hex of the identifier's UTF-8 bytes.
.report_anchor <- function(i, ids) {
  paste0("g", i, "-", vapply(as.character(ids), function(v) {
    paste(as.character(charToRaw(enc2utf8(v))), collapse = "")
  }, character(1), USE.NAMES = FALSE))
}

.report_p <- function(text) .report_el("p", text)

.report_style <- c(
  "body{font-family:sans-serif;margin:2em;max-width:80em}",
  "table{border-collapse:collapse;margin:0.5em 0}",
  "th,td{border:1px solid #bbb;padding:2px 6px;text-align:left;vertical-align:top}",
  "th{background:#f2f2f2}",
  "img{max-width:100%}",
  "code{background:#f2f2f2;padding:0 2px}")

.report_html <- function(x, groups, gsea, id, bound, audits, figs, themes, gene_view,
                         view, title, generated, listing) {
  h <- c("<!DOCTYPE html>", "<html>", "<head>", .report_open("meta", c(charset = "utf-8")),
         .report_el("title", title), "<style>", .report_style, "</style>", "</head>",
         "<body>", .report_el("h1", title),
         .report_p(paste0("Generated: ", generated)),
         .report_p(sprintf("richR %s; R %s", .richr_version(), as.character(getRversion()))),
         .report_p(paste("This report displays stored enrichment results.",
                         "It performs no new statistical test.")))
  pcol <- if (gsea) "padj" else "Padj"
  for (i in seq_along(x)) {
    o <- x[[i]]
    h <- c(h, .report_el("h2", sprintf("Group %d: %s", i, groups[i]),
                         c(id = paste0("group-", i))),
           .report_el("h3", "Analysis record"))
    if (is.null(o)) {
      h <- c(h, .report_p("Execution: failed."), .report_el("h3", "Figure"),
             .report_figure_html(figs[[i]]))
      next
    }
    info <- enrichmentInfo(o)
    h <- c(h, if (isTRUE(info$available)) .report_kv_table(.report_record(info))
              else .report_p(sprintf("Analysis record unavailable: %s.",
                                     paste(info$reason, collapse = " "))))
    res <- o@result
    h <- c(h, .report_el("h3", "Reported terms"),
           .report_p(sprintf("Reported terms: %d rows.", nrow(res))),
           .report_html_table(res, .report_anchor(i, res[[id]])),
           .report_el("h3", "Figure"), .report_figure_html(figs[[i]]),
           .report_el("h3", "Mapping audit"), .report_audit_html(audits[[i]], i))
  }
  if (length(x) >= 2L) h <- c(h, .report_comparison_html(x, groups, id, pcol, bound))
  if (!is.null(themes)) h <- c(h, .report_themes_html(themes, gene_view, groups, id))
  if (!is.null(view)) h <- c(h, .report_sensitivity_html(view))
  links <- vapply(listing, function(f) .report_el("a", f, c(href = f)), character(1))
  h <- c(h, .report_el("h2", "Reproduce", c(id = "reproduce")),
         .report_p("rerun.R regenerates this report from the saved results. It does not re-execute the enrichment."),
         .report_el("p", html = paste0(
           "Run ", .report_el("code", "Rscript rerun.R <new folder>"),
           " inside this folder; it writes the report again into a new, empty folder.")),
         .report_p("Files of this report:"),
         .report_html_rows("file", list(unname(links))),
         "</body>", "</html>")
  h
}

.report_figure_html <- function(f) {
  if (is.null(f$png)) {
    return(c("<figure>", .report_el("figcaption", f$caption), "</figure>"))
  }
  c("<figure>",
    .report_open("img", c(alt = f$alt, src = paste0("data:image/png;base64,", f$png))),
    .report_el("figcaption", f$caption), "</figure>")
}

.report_audit_html <- function(a, i) {
  if (is.null(a$table)) {
    return(.report_p(sprintf("Mapping audit unavailable: %s.", a$reason)))
  }
  tab <- a$table
  role <- as.character(tab$role)
  decision <- as.character(tab$decision)
  pairs <- unique(data.frame(role = role, decision = decision, stringsAsFactors = FALSE))
  pairs <- pairs[order(pairs$role, pairs$decision, method = "radix"), , drop = FALSE]
  key <- function(r, d) paste(ifelse(is.na(r), "\001", r), ifelse(is.na(d), "\001", d),
                              sep = "\002")
  rows <- tabulate(match(key(role, decision), key(pairs$role, pairs$decision)),
                   nbins = nrow(pairs))
  keep <- rows > 0L
  file <- sprintf("tables/audit_group_%d.tsv", i)
  c(.report_p(sprintf("Mapping audit: %d rows.", nrow(tab))),
    vapply(a$notes, function(n) .report_p(paste("Mapping audit note:", n)), character(1),
           USE.NAMES = FALSE),
    .report_html_rows(c("role", "decision", "rows"),
                      list(.report_escape(.report_values(pairs$role[keep])),
                           .report_escape(.report_values(pairs$decision[keep])),
                           .report_num(rows[keep]))),
    .report_el("p", html = paste0("Full audit: ", .report_el("a", file, c(href = file)))),
    .report_p("Shown as recorded; not recomputed by richReport()."))
}

.report_comparison_html <- function(x, groups, id, pcol, bound) {
  ids <- as.character(bound[[id]])
  grp <- as.character(bound$group)
  shown <- sort(unique(ids[bound$state %in% "reported"]), method = "radix")
  all_ids <- unique(ids[!is.na(ids)])
  cols <- lapply(groups, function(g) {
    if (is.null(x[[g]])) return(rep("failed", length(shown)))
    rows <- which(grp == g)
    at <- rows[match(shown, ids[rows])]
    state <- bound$state[at]
    with_p <- state %in% c("reported", "tested_not_reported")
    text <- .report_values(state)
    p <- if (pcol %in% names(bound)) bound[[pcol]][at] else rep(NA_real_, length(at))
    text[with_p] <- sprintf("%s (%s)", state[with_p], .report_num(p[with_p]))
    .report_escape(text)
  })
  c(.report_el("h2", "Comparison", c(id = "comparison")),
    .report_p(sprintf("%d identifiers reported in at least one group, of %d in tables/comparison.tsv.",
                      length(shown), length(all_ids))),
    .report_p("A term significant in one group and not in another is not evidence of a difference between the groups."),
    .report_html_rows(c("identifier", groups), c(list(.report_escape(shown)), cols)))
}

.report_themes_html <- function(themes, gene_view, groups, id) {
  s <- themes$settings
  h <- c(.report_el("h2", "Themes", c(id = "themes")),
         .report_p("Theme grouping summarizes the listed terms; it is not a statistical test, and no theme-level p-value exists."),
         .report_p("Evidence genes are ORA overlap genes or GSEA leading-edge genes; they are not shown to be causal drivers."),
         .report_p("Related or overlapping gene sets are not independent replications."),
         .report_kv_table(.report_record(list(method = s$method, linkage = s$linkage,
                                              cutoff = s$cutoff,
                                              background_size = s$background_size))))
  mem <- themes$members
  for (r in seq_len(nrow(themes$themes))) {
    t <- themes$themes$theme_id[r]
    rows <- mem[mem$theme_id == t, , drop = FALSE]
    ids <- unique(as.character(rows[[id]]))
    links <- vapply(ids, function(j) {
      rep_rows <- rows[rows[[id]] %in% j & rows$state %in% "reported", , drop = FALSE]
      g <- min(match(as.character(rep_rows$group), groups))
      .report_el("a", j, c(href = paste0("#", .report_anchor(g, j))))
    }, character(1), USE.NAMES = FALSE)
    rep_term <- themes$themes$representative_term[r]
    h <- c(h,
           .report_el("h3", sprintf("%s (%s)", .report_values(themes$themes$label[r]), t),
                      c(id = t)),
           .report_p(paste0("Representative: ", themes$themes$representative[r],
                            if (!is.na(rep_term)) paste0(" (", rep_term, ")"))),
           .report_el("p", html = paste0("Members: ", paste(links, collapse = ", "))),
           .report_html_table(gene_view[gene_view$theme_id == t, , drop = FALSE]))
  }
  h
}

.report_sensitivity_html <- function(view) {
  d <- view$design
  c(.report_el("h2", "Sensitivity", c(id = "sensitivity")),
    .report_kv_table(.report_record(list(backgrounds = d$background_names,
                                         background_sizes = d$background_sizes,
                                         maxGSSize = d$maxGSSize,
                                         criterion = d$criterion))),
    .report_p(d$interpretation),
    .report_p("Shown as supplied; not recomputed by richReport()."),
    .report_p(sprintf("Cells: %d rows.", nrow(view$cells))),
    .report_html_table(view$cells),
    .report_p(sprintf("Summary: %d rows.", nrow(view$summary))),
    .report_html_table(view$summary))
}

# Report-only renderer accounting. The private caller opts into stored-row
# eligibility; ordinary public plot methods retain their strict cutoffs.
.report_plot_selection <- function(selected, eligible, submitted, id) {
  ids <- as.character(selected[[id]])
  message(structure(list(message = "report plot selection", call = NULL,
                         k = as.integer(nrow(selected)), ids = ids,
                         eligible = as.integer(eligible),
                         omitted = as.integer(submitted - eligible),
                         truncated = as.integer(eligible - nrow(selected))),
                    class = c("richR_plot_selection", "message", "condition")))
}
