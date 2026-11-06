#' Stability of ORA conclusions over a prespecified grid of settings
#'
#' \code{richSensitivity()} re-runs \code{\link{enrich}()} once per cell of the
#' full cross product of user-named backgrounds and gene-set upper bounds
#' (\code{maxGSSize}), every other argument identical across cells, and
#' summarises how often each term is reported among the cells that tested it.
#' Every cell is executed independently. In this upstream engine, size bounds
#' are reporting filters applied after adjustment, so changing only a bound
#' can change reported rows without changing adjusted p-values.
#'
#' @section Cells:
#' Cells are ordered by background level, then bound, and named
#' \code{"<background>/<maxGSSize>"}. \code{cells} holds one row per cell:
#' \code{cell}, \code{background}, \code{maxGSSize}, \code{execution}
#' (\code{"success"}, \code{"empty"} or \code{"failed"}), \code{error} (the
#' engine's message for a failed cell), \code{warnings} (a list of the
#' messages the cell raised; each distinct message is signalled once by the
#' grid), \code{background_size} (unique non-missing, non-blank identifiers in
#' the level as passed, \code{NA} for \code{NULL}), \code{effective_N},
#' \code{effective_n} and \code{tested_count} from the cell's metadata
#' (\code{NA} for a failed cell). A failed cell does not stop the grid; a grid
#' whose every cell fails is an error. \code{results} holds each cell's
#' \code{richResult} (\code{NULL} when failed); \code{states} is
#' \code{\link{bindRich}(results)} with \code{background} and \code{maxGSSize}
#' columns after \code{group}.
#'
#' @section Summary:
#' One row per identifier of \code{states}, that is per identifier at least one
#' cell tested (an identifier no cell tests does not appear): \code{Annot};
#' \code{Term} (first
#' non-missing label in cell order); \code{cells}; \code{ran} (cells whose
#' execution is \code{success} or \code{empty}); \code{tested} (state
#' \code{reported} or \code{tested_not_reported}); \code{reported};
#' \code{not_tested}; \code{undetermined} (state
#' \code{untested_or_zero_overlap} or \code{unavailable}); \code{support} =
#' \code{reported / tested} (\code{NA} when \code{tested} is 0);
#' \code{padj_min}, \code{padj_max}, \code{fold_min}, \code{fold_max} over the
#' tested cells' recorded \code{Padj} and \code{FoldEnrichment}; and
#' \code{rank_best}, \code{rank_worst}, the term's \code{Pvalue} rank among each
#' tested cell's tested terms (ties \code{"min"}). Each state is the one
#' \code{bindRich()} assigns from the cell's own record, so an absent
#' identifier is \code{not_tested} in a cell whose \code{zero_overlap_count} is
#' 0 and undetermined where it is above 0. Values no cell records are \code{NA},
#' never 0.
#'
#' The summary is descriptive stability evidence across prespecified settings.
#' It is not a p-value, a confidence interval, an independent replication or a
#' combined test, and \code{richSensitivity()} computes no p-value, adjusted
#' p-value or combined statistic.
#'
#' @param x a gene vector (or a data frame of genes) as for \code{enrich()}; a
#'   numeric ranking is an error (the grid is ORA only), and so are unnamed
#'   numeric identifiers: pass them as character.
#' @param object the annotation passed to \code{enrich()}.
#' @param backgrounds a non-empty named list of background levels: \code{NULL}
#'   for the annotation-derived background or a character vector used as
#'   \code{universe}. Names must be unique, non-blank and free of \code{/}.
#' @param maxGSSize the gene-set upper bounds: unique positive whole numbers.
#' @param ... further arguments passed unchanged to every \code{enrich()} call
#'   (not \code{universe} or \code{maxGSSize}). \code{universe_policy} applies to
#'   every cell; \code{"measured"} needs every level to be a character vector.
#' @return a list of class \code{richSensitivity} with elements \code{cells},
#'   \code{results}, \code{states}, \code{summary} and \code{design}
#'   (\code{backgrounds} as passed, \code{background_sizes}, \code{maxGSSize},
#'   \code{engine}, \code{arguments} (the \code{...} list), \code{criterion},
#'   \code{definitions}, \code{richR_version}, \code{interpretation}). Re-running
#'   \code{richSensitivity(x, object, design$backgrounds, design$maxGSSize,
#'   ...)} with \code{design$arguments} reproduces the grid.
#' @seealso \code{\link{enrich}()}, \code{\link{bindRich}()},
#'   \code{\link{testedTerms}()}
#' @examples
#' set.seed(1)
#' genes <- sprintf("G%03d", 1:300)
#' annot <- data.frame(GeneID = c(sample(genes, 120), sample(genes, 30), sample(genes, 30)),
#'                     Term = rep(c("big", "s1", "s2"), c(120, 30, 30)))
#' query <- c(annot$GeneID[annot$Term == "s1"][1:8], sample(genes, 20))
#' g <- richSensitivity(query, annot, backgrounds = list(all = NULL, part = genes[1:250]),
#'                      maxGSSize = c(100, 200))
#' g$cells
#' g$summary
#' @export
richSensitivity <- function(x, object, backgrounds, maxGSSize = c(500, 1000), ...) {
  args <- list(...)
  if (is.numeric(x)) {
    if (is.null(names(x))) {
      stop("richSensitivity(): `x` holds numeric identifiers; pass gene identifiers as ",
           "character (for example as.character(x)).", call. = FALSE)
    }
    stop("richSensitivity(): `x` is a numeric ranking; the sensitivity grid is ORA only ",
         "and takes a gene vector.", call. = FALSE)
  }
  .grid_check_backgrounds(backgrounds)
  .grid_check_bounds(maxGSSize)
  fixed <- intersect(names(args), c("universe", "maxGSSize"))
  if (length(fixed)) {
    stop(sprintf("richSensitivity(): `%s` is set by the grid and cannot be passed in `...`.",
                 fixed[[1L]]), call. = FALSE)
  }
  if (identical(args$universe_policy, "measured")) {
    null_levels <- names(backgrounds)[vapply(backgrounds, is.null, logical(1))]
    if (length(null_levels)) {
      stop(sprintf(paste("richSensitivity(): universe_policy = \"measured\" needs an",
                         "explicit universe, but `backgrounds` level(s) %s are NULL",
                         "(the annotation-derived background)."),
                   paste0("\"", null_levels, "\"", collapse = ", ")), call. = FALSE)
    }
  }

  levels <- names(backgrounds)
  grid <- data.frame(background = rep(levels, each = length(maxGSSize)),
                     maxGSSize = rep(maxGSSize, times = length(levels)),
                     stringsAsFactors = FALSE)
  cell <- paste0(grid$background, "/", format(grid$maxGSSize, scientific = FALSE, trim = TRUE))
  n <- nrow(grid)
  results <- stats::setNames(vector("list", n), cell)
  errors <- rep(NA_character_, n)
  warns <- rep(list(character(0)), n)
  for (i in seq_len(n)) {
    call_args <- c(list(x, object, universe = backgrounds[[grid$background[i]]],
                        maxGSSize = grid$maxGSSize[i]), args)
    res <- withCallingHandlers(
      tryCatch(do.call(enrich, call_args),
               error = function(e) {
                 errors[[i]] <<- conditionMessage(e)
                 NULL
               }),
      warning = function(w) {
        warns[[i]] <<- c(warns[[i]], conditionMessage(w))
        invokeRestart("muffleWarning")
      })
    if (!is.null(res)) results[[i]] <- res
  }
  failed <- !is.na(errors)
  if (all(failed)) {
    stop("richSensitivity(): every cell failed: ",
         paste0(cell, ": ", errors, collapse = "; "), call. = FALSE)
  }
  for (msg in unique(unlist(warns))) warning(msg, call. = FALSE)

  testing <- lapply(results, function(r) if (!is.null(r)) enrichmentInfo(r)$testing)
  pick <- function(field) vapply(testing, function(t) {
    if (is.null(t)) NA_integer_ else as.integer(t[[field]])
  }, integer(1), USE.NAMES = FALSE)
  sizes <- .grid_background_sizes(backgrounds)
  cells <- data.frame(
    cell = cell, background = grid$background, maxGSSize = grid$maxGSSize,
    execution = ifelse(failed, "failed",
                       ifelse(vapply(results, function(r) !is.null(r) && nrow(r@result) > 0L,
                                     logical(1)), "success", "empty")),
    error = errors, stringsAsFactors = FALSE)
  cells$warnings <- warns
  cells$background_size <- unname(sizes[grid$background])
  cells$effective_N <- pick("effective_N")
  cells$effective_n <- pick("effective_n")
  cells$tested_count <- pick("tested_count")
  rownames(cells) <- NULL

  states <- bindRich(results)
  at <- match(as.character(states$group), cell)
  states$background <- grid$background[at]
  states$maxGSSize <- grid$maxGSSize[at]
  states <- states[c("group", "background", "maxGSSize",
                     setdiff(names(states), c("group", "background", "maxGSSize")))]

  summary <- .grid_summary(states, results, cells)
  first_ok <- enrichmentInfo(results[[which(!failed)[1L]]])$options
  design <- list(
    backgrounds = backgrounds,
    background_sizes = sizes,
    maxGSSize = maxGSSize,
    engine = "enrich",
    arguments = args,
    criterion = list(state = "reported", pvalue = first_ok$pvalue, padj = first_ok$padj,
                     padj.method = first_ok$padj.method,
                     family_policy = enrichmentInfo(results[[which(!failed)[1L]]])$testing$family_policy,
                     universe_policy = first_ok$universe_policy %||% "annotated"),
    definitions = list(
      ran = "cells whose execution is success or empty",
      tested = "cells in which the term's state is reported or tested_not_reported",
      support = "reported / tested; NA when tested is 0"),
    richR_version = .richr_version(),
    interpretation = paste(
      "The summary is descriptive stability evidence across prespecified settings;",
      "it is not a p-value, a confidence interval, an independent replication or a",
      "combined test."))
  structure(list(cells = cells, results = results, states = states,
                 summary = summary, design = design),
            class = "richSensitivity")
}

.grid_check_backgrounds <- function(backgrounds) {
  bad <- function(what) {
    stop("richSensitivity(): `backgrounds` ", what, call. = FALSE)
  }
  if (!is.list(backgrounds) || is.data.frame(backgrounds) || length(backgrounds) == 0L) {
    bad("must be a non-empty named list of background levels.")
  }
  nm <- names(backgrounds)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(trimws(nm)))) {
    bad("must name every level (names must be non-blank).")
  }
  if (anyDuplicated(nm)) {
    bad(sprintf("has the duplicated level name \"%s\".", nm[anyDuplicated(nm)]))
  }
  slash <- grepl("/", nm, fixed = TRUE)
  if (any(slash)) {
    bad(sprintf("level name \"%s\" contains \"/\" (cell ids are \"<background>/<maxGSSize>\").",
                nm[slash][[1L]]))
  }
  ok <- vapply(backgrounds, function(v) is.null(v) || (is.character(v) && is.null(dim(v))),
               logical(1))
  if (!all(ok)) {
    bad(sprintf("level \"%s\" is a %s; a level is NULL or a character vector.",
                nm[!ok][[1L]], class(backgrounds[[which(!ok)[1L]]])[[1L]]))
  }
}

.grid_check_bounds <- function(maxGSSize) {
  ok <- is.numeric(maxGSSize) && is.null(dim(maxGSSize)) && length(maxGSSize) > 0L &&
    !anyNA(maxGSSize) && all(is.finite(maxGSSize)) && all(maxGSSize > 0) &&
    all(maxGSSize == floor(maxGSSize)) && !anyDuplicated(maxGSSize)
  if (!ok) {
    stop(sprintf(paste("richSensitivity(): `maxGSSize` must be a non-empty vector of unique",
                       "positive whole numbers, not %s."),
                 paste(deparse(maxGSSize), collapse = " ")), call. = FALSE)
  }
}

.grid_background_sizes <- function(backgrounds) {
  vapply(backgrounds, function(v) {
    if (is.null(v)) return(NA_integer_)
    length(unique(v[!is.na(v) & nzchar(trimws(v))]))
  }, integer(1))
}

.grid_summary <- function(states, results, cells) {
  ids <- unique(as.character(states$Annot))
  ran_cells <- cells$cell[cells$execution != "failed"]
  ranks <- lapply(results, function(r) {
    if (is.null(r)) return(NULL)
    t <- testedTerms(r)
    stats::setNames(as.integer(rank(t$Pvalue, ties.method = "min")), t$Annot)
  })
  by_id <- split(seq_len(nrow(states)), factor(as.character(states$Annot), levels = ids))
  rows <- lapply(ids, function(id) {
    k <- by_id[[id]]
    st <- states$state[k]
    group <- as.character(states$group[k])
    tested <- st %in% c("reported", "tested_not_reported")
    term <- states$Term[k]
    term <- term[!is.na(term)]
    rng <- function(v, fun) {
      v <- v[tested]
      v <- v[!is.na(v)]
      if (length(v)) fun(v) else NA
    }
    rk <- unlist(lapply(group[tested], function(g) ranks[[g]][id]), use.names = FALSE)
    rk <- rk[!is.na(rk)]
    n_tested <- sum(tested)
    n_reported <- sum(st %in% "reported")
    list(Annot = id,
         Term = if (length(term)) as.character(term[[1L]]) else NA_character_,
         cells = length(k),
         ran = sum(group %in% ran_cells),
         tested = n_tested,
         reported = n_reported,
         not_tested = sum(st %in% "not_tested"),
         undetermined = sum(st %in% c("untested_or_zero_overlap", "unavailable")),
         support = if (n_tested == 0L) NA_real_ else n_reported / n_tested,
         padj_min = as.double(rng(states$Padj[k], min)),
         padj_max = as.double(rng(states$Padj[k], max)),
         fold_min = as.double(rng(states$FoldEnrichment[k], min)),
         fold_max = as.double(rng(states$FoldEnrichment[k], max)),
         rank_best = if (length(rk)) min(rk) else NA_integer_,
         rank_worst = if (length(rk)) max(rk) else NA_integer_)
  })
  column <- function(name, type) {
    vapply(rows, function(r) r[[name]], type, USE.NAMES = FALSE)
  }
  data.frame(
    Annot = column("Annot", character(1)), Term = column("Term", character(1)),
    cells = column("cells", integer(1)), ran = column("ran", integer(1)),
    tested = column("tested", integer(1)), reported = column("reported", integer(1)),
    not_tested = column("not_tested", integer(1)),
    undetermined = column("undetermined", integer(1)),
    support = column("support", double(1)),
    padj_min = column("padj_min", double(1)), padj_max = column("padj_max", double(1)),
    fold_min = column("fold_min", double(1)), fold_max = column("fold_max", double(1)),
    rank_best = column("rank_best", integer(1)), rank_worst = column("rank_worst", integer(1)),
    stringsAsFactors = FALSE)
}
