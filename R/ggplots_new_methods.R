# ============================================================
#  S4 generics + methods for new plotting functions
#
#  Primary names use rich* prefix.
#  Old gg* names are kept as aliases for backward compatibility.
# ============================================================

# ----------------------------------------------------------
#  richLollipop  (aliases: gglollipop, ggLollipop)
# ----------------------------------------------------------
#' Lollipop plot of enrichment significance
#'
#' @param object richResult, data.frame, or GSEAResult
#' @param top number of terms (default 20)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param usePadj colour by adjusted p-value?
#' @param low low-end gradient colour
#' @param high high-end gradient colour
#' @param point.size point size
#' @param short shorten term names
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   hsako <- buildAnnot(species="human", keytype="SYMBOL", anntype="KEGG")
#'   gene  <- sample(unique(hsako$GeneID), 1000)
#'   res   <- richKEGG(gene, kodata = hsako)
#'   richLollipop(res)
#' }
#' @rdname richLollipop-methods
#' @export
setGeneric("richLollipop", function(object, top = 20, pvalue = 0.05, padj = NULL,
                                     usePadj = TRUE, low = "lightpink", high = "red",
                                     point.size = 4,
                                     short = FALSE, filename = NULL, width = 10, height = 8, ...)
  standardGeneric("richLollipop"))

#' @rdname richLollipop-methods
#' @export
setMethod("richLollipop", signature(object = "richResult"), function(object, top = 20, pvalue = 0.05, padj = NULL,
                                     usePadj = TRUE, low = "lightpink", high = "red",
                                     point.size = 4,
                                     short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  gglollipop_internal(object@result, top = top, pvalue = pvalue, padj = padj,
                      usePadj = usePadj, low = low, high = high, point.size = point.size,
                      short = short, filename = filename, width = width, height = height)
})

#' @rdname richLollipop-methods
#' @export
setMethod("richLollipop", signature(object = "data.frame"), function(object, top = 20, pvalue = 0.05, padj = NULL,
                                     usePadj = TRUE, low = "lightpink", high = "red",
                                     point.size = 4,
                                     short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  gglollipop_internal(object, top = top, pvalue = pvalue, padj = padj,
                      usePadj = usePadj, low = low, high = high, point.size = point.size,
                      short = short, filename = filename, width = width, height = height)
})

#' @rdname richLollipop-methods
#' @export
setMethod("richLollipop", signature(object = "GSEAResult"), function(object, top = 20, pvalue = 0.05, padj = NULL,
                                     usePadj = TRUE, low = "lightpink", high = "red",
                                     point.size = 4,
                                     short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  df <- .gsea_to_ora_frame(object)
  gglollipop_internal(df, top = top, pvalue = pvalue, padj = padj,
                      usePadj = usePadj, low = low, high = high, point.size = point.size,
                      short = short, filename = filename, width = width, height = height)
})

# ----------------------------------------------------------
#  richVolcano  (aliases: ggvolcano, ggVolcano)
# ----------------------------------------------------------
#' Volcano plot of enrichment results
#'
#' For ORA (richResult): x = log2(Fold Enrichment), y = -log10(p), color = effect.
#' For GSEA (GSEAResult): x = NES, y = -log10(p), color = NES.
#' Point size encodes significance; labels auto-adjusted via ggrepel.
#'
#' @param object richResult, GSEAResult, or data.frame
#' @param top max terms (default 50)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param usePadj use adjusted p-value for y-axis
#' @param low low-end gradient colour (negative effect)
#' @param mid midpoint gradient colour
#' @param high high-end gradient colour (positive effect)
#' @param size.range size range for -log10(p) mapping
#' @param label.size label text size
#' @param label.top number of top terms to label (via ggrepel)
#' @param short shorten term names
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   resgo <- richGO(gene, godata = hsago, ontology = "BP")
#'   richVolcano(resgo)
#'   resgs <- richGSEA(genelist, gseadata)
#'   richVolcano(resgs)
#' }
#' @rdname richVolcano-methods
#' @export
setGeneric("richVolcano", function(object, top = 50, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "#2166ac", mid = "white", high = "#b2182b",
                                    size.range = c(2, 8), label.size = 3, label.top = 5,
                                    short = FALSE, filename = NULL, width = 10, height = 8, ...)
  standardGeneric("richVolcano"))

#' @rdname richVolcano-methods
#' @export
setMethod("richVolcano", signature(object = "richResult"), function(object, top = 50, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "#2166ac", mid = "white", high = "#b2182b",
                                    size.range = c(2, 8), label.size = 3, label.top = 5,
                                    short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggvolcano_internal(object@result, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, mid = mid, high = high,
                     size.range = size.range, label.size = label.size, label.top = label.top,
                     short = short, filename = filename, width = width, height = height)
})

#' @rdname richVolcano-methods
#' @export
setMethod("richVolcano", signature(object = "GSEAResult"), function(object, top = 50, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "#2166ac", mid = "white", high = "#b2182b",
                                    size.range = c(2, 8), label.size = 3, label.top = 5,
                                    short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggvolcano_internal(object@result, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, mid = mid, high = high,
                     size.range = size.range, label.size = label.size, label.top = label.top,
                     short = short, filename = filename, width = width, height = height)
})

#' @rdname richVolcano-methods
#' @export
setMethod("richVolcano", signature(object = "data.frame"), function(object, top = 50, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "#2166ac", mid = "white", high = "#b2182b",
                                    size.range = c(2, 8), label.size = 3, label.top = 5,
                                    short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggvolcano_internal(object, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, mid = mid, high = high,
                     size.range = size.range, label.size = label.size, label.top = label.top,
                     short = short, filename = filename, width = width, height = height)
})

# ----------------------------------------------------------
#  richCircle  (aliases: ggcircbar, ggCircBar)
# ----------------------------------------------------------
#' Circular barplot of enrichment results
#'
#' @param object richResult, data.frame, or GSEAResult
#' @param top number of terms (default 15)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param usePadj colour by adjusted p-value?
#' @param low low-end gradient colour
#' @param high high-end gradient colour
#' @param short shorten term names
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   res <- richKEGG(gene, kodata = hsako)
#'   richCircle(res)
#' }
#' @rdname richCircle-methods
#' @export
setGeneric("richCircle", function(object, top = 15, pvalue = 0.05, padj = NULL,
                                   usePadj = TRUE, low = "#fee0d2", high = "#b2182b",
                                   short = FALSE, filename = NULL, width = 10, height = 8, ...)
  standardGeneric("richCircle"))

#' @rdname richCircle-methods
#' @export
setMethod("richCircle", signature(object = "richResult"), function(object, top = 15, pvalue = 0.05, padj = NULL,
                                   usePadj = TRUE, low = "#fee0d2", high = "#b2182b",
                                   short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggcircbar_internal(object@result, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, high = high, short = short,
                     filename = filename, width = width, height = height)
})

#' @rdname richCircle-methods
#' @export
setMethod("richCircle", signature(object = "data.frame"), function(object, top = 15, pvalue = 0.05, padj = NULL,
                                   usePadj = TRUE, low = "#fee0d2", high = "#b2182b",
                                   short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggcircbar_internal(object, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, high = high, short = short,
                     filename = filename, width = width, height = height)
})

#' @rdname richCircle-methods
#' @export
setMethod("richCircle", signature(object = "GSEAResult"), function(object, top = 15, pvalue = 0.05, padj = NULL,
                                   usePadj = TRUE, low = "#fee0d2", high = "#b2182b",
                                   short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  df <- .gsea_to_ora_frame(object)
  ggcircbar_internal(df, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, high = high, short = short,
                     filename = filename, width = width, height = height)
})

# ----------------------------------------------------------
#  richNES  (aliases: ggNES, ggnes)
# ----------------------------------------------------------
#' Bidirectional NES barplot of GSEA results
#'
#' @param object GSEAResult or data.frame with NES column
#' @param top number of terms (default 20)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param up.color colour for positive NES
#' @param down.color colour for negative NES
#' @param short shorten term names
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   res <- richGSEA(genelist, gseadata)
#'   richNES(res)
#' }
#' @rdname richNES-methods
#' @export
setGeneric("richNES", function(object, top = 20, pvalue = 0.05, padj = NULL,
                                up.color = "#b2182b", down.color = "#2166ac",
                                short = FALSE, filename = NULL, width = 10, height = 8, ...)
  standardGeneric("richNES"))

#' @rdname richNES-methods
#' @export
setMethod("richNES", signature(object = "GSEAResult"), function(object, top = 20, pvalue = 0.05, padj = NULL,
                                up.color = "#b2182b", down.color = "#2166ac",
                                short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggNES_internal(object@result, top = top, pvalue = pvalue, padj = padj,
                 up.color = up.color, down.color = down.color, short = short,
                 filename = filename, width = width, height = height)
})

#' @rdname richNES-methods
#' @export
setMethod("richNES", signature(object = "data.frame"), function(object, top = 20, pvalue = 0.05, padj = NULL,
                                up.color = "#b2182b", down.color = "#2166ac",
                                short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggNES_internal(object, top = top, pvalue = pvalue, padj = padj,
                 up.color = up.color, down.color = down.color, short = short,
                 filename = filename, width = width, height = height)
})

# ----------------------------------------------------------
#  richScatter  (aliases: ggscatter, ggScatter)
# ----------------------------------------------------------
#' Enrichment funnel plot (MA-plot analogy)
#'
#' Plots pathway size vs fold enrichment, coloured by significance.
#' Analogous to an MA plot in differential expression: small pathways
#' (left) can have extreme fold enrichment by chance, while large
#' pathways (right) with high enrichment are the most robust findings.
#'
#' @param object richResult, data.frame, or GSEAResult
#' @param top number of terms (default 50)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param usePadj colour by adjusted p-value?
#' @param low low-end gradient colour
#' @param high high-end gradient colour
#' @param point.size point size
#' @param label.size label text size
#' @param label.top number of top terms to label (via ggrepel)
#' @param short shorten term names
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   res <- richKEGG(gene, kodata = hsako)
#'   richScatter(res)
#' }
#' @rdname richScatter-methods
#' @export
setGeneric("richScatter", function(object, top = 50, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "#fee0d2", high = "#b2182b",
                                    point.size = 3, label.size = 3, label.top = 5,
                                    short = FALSE, filename = NULL, width = 10, height = 8, ...)
  standardGeneric("richScatter"))

#' @rdname richScatter-methods
#' @export
setMethod("richScatter", signature(object = "richResult"), function(object, top = 50, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "#fee0d2", high = "#b2182b",
                                    point.size = 3, label.size = 3, label.top = 5,
                                    short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggscatter_internal(object@result, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, high = high, point.size = point.size,
                     label.size = label.size, label.top = label.top, short = short,
                     filename = filename, width = width, height = height)
})

#' @rdname richScatter-methods
#' @export
setMethod("richScatter", signature(object = "data.frame"), function(object, top = 50, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "#fee0d2", high = "#b2182b",
                                    point.size = 3, label.size = 3, label.top = 5,
                                    short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggscatter_internal(object, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, high = high, point.size = point.size,
                     label.size = label.size, label.top = label.top, short = short,
                     filename = filename, width = width, height = height)
})

#' @rdname richScatter-methods
#' @export
setMethod("richScatter", signature(object = "GSEAResult"), function(object, top = 50, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "#fee0d2", high = "#b2182b",
                                    point.size = 3, label.size = 3, label.top = 5,
                                    short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  df <- .gsea_to_ora_frame(object)
  ggscatter_internal(df, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, high = high, point.size = point.size,
                     label.size = label.size, label.top = label.top, short = short,
                     filename = filename, width = width, height = height)
})

# ----------------------------------------------------------
#  richECDF  (aliases: ggecdf, ggECDF)
# ----------------------------------------------------------
#' ECDF step plot of GSEA ranked gene sets
#'
#' @param object GSEAResult or data.frame with NES column
#' @param top number of terms (default 30)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param line.color line colour
#' @param line.size line size
#' @param point.size point size
#' @param short shorten term names
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   res <- richGSEA(genelist, gseadata)
#'   richECDF(res)
#' }
#' @rdname richECDF-methods
#' @export
setGeneric("richECDF", function(object, top = 30, pvalue = 0.05, padj = NULL,
                                 line.color = "#1b7837", line.size = 1.2, point.size = 3,
                                 short = FALSE, filename = NULL, width = 10, height = 8, ...)
  standardGeneric("richECDF"))

#' @rdname richECDF-methods
#' @export
setMethod("richECDF", signature(object = "GSEAResult"), function(object, top = 30, pvalue = 0.05, padj = NULL,
                                 line.color = "#1b7837", line.size = 1.2, point.size = 3,
                                 short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggecdf_internal(object@result, top = top, pvalue = pvalue, padj = padj,
                  line.color = line.color, line.size = line.size, point.size = point.size,
                  short = short, filename = filename, width = width, height = height)
})

#' @rdname richECDF-methods
#' @export
setMethod("richECDF", signature(object = "data.frame"), function(object, top = 30, pvalue = 0.05, padj = NULL,
                                 line.color = "#1b7837", line.size = 1.2, point.size = 3,
                                 short = FALSE, filename = NULL, width = 10, height = 8, ...) {
  ggecdf_internal(object, top = top, pvalue = pvalue, padj = padj,
                  line.color = line.color, line.size = line.size, point.size = point.size,
                  short = short, filename = filename, width = width, height = height)
})

# ----------------------------------------------------------
#  richTermSim  (aliases: ggtermsim, ggTermSim)
# ----------------------------------------------------------
#' Enrichment map — term similarity network
#'
#' Positions terms by MDS of gene-set overlap (Jaccard) and draws
#' edges between terms that share genes above a similarity cutoff.
#' Node size = gene count, colour = -log10(p). Connected clusters
#' reveal groups of functionally related pathways.
#'
#' @param object richResult, data.frame, or GSEAResult
#' @param top number of terms (default 20)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param usePadj colour by adjusted p-value?
#' @param low low-end gradient colour
#' @param high high-end gradient colour
#' @param label.size label size
#' @param sim.cutoff minimum Jaccard similarity to draw an edge (0-1)
#' @param size.range point size range for gene count mapping
#' @param short shorten term names
#' @param sep gene-ID separator
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   res <- richKEGG(gene, kodata = hsako)
#'   richTermSim(res)
#'   richTermSim(res, sim.cutoff = 0.3)
#' }
#' @rdname richTermSim-methods
#' @export
setGeneric("richTermSim", function(object, top = 20, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "grey80", high = "#b2182b",
                                    label.size = 3, sim.cutoff = 0.2, size.range = c(3, 10),
                                    short = FALSE, sep = ",",
                                    filename = NULL, width = 10, height = 8, ...)
  standardGeneric("richTermSim"))

#' @rdname richTermSim-methods
#' @export
setMethod("richTermSim", signature(object = "richResult"), function(object, top = 20, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "grey80", high = "#b2182b",
                                    label.size = 3, sim.cutoff = 0.2, size.range = c(3, 10),
                                    short = FALSE, sep = ",",
                                    filename = NULL, width = 10, height = 8, ...) {
  ggtermsim_internal(object@result, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, high = high, label.size = label.size,
                     sim.cutoff = sim.cutoff, size.range = size.range,
                     short = short, sep = object@sep, filename = filename,
                     width = width, height = height)
})

#' @rdname richTermSim-methods
#' @export
setMethod("richTermSim", signature(object = "data.frame"), function(object, top = 20, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "grey80", high = "#b2182b",
                                    label.size = 3, sim.cutoff = 0.2, size.range = c(3, 10),
                                    short = FALSE, sep = ",",
                                    filename = NULL, width = 10, height = 8, ...) {
  ggtermsim_internal(object, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, high = high, label.size = label.size,
                     sim.cutoff = sim.cutoff, size.range = size.range,
                     short = short, sep = sep, filename = filename,
                     width = width, height = height)
})

#' @rdname richTermSim-methods
#' @export
setMethod("richTermSim", signature(object = "GSEAResult"), function(object, top = 20, pvalue = 0.05, padj = NULL,
                                    usePadj = TRUE, low = "grey80", high = "#b2182b",
                                    label.size = 3, sim.cutoff = 0.2, size.range = c(3, 10),
                                    short = FALSE, sep = ",",
                                    filename = NULL, width = 10, height = 8, ...) {
  df <- .gsea_to_ora_frame(object)
  ggtermsim_internal(df, top = top, pvalue = pvalue, padj = padj,
                     usePadj = usePadj, low = low, high = high, label.size = label.size,
                     sim.cutoff = sim.cutoff, size.range = size.range,
                     short = short, sep = object@sep, filename = filename,
                     width = width, height = height)
})

# ----------------------------------------------------------
#  richGeneDot  (aliases: gggenedot, ggGeneDot)
# ----------------------------------------------------------
#' Gene-pathway association dotplot
#'
#' @param object richResult, data.frame, or GSEAResult
#' @param fc named numeric vector of fold-changes (names = gene IDs). When NULL, colour by -log10(Pvalue/Padj)
#' @param top number of terms (default 10)
#' @param max.genes max genes per term (default 30)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param usePadj logical; use adjusted p-value for colour gradient? (default TRUE)
#' @param low colour for negative FC (or low significance when fc=NULL)
#' @param mid colour for zero FC
#' @param high colour for positive FC (or high significance when fc=NULL)
#' @param short shorten term names
#' @param sep gene-ID separator
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   res <- richKEGG(gene, kodata = hsako)
#'   richGeneDot(res, fc = my_foldchanges)
#'   richGeneDot(res)  # no fc: colour by significance
#' }
#' @rdname richGeneDot-methods
#' @export
setGeneric("richGeneDot", function(object, fc = NULL, top = 10, max.genes = 30,
                                    pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                    low = "#2166ac", mid = "white", high = "#b2182b",
                                    short = FALSE, sep = ",",
                                    filename = NULL, width = 12, height = 8, ...)
  standardGeneric("richGeneDot"))

#' @rdname richGeneDot-methods
#' @export
setMethod("richGeneDot", signature(object = "richResult"), function(object, fc = NULL, top = 10, max.genes = 30,
                                    pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                    low = "#2166ac", mid = "white", high = "#b2182b",
                                    short = FALSE, sep = ",",
                                    filename = NULL, width = 12, height = 8, ...) {
  gggenedot_internal(object@result, fc = fc, top = top, max.genes = max.genes,
                     pvalue = pvalue, padj = padj, usePadj = usePadj, low = low, mid = mid, high = high,
                     short = short, sep = object@sep, filename = filename,
                     width = width, height = height)
})

#' @rdname richGeneDot-methods
#' @export
setMethod("richGeneDot", signature(object = "data.frame"), function(object, fc = NULL, top = 10, max.genes = 30,
                                    pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                    low = "#2166ac", mid = "white", high = "#b2182b",
                                    short = FALSE, sep = ",",
                                    filename = NULL, width = 12, height = 8, ...) {
  gggenedot_internal(object, fc = fc, top = top, max.genes = max.genes,
                     pvalue = pvalue, padj = padj, usePadj = usePadj, low = low, mid = mid, high = high,
                     short = short, sep = sep, filename = filename,
                     width = width, height = height)
})

#' @rdname richGeneDot-methods
#' @export
setMethod("richGeneDot", signature(object = "GSEAResult"), function(object, fc = NULL, top = 10, max.genes = 30,
                                    pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                    low = "#2166ac", mid = "white", high = "#b2182b",
                                    short = FALSE, sep = ",",
                                    filename = NULL, width = 12, height = 8, ...) {
  df <- .gsea_to_ora_frame(object)
  gggenedot_internal(df, fc = fc, top = top, max.genes = max.genes,
                     pvalue = pvalue, padj = padj, usePadj = usePadj, low = low, mid = mid, high = high,
                     short = short, sep = object@sep, filename = filename,
                     width = width, height = height)
})

# ----------------------------------------------------------
#  richGeneHeat  (aliases: gggeneheat, ggGeneHeat)
# ----------------------------------------------------------
#' Gene-term heatmap
#'
#' @param object richResult, data.frame, or GSEAResult
#' @param fc named numeric vector of fold-changes. When NULL, fill by -log10(Pvalue/Padj)
#' @param top number of terms (default 10)
#' @param max.genes max genes per term (default 30)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param usePadj logical; use adjusted p-value for colour gradient? (default TRUE)
#' @param low colour for negative FC (or white when fc=NULL)
#' @param mid colour for zero FC
#' @param high colour for positive FC (or high significance when fc=NULL)
#' @param na.fill tile fill colour for gene-term combinations with no data (default \code{"grey90"})
#' @param border.color tile border colour (default \code{"grey40"})
#' @param label.cutoff absolute FC above which gene names are shown
#' @param short shorten term names
#' @param sep gene-ID separator
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   res <- richKEGG(gene, kodata = hsako)
#'   richGeneHeat(res, fc = my_foldchanges)
#'   richGeneHeat(res)  # no fc: fill by significance
#' }
#' @rdname richGeneHeat-methods
#' @export
setGeneric("richGeneHeat", function(object, fc = NULL, top = 10, max.genes = 30,
                                     pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                     low = "#2166ac", mid = "white", high = "#b2182b",
                                     na.fill = "grey90", border.color = "grey40",
                                     label.cutoff = 1.5, short = FALSE, sep = ",",
                                     filename = NULL, width = 12, height = 8, ...)
  standardGeneric("richGeneHeat"))

#' @rdname richGeneHeat-methods
#' @export
setMethod("richGeneHeat", signature(object = "richResult"), function(object, fc = NULL, top = 10, max.genes = 30,
                                     pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                     low = "#2166ac", mid = "white", high = "#b2182b",
                                     na.fill = "grey90", border.color = "grey40",
                                     label.cutoff = 1.5, short = FALSE, sep = ",",
                                     filename = NULL, width = 12, height = 8, ...) {
  gggeneheat_internal(object@result, fc = fc, top = top, max.genes = max.genes,
                      pvalue = pvalue, padj = padj, usePadj = usePadj, low = low, mid = mid, high = high,
                      na.fill = na.fill, border.color = border.color,
                      label.cutoff = label.cutoff, short = short, sep = object@sep,
                      filename = filename, width = width, height = height)
})

#' @rdname richGeneHeat-methods
#' @export
setMethod("richGeneHeat", signature(object = "data.frame"), function(object, fc = NULL, top = 10, max.genes = 30,
                                     pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                     low = "#2166ac", mid = "white", high = "#b2182b",
                                     na.fill = "grey90", border.color = "grey40",
                                     label.cutoff = 1.5, short = FALSE, sep = ",",
                                     filename = NULL, width = 12, height = 8, ...) {
  gggeneheat_internal(object, fc = fc, top = top, max.genes = max.genes,
                      pvalue = pvalue, padj = padj, usePadj = usePadj, low = low, mid = mid, high = high,
                      na.fill = na.fill, border.color = border.color,
                      label.cutoff = label.cutoff, short = short, sep = sep,
                      filename = filename, width = width, height = height)
})

#' @rdname richGeneHeat-methods
#' @export
setMethod("richGeneHeat", signature(object = "GSEAResult"), function(object, fc = NULL, top = 10, max.genes = 30,
                                     pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                     low = "#2166ac", mid = "white", high = "#b2182b",
                                     na.fill = "grey90", border.color = "grey40",
                                     label.cutoff = 1.5, short = FALSE, sep = ",",
                                     filename = NULL, width = 12, height = 8, ...) {
  df <- .gsea_to_ora_frame(object)
  gggeneheat_internal(df, fc = fc, top = top, max.genes = max.genes,
                      pvalue = pvalue, padj = padj, usePadj = usePadj, low = low, mid = mid, high = high,
                      na.fill = na.fill, border.color = border.color,
                      label.cutoff = label.cutoff, short = short, sep = object@sep,
                      filename = filename, width = width, height = height)
})

# ----------------------------------------------------------
#  richGeneBar  (aliases: gggenebar, ggGeneBar)
# ----------------------------------------------------------
#' Barplot with gene names embedded inside bars
#'
#' @param object richResult, data.frame, or GSEAResult
#' @param top number of terms (default 10)
#' @param pvalue p-value cutoff
#' @param padj adjusted p-value cutoff
#' @param usePadj logical; use adjusted p-value for colour gradient? (default TRUE)
#' @param low low-end gradient colour (default \code{"lightpink"})
#' @param high high-end gradient colour (default \code{"red"})
#' @param max.genes max gene names inside each bar
#' @param bar.color bar colour (unused; bars coloured by significance)
#' @param bar.alpha bar transparency
#' @param text.size gene label size
#' @param short shorten term names
#' @param sep gene-ID separator
#' @param filename save path
#' @param width figure width
#' @param height figure height
#' @param ... additional arguments
#' @return ggplot2 object
#' @examples
#' \dontrun{
#'   res <- richKEGG(gene, kodata = hsako)
#'   richGeneBar(res)
#' }
#' @rdname richGeneBar-methods
#' @export
setGeneric("richGeneBar", function(object, top = 10, pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                    low = "lightpink", high = "red",
                                    max.genes = 5, bar.color = "#b2182b", bar.alpha = 0.85,
                                    text.size = 3, short = FALSE, sep = ",",
                                    filename = NULL, width = 10, height = 8, ...)
  standardGeneric("richGeneBar"))

#' @rdname richGeneBar-methods
#' @export
setMethod("richGeneBar", signature(object = "richResult"), function(object, top = 10, pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                    low = "lightpink", high = "red",
                                    max.genes = 5, bar.color = "#b2182b", bar.alpha = 0.85,
                                    text.size = 3, short = FALSE, sep = ",",
                                    filename = NULL, width = 10, height = 8, ...) {
  gggenebar_internal(object@result, top = top, pvalue = pvalue, padj = padj, usePadj = usePadj,
                     low = low, high = high,
                     max.genes = max.genes, bar.color = bar.color, bar.alpha = bar.alpha,
                     text.size = text.size, short = short, sep = object@sep,
                     filename = filename, width = width, height = height)
})

#' @rdname richGeneBar-methods
#' @export
setMethod("richGeneBar", signature(object = "data.frame"), function(object, top = 10, pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                    low = "lightpink", high = "red",
                                    max.genes = 5, bar.color = "#b2182b", bar.alpha = 0.85,
                                    text.size = 3, short = FALSE, sep = ",",
                                    filename = NULL, width = 10, height = 8, ...) {
  gggenebar_internal(object, top = top, pvalue = pvalue, padj = padj, usePadj = usePadj,
                     low = low, high = high,
                     max.genes = max.genes, bar.color = bar.color, bar.alpha = bar.alpha,
                     text.size = text.size, short = short, sep = sep,
                     filename = filename, width = width, height = height)
})

#' @rdname richGeneBar-methods
#' @export
setMethod("richGeneBar", signature(object = "GSEAResult"), function(object, top = 10, pvalue = 0.05, padj = NULL, usePadj = TRUE,
                                    low = "lightpink", high = "red",
                                    max.genes = 5, bar.color = "#b2182b", bar.alpha = 0.85,
                                    text.size = 3, short = FALSE, sep = ",",
                                    filename = NULL, width = 10, height = 8, ...) {
  df <- .gsea_to_ora_frame(object)
  gggenebar_internal(df, top = top, pvalue = pvalue, padj = padj, usePadj = usePadj,
                     low = low, high = high,
                     max.genes = max.genes, bar.color = bar.color, bar.alpha = bar.alpha,
                     text.size = text.size, short = short, sep = object@sep,
                     filename = filename, width = width, height = height)
})

# ============================================================
#  Helper: convert GSEAResult to ORA-like data.frame
# ============================================================

#' @noRd
.gsea_to_ora_frame <- function(object) {
  df <- object@result
  sep <- if (length(object@sep) && nchar(object@sep)) object@sep else ","
  if ("pval"    %in% names(df) && !"Pvalue" %in% names(df)) names(df)[names(df) == "pval"]    <- "Pvalue"
  if ("padj"    %in% names(df) && !"Padj"   %in% names(df)) names(df)[names(df) == "padj"]    <- "Padj"
  if ("pathway" %in% names(df) && !"Term"   %in% names(df)) names(df)[names(df) == "pathway"] <- "Term"
  if ("size"    %in% names(df) && !"Significant" %in% names(df)) names(df)[names(df) == "size"] <- "Significant"
  if (!"Annotated" %in% names(df) && "Significant" %in% names(df)) df$Annotated <- df$Significant
  if (!"RichFactor" %in% names(df) && "Significant" %in% names(df)) df$RichFactor <- 1
  if ("leadingEdge" %in% names(df)) {
    df$GeneID <- vapply(df$leadingEdge, function(x) {
      if (is.list(x)) paste(unlist(x), collapse = sep)
      else if (is.character(x)) paste(x, collapse = sep)
      else ""
    }, character(1))
  }
  if (!"Annot" %in% names(df) && "Term" %in% names(df)) df$Annot <- df$Term
  df
}
