# Backward-compatible aliases for rich* plotting functions
# Users can call either rich* or gg* versions
#
# roxygen silently drops @export for aliases pointing to S4 generics,
# so we force NAMESPACE exports via @rawNamespace.
#' @rawNamespace export(ggbar)
#' @rawNamespace export(ggdot)
#' @rawNamespace export(ggnetplot)
#' @rawNamespace export(ggnetwork)
#' @rawNamespace export(gglollipop)
#' @rawNamespace export(ggLollipop)
#' @rawNamespace export(ggvolcano)
#' @rawNamespace export(ggVolcano)
#' @rawNamespace export(ggcircbar)
#' @rawNamespace export(ggCircBar)
#' @rawNamespace export(ggNES)
#' @rawNamespace export(ggnes)
#' @rawNamespace export(ggscatter)
#' @rawNamespace export(ggScatter)
#' @rawNamespace export(ggecdf)
#' @rawNamespace export(ggECDF)
#' @rawNamespace export(ggtermsim)
#' @rawNamespace export(ggTermSim)
#' @rawNamespace export(gggenedot)
#' @rawNamespace export(ggGeneDot)
#' @rawNamespace export(gggeneheat)
#' @rawNamespace export(ggGeneHeat)
#' @rawNamespace export(gggenebar)
#' @rawNamespace export(ggGeneBar)
NULL

#' @rdname richLollipop-methods
#' @export
gglollipop <- richLollipop

#' @rdname richLollipop-methods
#' @export
ggLollipop <- richLollipop

#' @rdname richVolcano-methods
#' @export
ggvolcano <- richVolcano

#' @rdname richVolcano-methods
#' @export
ggVolcano <- richVolcano

#' @rdname richCircle-methods
#' @export
ggcircbar <- richCircle

#' @rdname richCircle-methods
#' @export
ggCircBar <- richCircle

#' @rdname richNES-methods
#' @export
ggNES <- richNES

#' @rdname richNES-methods
#' @export
ggnes <- richNES

#' @rdname richScatter-methods
#' @export
ggscatter <- richScatter

#' @rdname richScatter-methods
#' @export
ggScatter <- richScatter

#' @rdname richECDF-methods
#' @export
ggecdf <- richECDF

#' @rdname richECDF-methods
#' @export
ggECDF <- richECDF

#' @rdname richTermSim-methods
#' @export
ggtermsim <- richTermSim

#' @rdname richTermSim-methods
#' @export
ggTermSim <- richTermSim

#' @rdname richGeneDot-methods
#' @export
gggenedot <- richGeneDot

#' @rdname richGeneDot-methods
#' @export
ggGeneDot <- richGeneDot

#' @rdname richGeneHeat-methods
#' @export
gggeneheat <- richGeneHeat

#' @rdname richGeneHeat-methods
#' @export
ggGeneHeat <- richGeneHeat

#' @rdname richGeneBar-methods
#' @export
gggenebar <- richGeneBar

#' @rdname richGeneBar-methods
#' @export
ggGeneBar <- richGeneBar

# ============================================================
#  Aliases for renamed existing functions
# ============================================================

# S4 generics
#' @rdname richBar-methods
#' @export
ggbar <- richBar

#' @rdname richDot-methods
#' @export
ggdot <- richDot

#' @rdname richNetplot-method
#' @export
ggnetplot <- richNetplot

#' @rdname richNetwork-methods
#' @export
ggnetwork <- richNetwork

# Regular functions
#' @rdname richNetmap
#' @export
ggnetmap <- richNetmap

#' @rdname richGSEAplot
#' @export
ggGSEA <- richGSEAplot

#' @rdname richGSEAcurve
#' @export
plotGSEA <- richGSEAcurve

#' @rdname richClusterDot
#' @export
ggcluster <- richClusterDot

#' @rdname richUpset
#' @export
ggupset <- richUpset

#' @rdname richHeatmap
#' @export
ggheatmap <- richHeatmap

#' @rdname richCompareDot
#' @export
comparedot <- richCompareDot
