#' Central species lookup table used by all species-mapping functions.
#' Each row: common name, Bioconductor OrgDb, KEGG 3-letter code,
#'           msigdbr scientific name, Reactome scientific name.
#' @return data.frame with columns: species, dbname, kegg, msigdb, reactome
.species_table <- function() {
  data.frame(
    species  = c("anopheles","arabidopsis","bovine","celegans","canine",
                 "chicken","chipm","ecoli","ecsakai","fly",
                 "human","malaria","mouse","pig","rat",
                 "rhesus","streptomyces","toxoplasma","xenopus","yeast",
                 "zebrafish"),
    dbname   = c("org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db",
                 "org.Gg.eg.db","org.Pt.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db","org.Dm.eg.db",
                 "org.Hs.eg.db","org.Pf.plasmo.db","org.Mm.eg.db","org.Ss.eg.db","org.Rn.eg.db",
                 "org.Mmu.eg.db","org.Sco.eg.db","org.Tgondii.eg.db","org.Xl.eg.db","org.Sc.sgd.db",
                 "org.Dr.eg.db"),
    kegg     = c("aga","ath","bta","cel","cfa",
                 "gga","ptr","eco","ecs","dme",
                 "hsa","pfa","mmu","ssc","rno",
                 "mcc","sco",NA,"xla","sce",
                 "dre"),
    msigdb   = c(NA, NA, "Bos taurus","Caenorhabditis elegans","Canis lupus familiaris",
                 "Gallus gallus", NA, NA, NA, "Drosophila melanogaster",
                 "Homo sapiens", NA, "Mus musculus","Sus scrofa","Rattus norvegicus",
                 NA, NA, NA, NA, "Saccharomyces cerevisiae",
                 "Danio rerio"),
    reactome = c(NA, "Arabidopsis thaliana", "Bos taurus","Caenorhabditis elegans","Canis familiaris",
                 "Gallus gallus","Pan troglodytes", NA, NA, "Drosophila melanogaster",
                 "Homo sapiens","Plasmodium falciparum","Mus musculus","Sus scrofa","Rattus norvegicus",
                 "Macaca mulatta", NA, NA, "Xenopus tropicalis","Saccharomyces cerevisiae",
                 "Danio rerio"),
    stringsAsFactors = FALSE
  )
}

#' Look up a single column from the species table
#' @param species common species name
#' @param column column to retrieve
#' @return character value or NULL
#' @keywords internal
.species_lookup <- function(species, column) {
  tbl <- .species_table()
  idx <- match(species, tbl$species)
  if (is.na(idx)) return(NULL)
  val <- tbl[[column]][idx]
  if (is.na(val)) return(NULL)
  val
}

#' Get OrgDb package name for a species
#' @param species common species name
#' @return character OrgDb package name or NULL
#' @keywords internal
.getdb <- function(species) {
  .species_lookup(species, "dbname")
}

#' Get OrgDb package name with error handling
#' @param species common species name
#' @return character OrgDb package name (stops on failure)
#' @keywords internal
.getdbname <- function(species = "human") {
  dbname <- .getdb(species)
  if (is.null(dbname)) {
    cat("You must check if your request database is available by using showData(),
        If not you could make your database by using makeOwnGO and makeOwnKO
        and give a user defined database\n")
    stop("database error!")
  }
  dbname
}

#' Get KEGG 3-letter species code
#' @param species common species name
#' @return character KEGG code or NULL
#' @keywords internal
.getspecies <- function(species = "human") {
  .species_lookup(species, "kegg")
}

#' @rdname dot-getspecies
#' @usage NULL
#' @keywords internal
.getspeices <- .getspecies

#' Get msigdbr scientific name
#' @param species common species name
#' @return character scientific name or NULL
#' @keywords internal
.getmsig <- function(species = "human") {
  .species_lookup(species, "msigdb")
}

#' Get Reactome species name
#' @param species common species name
#' @return character Reactome species name or NULL
#' @keywords internal
.getrodbname <- function(species) {
  .species_lookup(species, "reactome")
}

#' Show available species and their OrgDb packages
#' @return data.frame with dbname and species columns
#' @export
#' @author Kai Guo
showData <- function() {
  tbl <- .species_table()
  data.frame(dbname = tbl$dbname, species = tbl$species, stringsAsFactors = FALSE)
}

#' Convert gene IDs between types using Bioconductor annotation packages
#' @param species species common name (see showData())
#' @param keys character vector of gene IDs to convert
#' @param fkeytype source key type (e.g. "SYMBOL", "ENTREZID")
#' @param tkeytype target key type (e.g. "ENTREZID", "SYMBOL")
#' @return named character vector of converted IDs
#' @examples
#' \dontrun{
#'   id <- idconvert(species="human", keys=c("TP53","BRCA1"),
#'                   fkeytype="SYMBOL", tkeytype="ENTREZID")
#' }
#' @export
#' @author Kai Guo
idconvert <- function(species, keys, fkeytype, tkeytype) {
  dbname <- .getdbname(species)
  if (!requireNamespace(dbname, quietly = TRUE)) stop(paste("Package", dbname, "is required"), call. = FALSE)
  dbname <- getExportedValue(dbname, dbname)
  unlist(AnnotationDbi::mapIds(dbname, keys = as.vector(keys),
                column = tkeytype,
                keytype = fkeytype,
                multiVals = "first"))
}
