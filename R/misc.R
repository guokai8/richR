##' @method as.data.frame Annot
##' @export
as.data.frame.Annot<-function(x,...){
  as.data.frame(x@annot)
}
##' @method as.data.frame richResult
##' @export
as.data.frame.richResult <- function(x, ...) {
  as.data.frame(x@result, ...)
}
##' @method as.data.frame GSEAResult
##' @export
as.data.frame.GSEAResult <- function(x, ...) {
  as.data.frame(x@result, ...)
}


#' S4 row- and column-related methods for Annot, richResult, GSEAResult
#'
#' These methods let you get rownames, colnames, names, etc. for your custom S4 classes.
#'
#' @param x An object of class Annot, richResult, or GSEAResult
#' @param ... Further arguments (ignored)
#'
#' @name s4-accessors
#' @rdname s4-accessors
NULL

# --- row.names() methods ---

#' @rdname s4-accessors
#' @aliases row.names,Annot-method
setMethod(
  "row.names",
  signature(x = "Annot"),
  function(x) {
    row.names(x@annot)
  }
)

#' @rdname s4-accessors
#' @aliases row.names,richResult-method
setMethod(
  "row.names",
  signature(x = "richResult"),
  function(x) {
    row.names(x@result)
  }
)

#' @rdname s4-accessors
#' @aliases row.names,GSEAResult-method
setMethod(
  "row.names",
  signature(x = "GSEAResult"),
  function(x) {
    # if GSEAResult extends richResult
    row.names(x@result)
  }
)


#' @rdname s4-accessors
#' @aliases rownames,Annot-method
setMethod(
  "rownames",
  signature(x = "Annot"),
  function(x) {
    rownames(x@annot)
  }
)

#' @rdname s4-accessors
#' @aliases rownames,richResult-method
setMethod(
  "rownames",
  signature(x = "richResult"),
  function(x) {
    rownames(x@result)
  }
)

#' @rdname s4-accessors
#' @aliases rownames,GSEAResult-method
setMethod(
  "rownames",
  signature(x = "GSEAResult"),
  function(x) {
    rownames(x@result)
  }
)

# --- names() methods ---

#' @rdname s4-accessors
#' @aliases names,Annot-method
setMethod(
  "names",
  signature(x = "Annot"),
  function(x) {
    names(x@annot)
  }
)

#' @rdname s4-accessors
#' @aliases names,richResult-method
setMethod(
  "names",
  signature(x = "richResult"),
  function(x) {
    names(x@result)
  }
)

#' @rdname s4-accessors
#' @aliases names,GSEAResult-method
setMethod(
  "names",
  signature(x = "GSEAResult"),
  function(x) {
    names(x@result)
  }
)

# --- colnames() methods ---

#' @rdname s4-accessors
#' @aliases colnames,Annot-method
setMethod(
  "colnames",
  signature(x = "Annot"),
  function(x) {
    colnames(x@annot)
  }
)

#' @rdname s4-accessors
#' @aliases colnames,richResult-method
setMethod(
  "colnames",
  signature(x = "richResult"),
  function(x) {
    colnames(x@result)
  }
)

#' @rdname s4-accessors
#' @aliases colnames,GSEAResult-method
setMethod(
  "colnames",
  signature(x = "GSEAResult"),
  function(x) {
    colnames(x@result)
  }
)

#' S4 methods for head, tail, dim on Annot, richResult, GSEAResult
#'
#' Defines how to call \code{head}, \code{tail}, and \code{dim} for these S4 classes.
#'
#' @param x An object of class \code{Annot}, \code{richResult}, or \code{GSEAResult}
#' @param n Number of rows to display (for \code{head} or \code{tail})
#' @param ... Further arguments passed to the underlying data.frame method
#'
#' @name utilities-methods
#' @rdname utilities-methods
NULL

# ---------------------------
#  HEAD METHODS
# ---------------------------
#' @rdname utilities-methods
#' @aliases head,Annot-method
setMethod(
  "head",
  signature = signature(x = "Annot"),
  function(x, n = 6L, ...) {
    cat(
      "=== species is:", x@species,
      "and Annotation is", x@anntype,
      "keytype is", x@keytype, "===\n"
    )
    utils::head(x@annot, n, ...)
  }
)

#' @rdname utilities-methods
#' @aliases head,richResult-method
setMethod(
  "head",
  signature = signature(x = "richResult"),
  function(x, n = 6L, ...) {
    cat("=== Total significant terms is:", dim(x@result), "===\n")
    utils::head(x@result, n, ...)
  }
)

#' @rdname utilities-methods
#' @aliases head,GSEAResult-method
setMethod(
  "head",
  signature = signature(x = "GSEAResult"),
  function(x, n = 6L, ...) {
    # GSEAResult extends richResult
    # We can either replicate the code, or call the richResult method
    # Example: call the same code
    cat("=== Total significant terms is:", dim(x@result), "===\n")
    utils::head(x@result, n, ...)
  }
)

# ---------------------------
#  TAIL METHODS
# ---------------------------
#' @rdname utilities-methods
#' @aliases tail,Annot-method
setMethod(
  "tail",
  signature = signature(x = "Annot"),
  function(x, n = 6L, ...) {
    cat(
      "=== species is:", x@species,
      "and Annotation is", x@anntype,
      "keytype is", x@keytype, "===\n"
    )
    utils::tail(x@annot, n, ...)
  }
)

#' @rdname utilities-methods
#' @aliases tail,richResult-method
setMethod(
  "tail",
  signature = signature(x = "richResult"),
  function(x, n = 6L, ...) {
    cat("=== Total significant terms is:", dim(x@result), "===\n")
    utils::tail(x@result, n, ...)
  }
)

#' @rdname utilities-methods
#' @aliases tail,GSEAResult-method
setMethod(
  "tail",
  signature = signature(x = "GSEAResult"),
  function(x, n = 6L, ...) {
    cat("=== Total significant terms is:", dim(x@result), "===\n")
    utils::tail(x@result, n, ...)
  }
)

# ---------------------------
#  DIM METHODS
# ---------------------------
#' @rdname utilities-methods
#' @aliases dim,Annot-method
setMethod(
  "dim",
  signature = signature(x = "Annot"),
  function(x) {
    dim(x@annot)
  }
)

#' @rdname utilities-methods
#' @aliases dim,richResult-method
setMethod(
  "dim",
  signature = signature(x = "richResult"),
  function(x) {
    dim(x@result)
  }
)

#' @rdname utilities-methods
#' @aliases dim,GSEAResult-method
setMethod(
  "dim",
  signature = signature(x = "GSEAResult"),
  function(x) {
    dim(x@result)
  }
)

#' S4 summary method for richResult objects
#'
#' Summarizes the content of a \code{richResult} object.
#'
#' @param object A \code{richResult} object
#' @param ... Not used
#'
#' @return Prints summary info to the console
#' @export
#' @aliases summary,richResult-method
setMethod(
  f = "summary",
  signature = "richResult",
  definition = function(object, ...) {
    cat("Total input genes is:", length(object@gene),
        "and significant biological term is:", nrow(object@result), "\n")
  }
)

#' S4 summary method for GSEAResult objects
#'
#' Summarizes the content of a \code{GSEAResult} object.
#'
#' @param object A \code{GSEAResult} object
#' @param ... Not used
#'
#' @return Prints summary info to the console
#' @export
#' @aliases summary,GSEAResult-method
setMethod(
  f = "summary",
  signature = "GSEAResult",
  definition = function(object, ...) {
    # If GSEAResult extends richResult, you could also callNextMethod().
    cat("Total significant biological term is:",
        table(object@result$padj < 0.05)[[2]], "\n")
  }
)

#' S4 subsetting methods for Annot, richResult, and GSEAResult
#'
#' These methods define `[` subsetting on each respective S4 class.
#'
#' @param x An S4 object (Annot, richResult, or GSEAResult)
#' @param i,j indices for subsetting
#' @param ... further arguments (usually ignored)
#' @param drop logical. Should dimensions be dropped?
#'
#' @name subset-methods
#' @rdname subset-methods
NULL

#' @rdname subset-methods
#' @aliases [,Annot-method
setMethod(
  f = "[",
  signature = signature(x = "Annot", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    x@annot[i, j, drop = drop]
  }
)

#' @rdname subset-methods
#' @aliases [,richResult-method
setMethod(
  f = "[",
  signature = signature(x = "richResult", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    x@result[i, j, drop = drop]
  }
)

#' @rdname subset-methods
#' @aliases [,GSEAResult-method
setMethod(
  f = "[",
  signature = signature(x = "GSEAResult", i = "ANY", j = "ANY", drop = "ANY"),
  definition = function(x, i, j, ..., drop = TRUE) {
    # if GSEAResult extends richResult, you might simply do:
    # callNextMethod() or x@result[i, j, drop=drop]
    x@result[i, j, drop = drop]
  }
)

#' S4 `$` accessors for Annot and richResult
#'
#' Provide `obj$name` access that extracts from the relevant slot columns.
#'
#' @param x An S4 object (Annot or richResult)
#' @param name column name to extract
#'
#' @name dollar-methods
#' @rdname dollar-methods
NULL

#' @rdname dollar-methods
#' @aliases $-Annot
setMethod(
  f = "$",
  signature = signature(x = "Annot"),
  definition = function(x, name) {
    # x@annot is presumably a data.frame or matrix, so
    # use standard R subsetting by column name
    x@annot[[name]]
  }
)

#' @rdname dollar-methods
#' @aliases $-richResult
setMethod(
  f = "$",
  signature = signature(x = "richResult"),
  definition = function(x, name) {
    x@result[[name]]
  }
)

#' @rdname dollar-methods
#' @aliases $-GSEAResult
setMethod(
  f = "$",
  signature = signature(x = "GSEAResult"),
  definition = function(x, name) {
    # If GSEAResult inherits from richResult
    x@result[[name]]
  }
)

#' S4 methods for result() and detail() on richResult and GSEAResult
#'
#' These methods return \code{as.data.frame(x@result)} or \code{x@detail},
#' depending on the class, effectively exposing object contents.
#'
#' @param x A \code{richResult} or \code{GSEAResult} object
#' @return A data frame
#'
#' @name result_detail_methods
#' @rdname result_detail_methods
NULL

#' @rdname result_detail_methods
#' @aliases result,richResult-method
#' @export
setMethod(
  f = "result",
  signature = "richResult",
  definition = function(x) {
    as.data.frame(x@result)
  }
)

#' @rdname result_detail_methods
#' @aliases result,GSEAResult-method
#' @export
setMethod(
  f = "result",
  signature = "GSEAResult",
  definition = function(x) {
    as.data.frame(x@result)
  }
)

#' @rdname result_detail_methods
#' @aliases detail,richResult-method
#' @export
setMethod(
  f = "detail",
  signature = "richResult",
  definition = function(x) {
    as.data.frame(x@detail)
  }
)

#' @rdname result_detail_methods
#' @aliases detail,GSEAResult-method
#' @export
setMethod(
  f = "detail",
  signature = "GSEAResult",
  definition = function(x) {
    as.data.frame(x@detail)
  }
)

##' get detail and integrate with the input gene information
##' @importFrom dplyr left_join
#' @param rese richResult or GSEAResult
#' @param resd dataframe with input gene as rownames
#' @param sep character string used to separate the genes when concatenating
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   hsako<-as.data.frame(hsako)
#'   gene=sample(unique(hsako$GeneID),1000)
#'   res<-richKEGG(gene,kodata = hsako)
#'   gened<-data.frame(lfc=rnorm(length(gene)))
#'   rownames(gened)<-gene
#'   head(getdetail(res,gened,sep=","))
#' }
#' @export
#' @author Kai Guo
getdetail<-function(rese,resd,sep){
  if(!is.data.frame(resd)){
    resd=data.frame(gene=resd)
  }
  if(!("gene"%in%colnames(resd))){
    resd$gene=rownames(resd)
  }
  gene<-strsplit(as.vector(rese$GeneID),split=sep)
  names(gene)<-rese$Annot
  gened<-data.frame("TERM"=rep(names(gene),times=unlist(lapply(gene,length))),
                    "Annot"=rep(rese$Term,times=unlist(lapply(gene,length))),
                    "GeneID"=unlist(gene),row.names=NULL,
                    "Pvalue"=rep(rese$Pvalue,times=unlist(lapply(gene,length))),
                    "Padj"=rep(rese$Padj,times=unlist(lapply(gene,length)))
  )
  gened$GeneID<-as.character(gened$GeneID)
  res<-left_join(gened,resd,by=c("GeneID"="gene"))
  return(res)
}

.color_scale <- function(c1="pink", c2="red") { #modified from DOSE
  pal <- colorRampPalette(c(c1, c2))
  colors <- pal(200)
  return(colors)
}
.getIdx <- function(v, MIN, MAX) { #modified from DOSE
  intervals <- seq(MIN, MAX, length.out=200)
  max(which(intervals <= v))
}

##' @importFrom AnnotationDbi keys
.get_go_dat<-function(ont="BP"){
  require(GO.db)
  key<-keys(GO.db)
  suppressMessages(go_dat<-AnnotationDbi::select(GO.db, keys=key, columns=c("TERM","ONTOLOGY"),keytype="GOID"))
  if(ont=="BP") res<-as.data.frame(subset(go_dat,ONTOLOGY=="BP"))
  if(ont=="CC") res<-as.data.frame(subset(go_dat,ONTOLOGY=="CC"))
  if(ont=="MF") res<-as.data.frame(subset(go_dat,ONTOLOGY=="MF"))
  rownames(res)<-res[,1]
  res<-res[, 2, drop = FALSE]
  colnames(res)<-"annotation"
  return(res)
}
##' @importFrom KEGGREST keggList
.get_kg_dat<-function(builtin=TRUE){
  if(isTRUE(builtin)){
    data(kegg)
    return(kegg.db)
  }else{
    pathway<-cbind(keggList('pathway'))
    rownames(pathway)<-sub('.*map','',rownames(pathway))
    colnames(pathway)<-"annotation"
    pathway<-as.data.frame(pathway)
    pathway$annotation<-as.vector(pathway$annotation)
    #data(pathway)
    if(is.na(pathway['04148',])){
      data(kegg)
      pathway<-kegg.db
    }
    return(pathway)
  }
}
##' @importFrom KEGGREST keggList
##'
.get_kgm.data <- function(){
  #module <-  cbind(keggList('module'))
  #rownames(module)<-sub('md:','',rownames(module))
  #colnames(module)<-"annotation"
  #module<-as.data.frame(module)
  #module$annotation<-as.vector(module$annotation)
  data(module)
  return(module)
}

##' build annotaion for kegg
##' @param ontype GO or KEGG
##' @examples
##' annot = getann("GO")
##' @author Kai Guo
getann<-function(ontype="GO"){
  if(ontype=="GO"){
    res<-rbind(.get_go_dat("BP"),.get_go_dat("MF"),.get_go_dat("CC"))
  }
  if(ontype=="KEGG"){
    res<-.get_kg_dat(builtin=F)
  }
  if(ontype=="Module"){
    res <-.get_kgm_dat()
  }
  return(res)
}

#' reverse List
#' @param lhs: list with names
#' @export
#' @author Kai Guo
# replace this with c++
reverseList_bk<-function(lhs){
  lhs_n<-rep(names(lhs),times=lapply(lhs,function(x)length(x)))
  res<-sf(as.data.frame(cbind(lhs_n,unlist(lhs))))
  return(res)
}
#' ovelap
overlap <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

#' Get all children terms of node
#' @param  node  input node of GO
#' @param  ontology ontology term of BP
#' @author Kai Guo
GO_child <- function(node = "GO:0008150", ontology = "BP") {
  #MF = "GO:0003674", node of MF
  #BP = "GO:0008150", node of BP
  #CC = "GO:0005575", node of CC
  if (ontology == "BP") res <- c(node,GOBPOFFSPRING[[node]])
  if (ontology == "CC") res <- c(node,GOCCOFFSPRING[[node]])
  if (ontology == "MF") res <- c(node,GOMFOFFSPRING[[node]])
  return(res[!is.na(res)])
}

#' Convert ID between ENTREZID to SYMBOL or other type ID based on bioconductor annotation package
#' @param species: you can check the support species by using showData()
#' @param fkeytype: the gene type you want to convert
#' @param tkeytype: the gene type you want to get
#' @examples
#' \dontrun{
#'   hsako<-buildAnnot(species="human",keytype="SYMBOL",anntype = "KEGG")
#'   hsako<-as.data.frame(hsako)
#'   gene=sample(unique(hsako$GeneID),1000)
#'   id<-idconvert(species="human",fkeytype="SYMBOL",tkeytype="ENTREZID")
#' }
#' @export
#' @author Kai Guo
#' @export
#' @author Kai Guo
idconvert<-function(species,keys,fkeytype,tkeytype){
  dbname<-.getdbname(species);
  suppressMessages(require(dbname,character.only = T))
  dbname<-eval(parse(text=dbname))
  unlist(mapIds(dbname,keys=as.vector(keys),
         column=tkeytype,
         keytype=fkeytype,
         multiVals="first"))
}
.getdbname<-function(species="human"){
  dbname=.getdb(species=species);
  if(is.null(dbname)){
    cat("You must check if your request database is avaliable by using showData,
        If not you could make your database by using makeOwnGO and makeOwnKO
        and give a user defined database\n")
    stop("databse error!")
  }
  return(dbname)
}
.getdb<-function(species=species){
  species=tryCatch(match.arg(species,c("anopheles","arabidopsis","bovine","celegans","canine","fly","zebrafish",
                                       "ecoli","ecsakai","chicken","human","mouse","rhesus","malaria","chipm","rat",
                                       "toxoplasma","streptomyces","pig","yeast","xenopus","warm")),
                   error=function(cond){return("unsupported")})
  if (species == "anopheles") {
    dbname <- "org.Ag.eg.db"
  } else if (species == "arabidopsis") {
    dbname <- "org.At.tair.db"
  } else if (species == "bovine") {
    dbname <- "org.Bt.eg.db"
  } else if (species == "canine") {
    dbname <- "org.Cf.eg.db"
  } else if (species == "worm" || species == "celegans") {
    dbname <- "org.Ce.eg.db"
  } else if (species == "chicken") {
    dbname <- "org.Gg.eg.db"
  } else if (species == "ecolik12") {
    dbname <- "org.EcK12.eg.db"
  } else if (species == "ecsakai") {
    dbname <- "org.EcSakai.eg.db"
  } else if (species == "fly") {
    dbname <- "org.Dm.eg.db"
  } else if (species == "human") {
    dbname <- "org.Hs.eg.db"
  } else if (species == "malaria") {
    dbname <- "org.Pf.plasmo.db"
  } else if (species == "chipm") {
    dbname <- "org.Pt.eg.db"
  }else if (species == "mouse") {
    dbname <- "org.Mm.eg.db"
  } else if (species == "pig") {
    dbname <- "org.Ss.eg.db"
  } else if (species == "rat") {
    dbname <- "org.Rn.eg.db"
  } else if (species == "rhesus") {
    dbname <- "org.Mmu.eg.db"
  } else if (species == "xenopus") {
    dbname <- "org.Xl.eg.db"
  } else if (species == "yeast") {
    dbname <- "org.Sc.sgd.db"
  } else if (species == "streptomyces") {
    dbname <- "org.Sco.eg.db"
  } else if (species == "zebrafish") {
    dbname <- "org.Dr.eg.db"
  } else if (species == "toxoplasma"){
    dbname<- "org.Tgondii.eg.db"
  } else {
    dbname <- NULL
  }
  return(dbname)
}
.getspeices<-function(species="human"){
  species=tryCatch(match.arg(species,c("anopheles","arabidopsis","bovine","celegans","canine","fly","zebrafish",
                                       "ecoli","ecsakai","chicken","human","mouse","rhesus","malaria","chipm","rat",
                                       "toxoplasma","sco","pig","yeast","xenopus","warm")),
                   error=function(cond){return("unsupported")})
  if (species == "anopheles") {
    species <- "aga"
  } else if (species == "arabidopsis") {
    species <- "ath"
  } else if (species == "bovine") {
    species <- "bta"
  } else if (species == "canine") {
    species <- "cfa"
  } else if (species == "chicken") {
    species <- "gga"
  } else if (species == "chipm") {
    species <- "ptr"
  } else if (species == "ecolik12") {
    species <- "eco"
  } else if (species == "ecsakai") {
    species <- "ecs"
  } else if (species == "fly") {
    species <- "dme"
  } else if (species == "human") {
    species <- "hsa"
  } else if (species == "malaria") {
    species <- "pfa"
  } else if (species == "mouse") {
    species <- "mmu"
  } else if (species == "pig") {
    species <- "ssc"
  } else if (species == "rat") {
    species <- "rno"
  } else if (species == "rhesus") {
    species <- "mcc"
  } else if (species == "worm" || species == "celegans") {
    species <- "cel"
  } else if (species == "xenopus") {
    species <- "xla"
  } else if (species == "yeast") {
    species <- "sce"
  } else if (species =="streptomyces"){
    species <- "sco"
  } else if (species == "zebrafish") {
    species <- "dre"
  } else {
    species <- NULL
  }
  return(species)
}
#' show avaliable data based on bioconductor annotation package
#' @export
#' @author Kai Guo
showData<-function(){
  species=c("anopheles","arabidopsis","bovine","celegans","canine","fly","zebrafish",
            "ecoli","ecsakai","chicken","human","mouse","rhesus","malaria","chipm","rat",
            "toxoplasma","sco","pig","yeast","xenopus")
  dbname=c("org.Ag.eg.db","org.At.tair.db","org.Bt.eg.db","org.Ce.eg.db","org.Cf.eg.db","org.Dm.eg.db",
           "org.Dr.eg.db","org.EcK12.eg.db","org.EcSakai.eg.db","org.Gg.eg.db","org.Hs.eg.db","org.Mm.eg.db",
           "org.Mmu.eg.db","org.Pf.plasmo.db","org.Pt.eg.db","org.Rn.eg.db","org.Sc.sgd.db","org.Sco.eg.db",
           "org.Ss.eg.db","org.Tgondii.eg.db","org.Xl.eg.db")
  dbdata<-data.frame(dbname=dbname,species=species)
  dbdata
}
##' vector to data.frame
vec_to_df<-function(x,name){
  dd<-data.frame(names(x),x)
  colnames(dd)<-name
  return(dd)
}

##' msigdb support species
##' @param species with common name
.getmsig<-function(species="human"){
  out<-NULL
  if(species=="human"){
    out<-"Homo sapiens"
  }else if(species=="mouse"){
    out<-"Mus musculus"
  }else if(species=="rat"){
    out<-"Rattus norvegicus"
  }else if(species=="celegans"){
    out<-"Caenorhabditis elegans"
  }else if(species=="fly"){
    out<-"rosophila melanogaster"
  }else if(species=="yeast"){
    out<-"Saccharomyces cerevisiae"
  }else if(species=="bovine"){
    out<-"Bos taurus"
  }else if(species=="canine"){
    out<-"Canis lupus familiaris"
  }else if(species=="pig"){
    out<-"Sus scrofa"
  }else if(species=="chicken"){
    out<-"Gallus gallus"
  }else if(species=="zebrafish"){
    out<-"Danio rerio"
  }else{
    out<-NULL
  }
}
##' Print MSIGDB infomation
##' @export
msigdbinfo <- function() {
  cat("#--------------------------------------------------------------#\n")
  cat("# Molecular Signatures Database                        v6.2.1  #\n")
  cat("#--------------------------------------------------------------#\n")
  cat("# Category | Subcategory # Details ----------------------------#\n")
  cat("# C1               # Positional (326)                          #\n")
  cat("# C2 | CGP         # Chemical and Genetic Perturbations (3433) #\n")
  cat("# C2 | CP          # Canonical Pathways (252)                  #\n")
  cat("# C2 | BIOCARTA # Canonical BIOCARTA (217)                     #\n")
  cat("# C2 | KEGG     # Canonical KEGG (186)                         #\n")
  cat("# C2 | CPREACTOME  # Canonical REACTOME (674)                  #\n")
  cat("# C3 | MIR         # Motif miRNA Targets (221)                 #\n")
  cat("# C3 | TFT         # Motif Transcription Factor Targets (615)  #\n")
  cat("# C4 | CGN         # Cancer Gene Neighborhoods (427)           #\n")
  cat("# C4 | CM          # Cancer Modules (431)                      #\n")
  cat("# C5 | BP          # GO Biological Process (4436)              #\n")
  cat("# C5 | CC          # GO Cellular Component (580)               #\n")
  cat("# C5 | MF          # GO Molecular Function (901)               #\n")
  cat("# C6               # Oncogenic Signatures (189)                #\n")
  cat("# C7               # Immunologic Signatures (4872)             #\n")
  cat("# H                # Hallmark (50)                             #\n")
  cat("#--------------------------------------------------------------#\n")
  cat("# Source: http://software.broadinstitute.org/gsea/msigdb       #\n")
  cat("#--------------------------------------------------------------#\n")
  listspe<-c("human","mouse","rat","celegans","fly","yeast","bovine","canine",
             "pig","chicken","zebrafish")
  cat("# Support species:                                             #\n")
  cat(sort(listspe),"\n")
}
.getrodbname<-function(species){
 # "Schizosaccharomyces pombe","Taeniopygia guttata",,"Mycobacterium tuberculosis"
  spe<-c("Homo sapiens","Dictyostelium discoideum","Plasmodium falciparum",
         "Saccharomyces cerevisiae","Caenorhabditis elegans",
         "Sus scrofa","Bos taurus","Canis familiaris","Mus musculus","Rattus norvegicus",
         "Xenopus tropicalis","Danio rerio",
         "Drosophila melanogaster","Arabidopsis thaliana","Oryza sativa","Gallus gallus")
  names(spe)<-c("human","dictyostelium","malaria","yeast","celegans","pig",
                    "bovine","dog","mouse","rat","xenopus","zebrafish","fly","arabidopsis","rice","chicken")
  return(spe[species])
}
##'
setAs(from = "data.frame", to = "Annot", def = function(from){
  keytype <- character()
  species <- character()
  anntype <- character()
  GeneID <- as.vector(from[,1])
  Term<-as.vector(from[,2])
  Annot=from$Annot
  annot <- data.frame(GeneID, Term, Annot)
  new("Annot",
      species = species,
      anntype = anntype,
      keytype = keytype,
      annot = annot
  )
})

##'
##'
setAs(from = "data.frame", to = "richResult", def = function(from){
  keytype <- character()
  organism <- character()
  ontology <- character()
  pvalueCutoff <- numeric()
  pAdjustMethod <-character()
  padjCutoff <- numeric()
  Annot <- from$Annot
  Term <- from$Term
  Annotated <- from$Annotated
  Significant <- from$Significant
  Pvalue <- from$Pvalue
  Padj <- from$Padj
  GeneID <- as.vector(from$GeneID)
  gene<-unique(unlist(strsplit(GeneID,",")))
  genenumber <- length(gene)
  resultFis <- data.frame(Annot, Term, Annotated, Significant, Pvalue, Padj, GeneID)
  rownames(resultFis) <- Annotated
  new("richResult",
      result=resultFis,
      detail=detail,
      pvalueCutoff   = pvalue,
      pAdjustMethod  = padj.method,
      padjCutoff   = padj,
      genenumber    = length(input),
      organism       = organism,
      ontology       = ontology,
      gene           = input,
      keytype        = keytype
  )
})

#' rbind generic function for richResult object
#'@importFrom S4Vectors bindROWS
#'@export
#'@author Kai Guo
rbind.richResult<-function(...){
    objects <- list(...)
    objects <- lapply(objects,as.data.frame)
    bindROWS(objects[[1L]],objects=objects[-1L])
}

#' rbind generic function for GSEAResult object
#'@importFrom S4Vectors bindROWS
#'@export
#'@author Kai Guo
rbind.GSEAResult<-function(...){
  objects <- list(...)
  objects <- lapply(objects,as.data.frame)
  bindROWS(objects[[1L]],objects=objects[-1L])
}

#' Insert newlines after every n words (defaults to 4)
#'
#' Splits a string into words, groups them in chunks of \code{n} words,
#' then rejoins with a newline. If a string has fewer than \code{n} words,
#' it remains unaltered.
#'
#' @param x Character vector to process.
#' @param n Integer. Number of words per line (default 4).
#' @return Character vector with inserted newlines.
#'
#' @examples
#' .paste.char("This is a long string we want to split after four words")
#' .paste.char(c("One phrase", "Another phrase with many words"), n=3)
.paste.char <- function(x, n = 4) {
  sapply(x, function(elem) {
    words <- strsplit(elem, "\\s+")[[1]]
    if (length(words) <= n) {
      return(elem)
    }
    groups <- split(words, ceiling(seq_along(words) / n))
    grouped_strings <- sapply(groups, paste, collapse = " ")
    paste(grouped_strings, collapse = "\n")
  }, USE.NAMES = FALSE)
}

#' Remove newlines
#'
#' Replaces any newline characters with a space.
#'
#' @param x Character vector.
#' @return Character vector without newlines.
#'
#' @examples
#' .clean.char("First line\nSecond line\nThird line")
.clean.char <- function(x) {
  gsub("\\n", " ", x)
}

#' remove the newlines
.clean.char<-function(x){
  return(gsub('\\\n',' ',x))
}

##'
setAs(from = "richResult", to = "data.frame", def = function(from){
  result <- as.data.frame(from@result)
  result
})
##'
setAs(from = "GSEAResult", to = "data.frame", def = function(from){
  result <- as.data.frame(from@result)
  result
})

##'
setAs(from = "Annot", to = "data.frame", def = function(from){
  result <- as.data.frame(from@annot)
  result
})

##' kappa function
.kappa<-function(x,y,geneall){
  x<-unlist(strsplit(x,","))
  y<-unlist(strsplit(y,","))
  if(length(intersect(x,y))==0){
    kab=0
  }else{
    tmp<-matrix(0,2,2)
    tmp[1,1]<-length(intersect(x,y))
    tmp[2,1]<-length(setdiff(x,y))
    tmp[1,2]<-length(setdiff(y,x))
    tmp[2,2]<-length(setdiff(geneall,union(x,y)))
    oab<-(tmp[1,1]+tmp[2,2])/sum(tmp)
    aab<-((tmp[1,1]+tmp[2,1])*(tmp[1,1]+tmp[1,2])+(tmp[1,2]+tmp[2,2])*(tmp[2,1]+tmp[2,2]))/(sum(tmp)*sum(tmp))
    if(aab==1){
      kab=0
    }else{
      kab<-(oab-aab)/(1-aab)
    }
  }
  return(kab)
}
##' calculate enrichment score
.calculate_Enrichment_Score<-function(x,df){
  pvalue <- df[x,"Pvalue"]
  esp = ifelse(pvalue==0,16,-log10(pvalue))
  es = sum(esp);
}

##' merge term
.merge_term<-function(x,overlap){
  ml <- x
  res<-list();
  for(i in names(ml)){
    lhs <- setdiff(names(ml),i)
    for(j in lhs){
      ov<-intersect(ml[[i]],ml[[j]])
      un<-union(ml[[i]],ml[[j]])
      ovl<-length(ov)/length(un)
      if(ovl > overlap){
        res[[i]]<-c(i,un)
        ml <- ml[setdiff(names(ml),j)]
      }else{
        res[[i]]<-c(i,ml[[i]])
      }
    }
  }
  return(res)
}

###
.calGSEA<-function (obj,sigpathway, fc, gseaParam = 1, ticksSize = 0.2){
  obj<-split(obj[,1],obj[,3])
  pathway <- obj[[sigpathway]]
  stats <- fc
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj/max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                                 returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  pathway<-data.frame(x = pathway,Group=sigpathway)
  toPlot$Group<-sigpathway
  return(list(toPlot=toPlot,pathway=pathway,tops=max(tops),bottoms=min(bottoms)))
}

