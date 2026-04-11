##' build annotation database
##' @name buildAnnot
##' @rdname buildAnnot-methods
##' @title make annotation database
##' @param species species for the annotation
##' @param keytype key type export
##' @param anntype annotation type
##' @param builtin use default database (TRUE or FALSE)
##' @param OP BP,CC,MF default use all
##' @examples
##' \dontrun{
##' annot<-buildAnnot(species="human",keytype="ENTREZID",anntype="GO",builtin=TRUE)
##' }
##' @export
##' @author Kai Guo
buildAnnot<-function(species="human",keytype="SYMBOL",anntype="GO",builtin=TRUE,OP=NULL){
  message("buildAnnot: building ", anntype, " annotation for ", species, " (keytype: ", keytype, ")...")
  annot <- switch(anntype,
    GO       = .makeGOdata(species = species, keytype = keytype, OP = OP),
    KEGG     = .makeKOdata(species = species, keytype = keytype, builtin = builtin),
    Reactome = .makeROdata(species = species, keytype = keytype),
    KEGGM    = .makeKOMdata(species = species, keytype = keytype),
    stop("Unsupported anntype '", anntype,
         "'. Use one of: GO, KEGG, Reactome, KEGGM")
  )
  n_genes <- length(unique(annot[, 1]))
  n_terms <- length(unique(annot[, 2]))
  message("buildAnnot: done. ", n_genes, " genes, ", n_terms, " terms.")
  new("Annot",
      species = species,
      anntype = anntype,
      keytype = keytype,
      annot = annot)
}
#' make GO annotation data function
#' @importFrom AnnotationDbi keys
#' @param species you can check the support species by using showData()
#' @param keytype the gene ID type
#' @param OP BP,CC,MF default use all
#' @importFrom dplyr distinct
#' @author Kai Guo
.makeGOdata<-function(species="human",keytype="ENTREZID",OP=NULL){
  dbname<-.getdbname(species);
  if (!requireNamespace(dbname, quietly = TRUE)) {
    stop("Package '", dbname, "' is required. Install it with:\n",
         "  BiocManager::install('", dbname, "')", call. = FALSE)
  }
  dbname<-getExportedValue(dbname, dbname)
  GO_FILE<-AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns=c("GOALL","ONTOLOGYALL"))
  colnames(GO_FILE)[1]<-"GeneID"
  #GO_FILE<-distinct_(GO_FILE,~GeneID, ~GOALL, ~ONTOLOGYALL)
  GO_FILE <- distinct(GO_FILE[,c('GeneID','GOALL','ONTOLOGYALL')])
  annot <- getann("GO")
  GO_FILE$Annot <- annot[GO_FILE[,2],"annotation"]
  if(!is.null(OP)){
    GO_FILE<-GO_FILE[GO_FILE$ONTOLOGYALL==OP,]
  }
  return(GO_FILE)
}
#' make KEGG annotation data function
#' @importFrom AnnotationDbi keys
#' @importFrom KEGGREST keggLink
#' @param species you can check the support species by using showData()
#' @param keytype the gene ID type
#' @param builtin use KEGG built-in annotation or not (default TRUE)
#' @author Kai Guo
.makeKOdata<-function(species="human",keytype="ENTREZID",builtin=TRUE){
  dbname<-.getdbname(species=species);
  if(builtin==TRUE){
    if (!requireNamespace(dbname, quietly = TRUE)) {
      stop("Package '", dbname, "' is required. Install it with:\n",
           "  BiocManager::install('", dbname, "')", call. = FALSE)
    }
    dbname<-getExportedValue(dbname, dbname)
    KO_FILE=AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns="PATH")
    KO_FILE<-na.omit(KO_FILE)
  }else{
    spe=.getspeices(species)
    tmp<-keggLink("pathway",spe)
    tmp<-substr(tmp,9,13)
    names(tmp)<-sub('.*:','',names(tmp))
    tmp<-vec_to_df(tmp,name=c(keytype,"PATH"))
    if(keytype!="ENTREZID"){
      tmp[,1]<-idconvert(species,keys=tmp[,1],fkeytype = "ENTREZID",tkeytype = keytype)
      tmp<-na.omit(tmp)
    }
    KO_FILE=tmp
  }
  annot<-getann("KEGG", builtin = builtin)
  KO_FILE[,1]<-as.vector(KO_FILE[,1])
  KO_FILE[,2]<-as.vector(KO_FILE[,2])
  KO_FILE$Annot<-annot[KO_FILE[,2],"annotation"]
  colnames(KO_FILE)[1]<-"GeneID"
  return(KO_FILE)
}

##' Download database from Msigdb and prepare for enrichment analysis
##' @name buildMSIGDB
##' @importFrom msigdbr msigdbr
##' @importFrom dplyr filter
##' @importFrom dplyr select
##' @importFrom rlang sym
##' @importFrom magrittr %>%
##' @param species the species for query
##' @param keytype the gene ID type
##' @param anntype anntotaion type of  gene set (GO,BP,CC,MF,KEGG,REACTOME,
##' BIOCARTA,HALLMARK)
##' @examples
##' \dontrun{
##' hsamsi<-buildMSIGDB(species = "human", keytype = "SYMBOL", anntype = "GO")
##' }
##' @export
##' @author Kai Guo
buildMSIGDB<-function(species="human",keytype="SYMBOL",anntype="GO"){
  flag <- 0
  ## Lookup table: anntype -> (category, subcategory)
  msig_map <- list(
    HALLMARK     = list(cat = "H",  sub = NULL),
    IMMUNESIGDB  = list(cat = "C7", sub = "IMMUNESIGDB"),
    WIKIPATHWAYS = list(cat = "C2", sub = "CP:WIKIPATHWAYS"),
    CGP          = list(cat = "C2", sub = NULL),
    CP           = list(cat = "C2", sub = NULL),
    KEGG         = list(cat = "C2", sub = "CP:KEGG"),
    REACTOME     = list(cat = "C2", sub = "CP:REACTOME"),
    BIOCARTA     = list(cat = "C2", sub = "CP:BIOCARTA"),
    MIR          = list(cat = "C3", sub = NULL),
    TFT          = list(cat = "C3", sub = NULL),
    CGN          = list(cat = "C4", sub = NULL),
    CM           = list(cat = "C4", sub = NULL),
    GO           = list(cat = "C5", sub = NULL),
    BP           = list(cat = "C5", sub = "GO:BP"),
    CC           = list(cat = "C5", sub = "GO:CC"),
    MF           = list(cat = "C5", sub = "GO:MF")
  )
  if (!anntype %in% names(msig_map)) {
    stop("Unsupported anntype '", anntype,
         "'. Use one of: ", paste(names(msig_map), collapse = ", "))
  }
  info <- msig_map[[anntype]]
  category <- info$cat
  anntypes <- info$sub
  if (anntype == "HALLMARK") anntype <- ""
  mspe<-.getmsig(species)
  if(is.null(mspe)){
    stop(cat("can't find support species!\n"))
  }
  cat("Downloading msigdb datasets ...\n")
  res <- msigdbr(species = mspe)
  ## Handle column name changes in msigdbr >= 10.0
  ## Old: gs_cat, gs_subcat, gene_symbol, entrez_gene
  ## New: gs_collection, gs_subcollection, gene_symbol, ncbi_gene
  col_cat    <- if ("gs_cat"    %in% colnames(res)) "gs_cat"    else "gs_collection"
  col_subcat <- if ("gs_subcat" %in% colnames(res)) "gs_subcat" else "gs_subcollection"
  col_entrez <- if ("entrez_gene" %in% colnames(res)) "entrez_gene" else "ncbi_gene"
  if(keytype == "SYMBOL"){
    key = "gene_symbol"
  }else if(keytype == "ENTREZID"){
    key = col_entrez
  }else{
    key = col_entrez
    flag = 1
  }
  res <- res[res[[col_cat]] == category, ]
  if(!is.null(anntypes)){
    res <- res[res[[col_subcat]] == anntypes, ]
  }
  if(anntype%in%c("GO", "BP", "CC", "MF")){
    res <- res[, c(key, "gs_exact_source", "gs_name")]
    res <- as.data.frame(res)
    colnames(res) <- c("GeneID","GOALL","Annot")
    res$Annot <- sub('.*@','',sub('_','@',res$Annot))
  }else if(anntype == "KEGG"){
    res <- res[,c(key,"gs_exact_source","gs_name")]
    res<-as.data.frame(res)
    colnames(res)<-c("GeneID","PATH","Annot")
    res$PATH<-substr(res$PATH,4,8)
  }else{
    res<-res[,c(key,"gs_name")]
    res<-as.data.frame(res)
    colnames(res)<-c("GeneID","Term")
    res$Term<-sub('.*@','',sub('_','@',res$Term))
    res$Annot<-res[,2]
  }
  if(flag==1){
    res[,1]<-idconvert(species, keys=res[,1], fkeytype="ENTREZID",
                      tkeytype=keytype)
    res<-na.omit(res)
  }
  result<-new("Annot",
              species = species,
              anntype = anntype,
              keytype = keytype,
              annot = res)
  return(result)
}

#' make Reactome annotation data function
#' @importFrom AnnotationDbi as.list
#' @importFrom dplyr left_join
#' @param species you can check the supported species by using showAvailableRO
#' @param keytype key type export
#' @author Kai Guo
.makeROdata<-function(species="human",keytype="SYMBOL"){
  dbname<-.getrodbname(species=species);
  if (!requireNamespace("reactome.db", quietly = TRUE)) {
    stop("Package 'reactome.db' is required. Install it with:\n",
         "  BiocManager::install('reactome.db')", call. = FALSE)
  }
  dbname=sapply(strsplit(dbname,"_"),'[[',1)
  lhs<-as.list(reactomePATHNAME2ID)
  lhs<-lhs[grep(dbname,names(lhs))]
  roid<-as.list(reactomePATHID2EXTID)[unique(as.vector(unlist(lhs)))]
  roid<-lapply(roid, function(x)unique(x))
  roid<-data.frame("GeneID"=unlist(roid),"Term"=rep(names(roid),times=lapply(roid, length)),row.names=NULL)
  ll<-lapply(lhs,function(x)unique(x))
  roan<-data.frame("Term"=unlist(ll),"Annot"=rep(names(ll),times=lapply(ll,length)),row.names=NULL)
  roan$Annot<-sub('.*@ ','',sub(':','@',roan$Annot))
  res<-left_join(roid,roan,by=c("Term"="Term"))
  if(keytype!="ENTREZID"){
    keys = idconvert(species = species,keys = res$GeneID,fkeytype = "ENTREZID",tkeytype = keytype)
    res$GeneID = keys
    res <- na.omit(res)
  }
  return(res[,c(1,3)])
}

###
#' make KEGG module annotation data function
#' @importFrom AnnotationDbi keys
#' @importFrom KEGGREST keggLink
#' @param species you can check the support species by using showData()
#' @param keytype the gene ID type
#' @param builtin use KEGG built-in annotation or not (default TRUE)
#' @author Kai Guo
.makeKOMdata<-function(species="human",keytype="ENTREZID",builtin=TRUE){
    spe=.getspeices(species)
    tmp<-keggLink("module",spe)
    tmp<-substr(tmp,8,13)
    names(tmp)<-sub('.*:','',names(tmp))
    tmp<-vec_to_df(tmp,name=c(keytype,"Module"))
    if(keytype!="ENTREZID"){
      tmp[,1]<-idconvert(species,keys=tmp[,1],fkeytype = "ENTREZID",tkeytype = keytype)
      tmp<-na.omit(tmp)
    }
    MO_FILE=tmp
    annot<-.get_kgm.data()
    MO_FILE[,1]<-as.vector(MO_FILE[,1])
    MO_FILE[,2]<-as.vector(MO_FILE[,2])
    MO_FILE$Annot<-annot[MO_FILE[,2],"annotation"]
    colnames(MO_FILE)[1]<-"GeneID"
    return(MO_FILE)
}



