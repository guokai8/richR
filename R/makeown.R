#' make annotation database using bioAnno results
#' @importFrom AnnotationDbi keys
#' @importFrom dplyr distinct
#' @param dbname database name from bioAnno
#' @param anntype GO or KEGG
#' @param OP BP,CC,MF default use all
#' @param species species name
#' @param keytype gene ID type (e.g. "ENTREZID", "SYMBOL")
#' @examples
#' \dontrun{
#' fromKEGG(species="ath")
#' athgo<-buildOwn(dbname="org.ath.eg.db",anntype="GO")
#' }
#' @author Kai Guo
#' @export
buildOwn<-function(dbname,anntype="GO",OP=NULL,species="",keytype=""){
  provider_package <- dbname
  provider_version <- .package_version_or_na(provider_package)
  if (!requireNamespace(dbname,quietly=TRUE)){
    stop("Please give the package name")
  }
  dbname<-getExportedValue(dbname, dbname)
  if(anntype=="GO"){
    annof<-AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns=c("GOALL","ONTOLOGYALL"))
    colnames(annof)[1]<-"GeneID"
  #  annof<-distinct_(annof,~GeneID, ~GOALL, ~ONTOLOGYALL)
    annof<-annof[,c("GeneID", "GOALL", "ONTOLOGYALL")]
    annot <- getann("GO")
    annof$Annot <- annot[annof[,2],"annotation"]
    if(!is.null(OP)){
      annof<-annof[annof$ONTOLOGYALL==OP,]
    }
  }else if(anntype=="KEGG"){
    annof=AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns="PATH")
    annot<-getann("KEGG")
    annof[,1]<-as.vector(annof[,1])
    annof[,2]<-as.vector(annof[,2])
    annof$Annot<-annot[annof[,2],"annotation"]
  }else if(anntype=="KEGGM"){
    annof=AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns="KEGGM")
    annot<-.get_kgm.data()
    annof[,1]<-as.vector(annof[,1])
    annof[,2]<-as.vector(annof[,2])
    annof$Annot<-annot[annof[,2],"annotation"]
  }else{
    annof=AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns=anntype)
    annof[,1]<-as.vector(annof[,1])
    annof[,2]<-as.vector(annof[,2])
    annof$Annot<-annof[,2]
  }
  normalization <- .normalize_annotation(annof, func_name = "buildOwn")
  annof <- normalization$data
  result<-new("Annot",
              species = species,
              anntype = anntype,
              keytype = keytype,
              annot = annof

  )
  provenance <- .annotation_provenance(
    constructor = "buildOwn", species = species, keytype = keytype,
    anntype = anntype, source_type = "package",
    source_name = provider_package, source_version = NA_character_,
    provider_package = provider_package, provider_version = provider_version,
    constructor_normalization = normalization$counts
  )
  .attach_annotation_metadata(result, normalization, provenance,
                              options = list(ontology = OP))
}


