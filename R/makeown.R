#' make annotation database using bioAnno results
#' @importFrom AnnotationDbi keys
#' @importFrom dplyr distinct
#' @param dbname database name from bioAnno
#' @param anntype GO or KEGG
#' @param OP BP,CC,MF default use all
#' @examples
#' \dontrun{
#' fromKEGG(species="ath")
#' athgo<-buildOwn(dbname="org.ath.eg.db",anntype="GO")
#' }
#' @author Kai Guo
#' @export
buildOwn<-function(dbname,anntype="GO",OP=NULL,species="",keytype=""){
  if (!require(dbname,character.only=TRUE)){
    stop("Please give the package name")
  }else{
    suppressMessages(require(dbname,character.only = T,quietly = T))
  }
  dbname<-eval(parse(text=dbname))
  if(anntype=="GO"){
    annof<-AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns=c("GOALL","ONTOLOGYALL"))
    colnames(annof)[1]<-"GeneID"
  #  annof<-distinct_(annof,~GeneID, ~GOALL, ~ONTOLOGYALL)
    annof<-distinct(annof[,c("GeneID", "GOALL", "ONTOLOGYALL")])
    annot <- getann("GO")
    annof$Annot <- annot[annof[,2],"annotation"]
    if(!is.null(OP)){
      annof<-annof[annof$ONTOLOGYALL==OP,]
    }
  }else if(anntype=="KEGG"){
    annof=AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns="PATH")
    annof<-na.omit(annof)
    annot<-getann("KEGG")
    annof[,1]<-as.vector(annof[,1])
    annof[,2]<-as.vector(annof[,2])
    annof$Annot<-annot[annof[,2],"annotation"]
  }else if(anntype=="KEGGM"){
    annof=AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns="KEGGM")
    annof<-na.omit(annof)
    annot<-.get_kgm.data()
    annof[,1]<-as.vector(annof[,1])
    annof[,2]<-as.vector(annof[,2])
    annof$Annot<-annot[annof[,2],"annotation"]
  }else{
    annof=AnnotationDbi::select(dbname,keys=keys(dbname,keytype=keytype),keytype=keytype,columns=anntype)
    annof<-na.omit(annof)
    annof[,1]<-as.vector(annof[,1])
    annof[,2]<-as.vector(annof[,2])
    annof$Annot<-annof[,2]
  }
  annof <- na.omit(annof)
  result<-new("Annot",
              species = species,
              anntype = anntype,
              keytype = keytype,
              annot = annof

  )
  return(result)
}


