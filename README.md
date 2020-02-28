# richR  [![](https://img.shields.io/badge/devel%20version-0.0.1-green.svg)](https://github.com/guokai8/richR) 
## Description
_richR_ provide functions _richGO_, _richKEGG_,and _enrich_ to do enrichment analysis. 
## Installation
```
library(devtools)
install_github("guokai8/richR")
``` 
## Quick tour
```{r}
set.seed(123)   
library(richR)   
# To check if your the current species if supported !!!
showData()   
# Make the GO and KEGG Pathway data for your analysis
# find suitable species name by using showensemble()    
hsa_go<-buildAnnot(species="human",keytype="SYMBOL",anntype = "GO")
hsa_ko<-buildAnnot(species = "human",keytype="SYMBOL", anntype = "KEGG")
```         
____   

### You can make annotation data from ensemble   
* If had bioAnno installed, you can build annotation package with it  
```{r}
# library(bioAnno)
# fromKEGG(species="ath")
# athgo<-buildOwn(dbname="org.ath.eg.db",anntype="GO")  
# athko<-buildOwn(dbname="org.ath.eg.db",anntype="KEGG") 
# Please go over the bioAnno package webpage ("https://github.com/guokai8/bioAnno") to learn more
```   
----

### Simple example for enrichment analysis

```{r,fig.height=6,fig.width=6,fig.align="center",dpi=100}
gene=sample(unique(hsa_go$GeneID),1000)
res<-richGO(gene,godata = hsa_go,ontology ="BP")
head(res)
ggbar(res,top = 20,usePadj = F)
```       

```{r,fig.height=6,fig.width=6,fig.align="center"}
resko<-richKEGG(gene,hsa_ko,pvalue=0.05)
head(resko)
ggdot(resko,top=10,usePadj = F)
```    
You can also get network graphic for any type of enrichment analysis result and also combine different enrichment result
```{r,fig.height=6,fig.width=6,fig.align="center",dpi=100}
ggnetplot(res,top=20)
ggnetwork(res,top=20,weightcut = 0.01)
ggnetmap(list(res,resko),top=50,visNet=TRUE,smooth=FALSE)
```    
