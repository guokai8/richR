# Independent-audit regressions: actual report renderer selection and provenance.
audit_plot_fixture <- function(gsea=FALSE) {
  ids <- LETTERS[1:7]
  if (gsea) {
    tab <- data.frame(pathway=ids, pval=c(1,.01,1,1,NA,1,1),
                      padj=c(1,.02,1,1,1,1,1),NES=c(4,2,Inf,NA,-Inf,2,1))
    tab$pathway[6] <- NA_character_
    new('GSEAResult',result=tab,sep=',')
  } else {
    tab <- data.frame(Annot=ids,Term=ids,Pvalue=c(1,.01,1,1,NA,1,1),
                      Padj=c(1,.02,NA,1,1,1,1),Annotated=c(10,10,10,0,10,10,10),
                      Significant=c(2,3,2,2,Inf,2,2),GeneID='gene')
    tab$Term[6] <- NA_character_
    new('richResult',result=tab,sep=',')
  }
}

test_that('ORA and GSEA reports separate genuine drawing omissions from top truncation', {
  for (gsea in c(FALSE,TRUE)) {
    obj <- audit_plot_fixture(gsea)
    before <- serialize(obj,NULL)
    for (top in c(Inf,1)) {
      ids <- if(is.infinite(top)) c('A','B','G') else 'A'
      f <- richR:::.report_figure(obj,'audit',gsea,TRUE,top)
      expect_identical(f$record$ids,ids)
      expect_identical(f$record$k,as.integer(length(ids)))
      expect_identical(f$record$n,7L)
      expect_identical(f$record$eligible,3L)
      expect_identical(f$record$omitted,4L)
      expect_identical(f$record$truncated,as.integer(3-length(ids)))
      expect_match(f$caption,'4 omitted for missing, blank or non-finite drawing values',fixed=TRUE)
      expect_match(f$caption,paste0(3-length(ids),' truncated by top'),fixed=TRUE)
      expect_match(f$record$omission_reason,'missing, blank or non-finite',fixed=TRUE)
      # Independently inspect the plot's actual built layer, not only the signal.
      signal <- NULL
      plot <- withCallingHandlers(
        if(gsea) richR:::ggNES_internal(obj@result,top=top,.report_all=TRUE)
        else richR:::ggdot_internal(obj@result,top=top,.report_all=TRUE),
        richR_plot_selection=function(c) {signal <<- c;invokeRestart('muffleMessage')})
      id <- if(gsea) 'pathway' else 'Annot'
      expect_identical(as.character(plot$data[[id]]),ids)
      built <- ggplot2::ggplot_build(plot)$data[[1]]
      expect_identical(nrow(built),as.integer(length(ids)))
      expect_true(all(is.finite(built$x) & is.finite(built$y)))
      expect_identical(signal$ids,ids)
    }
    # Ordinary public paths retain strict P < 1 and select only B.
    p <- if(gsea) richNES(obj,pvalue=1,padj=NULL,top=Inf)
      else richDot(obj,pvalue=1,padj=NULL,top=Inf)
    id <- if(gsea) 'Term' else 'Annot'
    # Legacy subsetting retains the NA placeholder from an NA raw p-value.
    expect_identical(as.character(p$data[[id]]),c('B',NA_character_))
    expect_identical(serialize(obj,NULL),before)
  }
})

test_that('full annotation fingerprints reflect actual GO labels without changing membership or statistics', {
  a <- data.frame(GeneID=paste0('G',1:20),GOALL=rep(c('GO:0000001','GO:0000002'),each=10),
                  Annot=rep(c('label A','label B'),each=10),ONTOLOGYALL='BP')
  run <- function(d) richGO(c('G1','G2'),d,pvalue=1,padj=NULL,minSize=1,minGSSize=1,keepRich=FALSE)
  x <- run(a);a$Annot[a$Annot=='label A'] <- 'changed A';y <- run(a)
  f <- function(r) enrichmentInfo(r)$provenance$annotation$fingerprint
  expect_false(identical(x@result$Term,y@result$Term))
  expect_identical(f(x)$membership,f(y)$membership)
  expect_false(identical(f(x)$full,f(y)$full))
  cols <- setdiff(names(x@result),'Term')
  expect_identical(x@result[cols],y@result[cols])
  expect_identical(x@tested[setdiff(names(x@tested),'Term')],y@tested[setdiff(names(y@tested),'Term')])
  expect_identical(enrichmentInfo(x)$testing,enrichmentInfo(y)$testing)
  blank <- a;blank$Annot[blank$GOALL=='GO:0000001'] <- ' '
  missing <- a;missing$Annot[missing$GOALL=='GO:0000001'] <- NA_character_
  fallback <- a;fallback$Annot[fallback$GOALL=='GO:0000001'] <- 'GO:0000001'
  b <- run(blank);m <- run(missing);d <- run(fallback)
  expect_identical(b@result$Term,'GO:0000001')
  expect_identical(f(b)$full,f(m)$full)
  expect_identical(f(b)$full,f(d)$full)
  expect_identical(f(b)$membership,f(x)$membership)
})
