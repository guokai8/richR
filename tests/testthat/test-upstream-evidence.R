# Integrated evidence contracts preserve the upstream statistical policy.
evidence_annotation <- function() {
  rbind(data.frame(gene=paste0('g',1:10), term='GOOD'),
        data.frame(gene=paste0('g',1:2), term='SMALL'),
        data.frame(gene=paste0('g',1:600), term='BIG'),
        data.frame(gene=paste0('g',601:1000), term='ZERO'))
}
evidence_result <- function(...) enrich(paste0('g',1:4), evidence_annotation(), ...)

test_that('all eight integrated APIs resolve', {
  for (n in c('testedTerms','bindRich','writeRich','readRich','richThemes',
              'themeGenes','richSensitivity','richReport'))
    expect_true(n %in% getNamespaceExports('richR'), info=n)
})

test_that('ORA captures the actual old family before size and strict reporting filters', {
  x <- evidence_result(keepRich=FALSE)
  tt <- testedTerms(x)
  expect_setequal(tt$Annot,c('GOOD','SMALL','BIG'))
  expect_identical(tt$reported,tt$Annot %in% x@result$Annot)
  expect_equal(tt$Pvalue,phyper(tt$Significant-1,tt$Annotated,1000-tt$Annotated,4,lower.tail=FALSE))
  expect_equal(tt$Padj,p.adjust(tt$Pvalue,'BH'))
  m <- enrichmentInfo(x)
  expect_identical(m$testing$zero_overlap_count,0L)
  expect_identical(m$testing$tested_count,3L)
  expect_identical(m$testing$family_policy,'positive_overlap_before_reporting_filters')
  expect_identical(x@result$Annot,'GOOD')
  expect_identical(evidence_result(pvalue=1e-12,padj=1,keepRich=FALSE)@result$Annot,'GOOD')
  expect_setequal(evidence_result(keepRich=TRUE)@result$Annot,c('GOOD','SMALL'))
})

test_that('mapping records survive analysis coercion and caller inputs are unchanged', {
  a <- evidence_annotation()
  mp <- data.frame(from=c('a','b','c'), to=c('g1','g2','g3'))
  q <- mapGenes(c('a','b','a','X','c'),mp,from='from',to='to')
  before <- q
  r <- enrich(q,a,pvalue=1,minSize=1,minGSSize=1)
  expect_identical(q,before)
  expect_identical(enrichmentInfo(r)$input$mapping,attr(q,'richR_mapping'))
  au <- mappingAudit(r)
  expect_identical(sum(au$role=='query'),5L)
  expect_equal(enrichmentInfo(r)$input$occurrences,c(1L,1L,1L))
  expect_identical(attr(q,'richR_mapping')$decisions$input_id,c('a','b','a','X','c'))
  expect_true(is.character(enrichmentInfo(r)$provenance$annotation$fingerprint$membership))
})

test_that('binding preserves failed and unavailable evidence without fabricated values', {
  r <- evidence_result(keepRich=FALSE)
  old <- r
  old@tested <- data.frame()
  attr(old,'richR_metadata') <- NULL
  b <- bindRich(list(ok=r,old=old,failed=NULL,empty=evidence_result(pvalue=1e-30)))
  expect_true(all(b$execution[b$group=='failed']=='failed'))
  expect_true(all(is.na(b$Pvalue[b$group=='failed'])))
  expect_true('unavailable' %in% b$state[b$group=='old'])
  expect_true('tested_not_reported' %in% b$state[b$group=='ok'])
})

test_that('typed IO roundtrips original tables and rejects corruption', {
  r <- evidence_result()
  r@result$note <- c('NA',NA_character_)[seq_len(nrow(r@result))]
  r@detail$flag <- rep(NA_integer_,nrow(r@detail))
  p <- tempfile(); on.exit(unlink(p,recursive=TRUE))
  writeRich(r,p)
  z <- readRich(p)
  expect_identical(z@result,r@result)
  expect_identical(z@tested,r@tested)
  expect_identical(z@detail,r@detail)
  expect_identical(enrichmentInfo(z),enrichmentInfo(r))
  cat('tamper\n',file=file.path(p,'result.tsv'),append=TRUE)
  expect_error(readRich(p),'checksum|hash|inconsistent|fields found')
})

test_that('GSEA evidence and typed IO preserve legacy schema and spaced identifiers', {
  genes <- paste('gene',1:60)
  x <- setNames(seq(-3,3,length.out=60),genes)
  a <- rbind(data.frame(gene=genes[1:12],term='LOW'),
             data.frame(gene=genes[49:60],term='HIGH'),
             data.frame(gene=genes[20:35],term='MIDDLE'))
  args <- list(x=x,object=a,seed=7,nproc=1,eps=0,minSize=2,maxSize=50,pvalue=1)
  joined <- do.call(richGSEA,c(args,list(table=TRUE,sep='|')))
  raw <- do.call(richGSEA,c(args,list(table=FALSE,sep='|')))
  expect_false('Term' %in% names(joined@result))
  expect_identical(joined@result[setdiff(names(joined@result),'leadingEdge')],raw@result[setdiff(names(raw@result),'leadingEdge')])
  expect_identical(strsplit(joined@result$leadingEdge,'|',fixed=TRUE),unname(raw@result$leadingEdge))
  expect_true(all(unlist(raw@result$leadingEdge) %in% names(x)))
  tt <- testedTerms(joined)
  expect_identical(tt$reported,tt$pathway %in% joined@result$pathway)
  expect_equal(tt$padj,p.adjust(tt$pval,'BH'))
  p <- tempfile(); on.exit(unlink(p,recursive=TRUE))
  writeRich(joined,p); z <- readRich(p)
  expect_identical(z@result,joined@result)
  expect_identical(z@input,joined@input)
  expect_identical(bindRich(list(g=z))$Term,z@result$pathway)
  th <- richThemes(list(g=z))
  expect_true(all(themeGenes(th)$gene %in% genes))
})

test_that('sensitivity respects original family and independent backgrounds', {
  a <- evidence_annotation(); q <- paste0('g',1:4); before <- q
  g <- richSensitivity(q,a,backgrounds=list(all=paste0('g',1:1000),small=paste0('g',1:700)),
                       maxGSSize=c(20,800),minGSSize=1,keepRich=FALSE,pvalue=1)
  expect_identical(q,before)
  expect_identical(nrow(g$cells),4L)
  for (r in g$results) {
    tt <- testedTerms(r); m <- enrichmentInfo(r)
    expect_equal(tt$Pvalue,phyper(tt$Significant-1,tt$Annotated,m$testing$effective_N-tt$Annotated,m$testing$effective_n,lower.tail=FALSE))
    expect_equal(tt$Padj,p.adjust(tt$Pvalue,'BH'))
  }
  expect_identical(g$results[[1]]@tested$Padj,g$results[[2]]@tested$Padj)
})

test_that('measured universe is explicit and preserves the annotated default', {
  a <- data.frame(gene=c('a','b','c'),term=c('T','T','OTHER'))
  q <- c('a','outside'); u <- c('a','b','c','outside','other')
  plain <- enrich(q,a,universe=u,pvalue=1,minGSSize=1,minSize=1)
  measured <- enrich(q,a,universe=u,universe_policy='measured',pvalue=1,minGSSize=1,minSize=1)
  expect_identical(enrichmentInfo(plain)$testing$effective_N,3L)
  expect_identical(enrichmentInfo(measured)$testing$effective_N,5L)
  expect_identical(enrichmentInfo(measured)$testing$effective_n,2L)
  expect_equal(testedTerms(measured)$Pvalue,phyper(0,2,3,2,lower.tail=FALSE))
  expect_error(enrich(q,a,universe_policy='measured'),'explicit|universe')
})

test_that('report-only selection includes exact P=1 and records real rendering', {
  r <- evidence_result(pvalue=1,minGSSize=1,minSize=1)
  r@result$Pvalue[] <- 1
  r@result$Padj[] <- 1
  before <- r
  draw <- richR:::.report_draw(r,FALSE,1)
  expect_identical(draw$k,1L)
  expect_identical(length(draw$ids),1L)
  expect_gt(nchar(draw$png),100L)
  p <- tempfile(); on.exit(unlink(p,recursive=TRUE))
  richReport(list(one=r),p,figure.top=1)
  expect_true(all(file.exists(file.path(p,c('report.html','manifest.json','rerun.R')))))
  expect_silent(parse(file.path(p,'rerun.R')))
  expect_identical(r,before)
  expect_identical(nrow(suppressMessages(richR:::ggNES_internal(data.frame(pathway='T',pval=1,padj=1,NES=1),pvalue=1))$data),NULL)
})

test_that('binding malformed GSEA identifiers keeps contextual validation', {
  x <- new('GSEAResult',result=data.frame(pval=.1,NES=1),sep=',')
  expect_error(bindRich(list(broken=x)),"broken.*lacks column|lacks column.*pathway")
})

test_that('the old positional Annot enrich signature retains builtin position', {
  a <- data.frame(GeneID=paste0('G',1:20),Term=rep(c('A','B'),each=10),Annot=rep(c('A','B'),each=10))
  obj <- new('Annot',species='human',anntype='custom',keytype='SYMBOL',annot=a)
  r <- enrich(c('G1','G2'),obj,'',1,NULL,NULL,'',NULL,1,500,1,500,FALSE,'BH',',',NULL,FALSE)
  expect_identical(nrow(r@result),1L)
})

test_that('versioned table checksums detect valid-shaped tampering', {
  r <- evidence_result(); p <- tempfile(); on.exit(unlink(p,recursive=TRUE))
  writeRich(r,p)
  m <- jsonlite::fromJSON(file.path(p,'manifest.json'),simplifyVector=FALSE)
  expect_identical(m$format_version,2L)
  result <- which(vapply(m$tables,function(e) identical(e$slot,'result'),logical(1)))
  expect_match(m$tables[[result]]$sha256,'^[a-f0-9]{64}$')
  lines <- readLines(file.path(p,'result.tsv'))
  # Alter one label without changing columns, rows or any numeric syntax.
  fields <- strsplit(lines[2],'\t',fixed=TRUE)[[1]]
  fields[2] <- 'changed-label'; lines[2] <- paste(fields,collapse='\t')
  writeLines(lines,file.path(p,'result.tsv'),useBytes=TRUE)
  expect_error(readRich(p),'checksum.*result.tsv|result.tsv.*checksum')
})
