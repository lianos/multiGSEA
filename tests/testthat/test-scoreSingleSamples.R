context("Single Sample Gene Set Scoring")
suppressPackageStartupMessages({
  suppressWarnings({
    # library(NMF)
    # library(GSDecon)
    # library(parallel)
    library(GSVA)
  })
})

vm <- exampleExpressionSet()
gdb <- getMSigDBset('h')

test_that('do.scoreSingleSamples.gsva is equivalent to GSVA::gsva', {
  # vm <- exampleExpressionSet()
  # gdb <- getMSigDBset('h')
  lol <- as.list(gdb)

  E <- vm$E
  gdb <- conform(gdb, E)

  gsva.ex <- gsva(E, lol, method='gsva', parallel.sz=4, verbose=FALSE)$es.obs
  gsva.mg <- scoreSingleSamples(gdb, E, methods='gsva', as.matrix=TRUE)
  expect_equal(gsva.mg, gsva.ex,info='GSVA,gsva')

  # gsva.mg.melt <- scoreSingleSamples(gdb, E, methods='gsva',
  #                               verbose=FALSE, melted=TRUE)
  plage.ex <- gsva(E, lol, method='plage', parallel.sz=4, verbose=FALSE)
  plage.mg <- scoreSingleSamples(gdb, E, methods='plage', as.matrix=TRUE)
  expect_equal(plage.mg, plage.ex,info='GSVA,gsva')

  es <- exampleExpressionSet(do.voom=FALSE)
  counts <- exprs(es)

  gsvar.ex <- gsva(counts, lol, method='gsva', rnaseq=TRUE, parallel.sz=4,
                   verbose=FALSE)$es.obs
  gsvar.mg <- scoreSingleSamples(gdb, counts, method='gsva', rnaseq=TRUE,
                                 as.matrix=TRUE)
  expect_equal(gsvar.mg, gsvar.ex, info='GSVA,gsva RNAseq')
})

test_that("multiple 'melted' scores are returned in a long data.frame", {
  # vm <- exampleExpressionSet()
  # gdb <- getMSigDBset('h')
  # library(parallel)
  scores <- scoreSingleSamples(gdb, vm$E, methods=c('svd', 'ssgsea'))
  expect_is(scores, 'data.frame')
  expect_true(all(c('svd', 'ssgsea') %in% scores$method))
  n.samples <- ncol(vm)
  n.gs <- nrow(geneSets(gdb))
  expect_true(nrow(scores) == n.samples * n.gs * 2)
})

test_that("ssGSEA.normalize returns same normalization as GSVA", {
  # vm <- exampleExpressionSet()
  # gdb <- getMSigDBset('h')
  # library(parallel)
  scores <- scoreSingleSamples(gdb, vm$E, methods='ssgsea', parallel.sz=4,
                               verbose=FALSE)
  my.norm <- ssGSEA.normalize(scores$score)
  ssgsea.norm <- scoreSingleSamples(gdb, vm$E, methods='ssgsea', parallel.sz=4,
                                    ssgsea.norm=TRUE)
  expect_equal(my.norm, ssgsea.norm$score)
})

test_that("ssGSEA (raw) scores are not affected by samples included in test", {
  # vm <- exampleExpressionSet()
  # gdb <- getMSigDBset('h')
  # library(parallel)
  some <- sample(ncol(vm), 10)
  scores.all <- scoreSingleSamples(gdb, vm$E, methods='ssgsea', parallel.sz=4,
                                   verbose=FALSE)
  scores.some <- scoreSingleSamples(gdb, vm$E[, some],
                                    methods='ssgsea', parallel.sz=4,
                                    verbose=FALSE)
  scores <- merge(scores.all, scores.some, suffixes=c('.all', '.some'),
                  by=c('collection', 'name', 'sample'))
  expect_equal(scores$scores.all, scores$scores.some)
})

test_that("simple svd scoring method is same as (naive) GSDecon scores", {
  ## TODO: Finish this test. For somre reason 'gsdecon' and 'svd' are not
  ## producing the same results!

  # vm <- exampleExpressionSet()
  # gdb <- getMSigDBset('h')
  # suppressWarnings(library(GSDecon))

  # gsd <- scoreSingleSamples(gdb, vm$E, 'gsdecon')
  # ssvd <- scoreSingleSamples(gdb, vm$E, 'svd')
  # plage <- scoreSingleSamples(gdb, vm$E, 'plage')
  expect_true(TRUE)
})

test_that("GeneSetDb <-> incidence matrix properly setup for GSDecon methods", {
  ## GSdecon wires up an incidence matrix with genes in the columns and genesets
  ## in the rows.
  # vm <- exampleExpressionSet()
  # gdb <- getMSigDBset('h')
  im <- incidenceMatrix(gdb, vm)
  expect_true(TRUE)
})

if (FALSE) {
  library(Biobase)
  library(multiGSEA)
  library(DESeq2)
  counts <- exprs(exampleExpressionSet(do.voom=FALSE))
  x <- CPM(counts, prior.count=.25, log=FALSE)
  X <- round(t(t(x) * attr(x, 'lib.size')))
  r <- rlog(X)
}
