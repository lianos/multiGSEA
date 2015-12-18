context("multiGSEA")

test_that("multiGSEA fails on not-full-rank design matrix", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  design <- vm$design

  ## Add `extra` column to design matix which is linear combination
  ## other columns: not full rank
  design <- cbind(design, extra=design[, 1] + design[, 2])

  expect_error(suppressWarnings(multiGSEA(gsd, vm, design)))
})

test_that("multiGSEA wrapper generates same results as individual do.*", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)

  methods <-c('camera', 'hyperGeometricTest', 'roast', 'geneSetTest')
  min.logFC <- log2(1.25)
  max.padj <- 0.10
  mg <- multiGSEA(gsd, vm, vm$design, methods=methods,
                  nrot=250, nsim=500, split.updown=FALSE,
                  feature.min.logFC=min.logFC, feature.max.padj=max.padj)

  gsc <- conform(gsd, vm)
  do <- sapply(methods, function(m) {
    fn <- getFunction(paste0('do.', m), where=getNamespace('multiGSEA'))
    fn(gsc, vm, vm$design, nrot=250, nsim=500, split.updown=FALSE,
       feature.min.logFC=min.logFC, feature.max.padj=max.padj)
  }, simplify=FALSE)

  ## Some GSEA results use sampling and their outputs only converge under higher
  ## iterations, which will slow down testing. To avoid that we just use methods
  ## that are deterministic.
  no.random <- c('camera', 'hyperGeometricTest')
  for (m in no.random) {
    did.x <- do[[m]]
    res.x <- mg@results[[m]]
    expect_true(all(names(did.x) %in% names(res.x)),
                info=paste('colnames check for', m))
    res.x <- res.x[, names(did.x), with=FALSE]
    expect_equal(did.x, res.x, info=paste('values check for', m),
                 check.attributes=FALSE)
  }
})
