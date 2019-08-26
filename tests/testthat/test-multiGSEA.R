context("multiGSEA")

test_that("multiGSEA fails on not-full-rank design matrix", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  design <- vm$design

  ## Add `extra` column to design matix which is linear combination
  ## other columns: not full rank
  design <- cbind(design, extra=design[, 1] + design[, 2])

  expect_error(suppressWarnings(multiGSEA(gsd, vm, design)))
})

test_that("multiGSEA wrapper generates same results as individual do.*", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)

  methods <- c('camera', 'cameraPR')
  min.logFC <- log2(1.25)
  max.padj <- 0.10
  mg <- multiGSEA(gsd, vm, vm$design, methods=methods,
                  nrot=250, nsim=500, split.updown=FALSE,
                  feature.min.logFC=min.logFC, feature.max.padj=max.padj)
  lfc <- logFC(mg)
  gsc <- conform(gsd, vm)
  do <- sapply(methods, function(m) {
    fn <- getFunction(paste0('do.', m), where=getNamespace('multiGSEA'))
    fn(gsc, vm, vm$design, nrot=250, nsim=500, split.updown=FALSE,
       feature.min.logFC=min.logFC, feature.max.padj=max.padj)
  }, simplify=FALSE)

  ## Some GSEA results use sampling and their outputs only converge under higher
  ## iterations, which will slow down testing. To avoid that we just use methods
  ## that are deterministic.
  no.random <- c('camera', 'cameraPR')
  for (m in no.random) {
    do.x <- do[[m]]
    mg.x <- mg@results[[m]]
    expect_equal(do.x, mg.x, info=m)
  }
})
