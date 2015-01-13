context("multiGSEA")

test_that("multiGSEA wrapper generates same results as individual do.*", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)

  methods <-c('camera', 'hyperGeometricTest', 'roast', 'geneSetTest')

  mg <- multiGSEA(gsd, vm, vm$design, methods=methods,
                  nrot=250, nsim=500, use.cache=FALSE)

  gsc <- conform(gsd, vm)
  do <- sapply(methods, function(m) {
    fn <- getFunction(paste0('do.', m), where=getNamespace('multiGSEA'))
    fn(gsc, vm, vm$design, nrot=250, nsim=500, use.cache=FALSE)
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
    expect_equal(did.x, res.x, info=paste('values check for', m))
  }
})

if (FALSE) {
test_that("multiGSEA camera+roast run matches default runs of each", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  ## Generate camera and roast results to check against.
  photo <- camera(vm, gsets, vm$design, ncol(vm$design))
  set.seed(123)
  roasted <- mroast(vm, gsets, vm$design, ncol(vm$design), nrot=500,
                    sort='none')

  m <- multiGSEA(vm, gsets.lol, vm$design, methods=c('camera', 'roast'),
                 nrot=500, .seed=123, use.cache=FALSE)

  ## Columns of camera output are NGenes, Correlation, Direction, PValue, FDR
  ## make `my` look like that, and test for equality
  my.photo <- local({
    grab <- c('n', 'Correlation.camera', 'Direction.camera', 'pval.camera',
              'padj.camera')
    out <- m[, grab, with=FALSE]
    setnames(out, names(photo))
    out <- as.data.frame(out)
    rownames(out) <- paste(m$group, m$id, sep='.')
    out[rownames(photo),]
  })

  my.roast <- local({
    ## Temporarily not returning the "mixed" pvals from the roast result
    ## out <- m[, list(n, PropDown.roast, PropUp.roast, Direction.roast,
    ##                 pval.roast, padj.roast,
    ##                 pval.mixed.roast, padj.mixed.roast)]
    ## setnames(out, names(roasted))
    out <- m[, list(n, PropDown.roast, PropUp.roast, Direction.roast,
                    pval.roast, padj.roast)]
    setnames(out, head(names(roasted), -2))
    out <- as.data.frame(out)
    rownames(out) <- paste(m$group, m$id, sep='.')
    out[rownames(roasted),]
  })

  expect_equal(photo, my.photo, info='camera result from multiGSEA')
  expect_equal(roasted[, names(my.roast)],
               my.roast, info='roast result from multiGSEA',
               tolerance=0.1)
})

test_that('multiGSEA camera+hyperG looks reasonable', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  ## Generate camera and roast results to check against.
  photo <- camera(vm, gsets, vm$design, ncol(vm$design))

  m <- multiGSEA(vm, gsets.lol, vm$design,
                 methods=c('camera', 'hyperGeometricTest'),
                 force.reeval=TRUE)

  ## I actually don't know what to really test here, let's just make sure
  ## we get pvals in the right places
  p <- cbind(p.matrix(m, 'nominal'), p.matrix(m, 'adjusted'))
  expect_false(any(is.na(p)))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("feature.id's returned from multiGSEA are correct", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  m <- multiGSEA(vm, gsets.lol, vm$design, methods='camera')

  ## Test that feature.id's are correct as they are returned from multiGSEA
  for (i in 1:nrow(m)) {
    id <- m$id[[i]]
    group <- m$group[[i]]
    ifeatures <- m$feature.id[[i]]

    ## expected features
    efeatures <- intersect(gsets.lol[[group]][[id]], rownames(vm))
    expect_true(setequal(efeatures, ifeatures),
                info=sprintf('%s (intersected entrez)', id))

    ## These are index vectors into `vm` that were previously generated
    remapped <- gsets[[paste(group,id,sep='.')]]
    rfeatures <- rownames(vm)[remapped]
    expect_true(setequal(rfeatures, ifeatures),
                info=sprintf('%s (remapped)', id))
  }
})

test_that("plotting does something reasonable", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  tmp.dir <- sub('/file', '/dir', tempfile())

  ## Get initial results first
  mo <- multiGSEA(vm, gsets.lol, vm$design, methods='camera')

  ## Ensure that the expected number of plots are generated at an FDR
  ## threshold of 0.3
  FDR <- 0.3
  is.sig <- mo$padj.camera <= FDR | mo$padj.by.group.camera <= FDR
  n.expected <- sum(is.sig)

  m <- multiGSEA(vm, gsets.lol, vm$design, methods='camera',
                 outdir=tmp.dir, plots.generate=TRUE,
                 plots.padj.threshold=FDR, score.by='logFC')

  m.plotted <- m[is.na(img.path) == FALSE]
  expect_equal(n.expected, nrow(m.plotted))
  expect_true(setequal(mo$id[is.sig], m.plotted$id))
})
}
