context("multiGSEA")

test_that("multiGSEA camera run matches default camera", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  cres <- camera(vm, gsets, vm$design, ncol(vm$design))

  m <- multiGSEA(vm, gsets.lol, vm$design, methods='camera')

  ## Columns of camera output are NGenes, Correlation, Direction, PValue, FDR
  ## make `my` look like that, and test for equality
  my.camera <- local({
    grab <- c('n', 'Correlation.camera', 'Direction.camera', 'pval.camera',
              'padj.camera')
    out <- m[, grab, with=FALSE]
    setnames(out, c('NGenes', 'Correlation', 'Direction', 'PValue', 'FDR'))
    out <- as.data.frame(out)
    rownames(out) <- paste(m$group, m$id, sep='.')
    out[rownames(cres),]
  })

  expect_equal(cres, my.camera)
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
