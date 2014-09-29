context("multiGSEA")

test_that("multiGSEA camera run matches default camera", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsets.lol <- exampleGeneSets('lol')
  gsets <- exampleGeneSets('limma')

  m <- multiGSEA(vm, gsets.lol, vm$design)

  cres <- camera(vm, gsets, vm$design, ncol(vm$design))

  ## Create a data.frame that should be like the one returned from
  ## camera
  camera.df <- data.frame(
    NGenes=as.numeric(m$n),
    Correlation=m$Correlation.camera,
    Direction=m$Direction.camera,
    PValue=m$pval.camera,
    FDR=m$padj.camera,
    stringsAsFactors=FALSE)
  rownames(camera.df) <- paste(m$group, m$id, sep='.')
  camera.df <- camera.df[rownames(cres),]

  expect_equal(cres, camera.df)

  ## Test that feature.id's are correct
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

