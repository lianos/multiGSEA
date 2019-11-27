context("geneSetSummaryByGenes")

test_that("geneSetSummaryByGenes,GeneSetDb returns a legit result", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  features <- sample(featureIds(gdb), 10)

  res <- geneSetSummaryByGenes(gdb, features, with.features=TRUE)

  # 1. Ensure that geneset <-> geneset membership is legit. We do this by
  #    manipulating the result into a long form data.table that looks like
  #    what is stored in a GeneSetDb@db. We then filter the gdb.sub@db
  #    table to only include the queried features, then compare the two.
  gdb.sub <- subsetByFeatures(gdb, features)
  db.expect <- gdb.sub@db %>%
    copy %>%
    subset(featureId %in% features) %>%
    setkeyv(c('collection', 'name', 'featureId'))
  db.result <- res %>%
    dplyr::select(collection, name, starts_with('featureId_')) %>%
    reshape2::melt(c('collection', 'name')) %>%
      dplyr::rename(featureId=variable, present=value) %>%
    dplyr::mutate(featureId=sub('featureId_', '', featureId)) %>%
    dplyr::filter(present) %>%
    dplyr::select(-present) %>%
    setDT() %>%
    setkeyv(key(db.expect))
  expect_equal(db.result, db.expect)
})

test_that("geneSetSummaryByGenes,MultiGSEAResult returns a legit result", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  mg <- multiGSEA(gdb, vm, vm$design, ncol(vm$design), method='camera')
  mgdb <- geneSetDb(mg)
  features <- sample(featureIds(mg), 10)

  res <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                               feature.rename=FALSE)

  ## Check that logFC for each feature is accurate
  lfc <- reshape2::melt(as.matrix(res[, features, drop=FALSE])) %>%
    dplyr::transmute(featureId=as.character(Var2), logFC=value) %>%
    dplyr::filter(logFC != 0) %>%
    dplyr::distinct(featureId, .keep_all=TRUE) %>%
    dplyr::arrange(featureId)

  lfc.ex <- logFC(mg) %>%
    dplyr::select(featureId, logFC) %>%
    dplyr::filter(featureId %in% features) %>%
    dplyr::arrange(featureId)
  expect_equal(lfc, lfc.ex, check.attributes = FALSE)

  ## check that symbol remapping works, too
  res.s <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                 feature.rename='symbol', as.dt=TRUE)

  lfc.ex <- logFC(mg) %>%
    dplyr::filter(featureId %in% features) %>%
    dplyr::transmute(renamed=ifelse(!is.na(symbol), symbol, paste0('featureId_', featureId)),
                     logFC) %>%
    dplyr::arrange(renamed)
  expect_true(all(lfc.ex$renamed %in% colnames(res.s)))

  lfc.s <- res.s %>%
    as.data.frame() %>%
    dplyr::select(!!lfc.ex$renamed) %>%
    as.matrix() %>%
    reshape2::melt() %>%
    dplyr::transmute(renamed=as.character(Var2), logFC=value) %>%
    dplyr::filter(logFC != 0) %>%
    dplyr::distinct(renamed, .keep_all=TRUE) %>%
    dplyr::arrange(renamed)
  expect_equal(lfc.s, lfc.ex)
})

test_that("geneSetSummary,MultiGSEAResult properly filters significant genesets", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  mg <- multiGSEA(gdb, vm, vm$design, ncol(vm$design), method='camera')
  p.thresh <- 0.20
  camera.sig <- dplyr::filter(result(mg, 'camera'), padj <= p.thresh)

  features <- sample(featureIds(mg), 10)
  res.all <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                   feature.rename='symbol', as.dt=TRUE)
  res.sig <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                   feature.rename='symbol', as.dt=TRUE,
                                   method='camera', max.p=p.thresh)
  expect_true(all(res.sig$name %in% res.all$name))
  expect_true(all(res.sig$name %in% camera.sig$name))
})
