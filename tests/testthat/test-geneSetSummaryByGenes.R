context("geneSetSummaryByGenes")

test_that("geneSetSummaryByGenes,GeneSetDb returns a legit result", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet(do.voom=TRUE)
  gdb <- exampleGeneSetDb()
  features <- sample(featureIds(gdb), 10)

  res <- geneSetSummaryByGenes(gdb, features, with.features=TRUE)

  ## 1. Ensure that geneset <-> geneset membership is legit. We do this by
  ##    manipulating the result into a long form data.table that looks like
  ##    what is stored in a GeneSetDb@db. We then filter the gdb.sub@db
  ##    table to only include the queried features, then compare the two.
  gdb.sub <- subsetByFeatures(gdb, features)
  db.expect <- gdb.sub@db %>%
    copy %>%
    filter(featureId %in% features) %>%
    setkeyv(c('collection', 'name', 'featureId'))
  db.result <- res %>%
    select(collection, name, starts_with('featureId_')) %>%
    melt(c('collection', 'name')) %>%
    setDT %>%
    rename(featureId=variable, present=value) %>%
    mutate(featureId=sub('featureId_', '', featureId)) %>%
    filter(present) %>%
    select(-present) %>%
    setkeyv(key(db.expect))
  expect_equal(db.result, db.expect)
})

test_that("geneSetSummaryByGenes,MultiGSEAResult returns a legit result", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet(do.voom=TRUE)
  gdb <- exampleGeneSetDb()
  mg <- multiGSEA(gdb, vm, vm$design, ncol(vm$design), method='camera')
  mgdb <- geneSetDb(mg)
  features <- sample(featureIds(mg), 10)

  res <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                               feature.rename=FALSE)

  ## Check that logFC for each feature is accurate
  lfc <- melt(as.matrix(res[, features, drop=FALSE])) %>%
    transmute(featureId=as.character(Var2), logFC=value) %>%
    filter(logFC != 0) %>%
    distinct(featureId, .keep_all=TRUE) %>%
    arrange(featureId)

  lfc.ex <- logFC(mg) %>%
    select(featureId, logFC) %>%
    filter(featureId %in% features) %>%
    arrange(featureId)
  expect_equal(lfc, lfc.ex)

  ## check that symbol remapping works, too
  res.s <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                 feature.rename='symbol', .external=FALSE)
  lfc.s <- res.s %>%
    select(-(collection:N)) %>%
    as.matrix %>%
    melt %>%
    transmute(renamed=as.character(Var2), logFC=value) %>%
    filter(logFC != 0) %>%
    distinct(renamed, .keep_all=TRUE) %>%
    arrange(renamed)

  lfc.ex <- logFC(mg) %>%
    filter(featureId %in% features) %>%
    transmute(renamed=ifelse(!is.na(symbol), symbol, paste0('featureId_', featureId)),
              logFC) %>%
    arrange(renamed)
  expect_equal(lfc.s, lfc.ex)
})

test_that("geneSetSummary,MultiGSEAResult properly filters significant genesets", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet(do.voom=TRUE)
  gdb <- exampleGeneSetDb()
  mg <- multiGSEA(gdb, vm, vm$design, ncol(vm$design), method='camera')
  p.thresh <- 0.20
  camera.sig <- filter(result(mg, 'camera'), padj <= p.thresh)

  res.all <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                   feature.rename='symbol', .external=FALSE)
  res.sig <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                   feature.rename='symbol', .external=FALSE,
                                   method='camera', max.p=p.thresh)
  expect_true(all(res.sig$name %in% res.all$name))
  expect_true(all(res.sig$name %in% camera.sig$name))
})
