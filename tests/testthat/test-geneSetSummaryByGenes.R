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
    dplyr::filter(featureId %in% features) %>%
    setkeyv(c('collection', 'name', 'featureId'))
  db.result <- res %>%
    dplyr::select(collection, name, starts_with('featureId_')) %>%
    melt(c('collection', 'name')) %>%
    setDT %>%
    dplyr::rename(featureId=variable, present=value) %>%
    dplyr::mutate(featureId=sub('featureId_', '', featureId)) %>%
    dplyr::filter(present) %>%
    dplyr::select(-present) %>%
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
    dplyr::transmute(featureId=as.character(Var2), logFC=value) %>%
    dplyr::filter(logFC != 0) %>%
    dplyr::distinct(featureId, .keep_all=TRUE) %>%
    dplyr::arrange(featureId)

  lfc.ex <- logFC(mg) %>%
    dplyr::select(featureId, logFC) %>%
    dplyr::filter(featureId %in% features) %>%
    dplyr::arrange(featureId)
  expect_equal(lfc, lfc.ex)

  ## check that symbol remapping works, too
  res.s <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                 feature.rename='symbol', .external=FALSE)
  
  lfc.ex <- logFC(mg) %>%
    dplyr::filter(featureId %in% features) %>%
    dplyr::transmute(renamed=ifelse(!is.na(symbol), symbol, paste0('featureId_', featureId)),
                     logFC) %>%
    dplyr::arrange(renamed)
  expect_true(all(lfc.ex$renamed %in% colnames(res.s)))
  
  lfc.s <- res.s %>% 
    dplyr::select_(.dots=lfc.ex$renamed) %>% 
    as.matrix %>%
    melt %>%
    dplyr::transmute(renamed=as.character(Var2), logFC=value) %>%
    dplyr::filter(logFC != 0) %>%
    dplyr::distinct(renamed, .keep_all=TRUE) %>%
    dplyr::arrange(renamed)
  expect_equal(lfc.s, lfc.ex)
})

test_that("geneSetSummary,MultiGSEAResult properly filters significant genesets", {
  set.seed(0xBEEF)
  vm <- exampleExpressionSet(do.voom=TRUE)
  gdb <- exampleGeneSetDb()
  mg <- multiGSEA(gdb, vm, vm$design, ncol(vm$design), method='camera')
  p.thresh <- 0.20
  camera.sig <- dplyr::filter(result(mg, 'camera'), padj <= p.thresh)
  
  features <- sample(featureIds(mg), 10)
  res.all <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                   feature.rename='symbol', .external=FALSE)
  res.sig <- geneSetSummaryByGenes(mg, features, with.features=TRUE,
                                   feature.rename='symbol', .external=FALSE,
                                   method='camera', max.p=p.thresh)
  expect_true(all(res.sig$name %in% res.all$name))
  expect_true(all(res.sig$name %in% camera.sig$name))
})

if (FALSE) {
  fn <- '~/web/rutz/BAP1/johnnycache/multiGSEA-NGS429-joint-nointeraction.rds'
  mg <- readRDS(fn)
  fids <- sample(featureIds(mg), 10)
  fids <- c(Jag1='16449', Pdgfa='18590', App='11820', Vcan='13003', Nrp1='18186')
  s <- geneSetSummaryByGenes(mg, fids, feature.rename='symbol')
}
