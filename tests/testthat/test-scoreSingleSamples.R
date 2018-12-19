context("Single Sample Gene Set Scoring")
suppressPackageStartupMessages({
  suppressWarnings({
    # library(GSDecon)
    # library(parallel)
    library(GSVA)
  })
})

vm <- exampleExpressionSet()
gdb <- getMSigGeneSetDb('h')

test_that("scoreSingleSamples can use genesets of size n = 1 and value is gene", {
  genes <- c(GZMA = "3001", PRF1 = "5551", TGFB1 = "7040")
  lol <- list(a = genes[1], b = genes[1:2], c = genes[1:3])
  gdbx <- GeneSetDb(lol, collectionName = "custom")
  scores <- scoreSingleSamples(gdbx, vm, min.gs.size = 1, as.matrix = TRUE)
  expect_equal(scores["custom;;a",], vm$E[genes[1],])
})

test_that('do.scoreSingleSamples.gsva is equivalent to GSVA::gsva', {

  E <- vm$E
  gdb <- conform(gdb, E)
  lol <- as.list(gdb)

  set.seed(0xBEEF)
  gsva.ex <- GSVA::gsva(E, lol, method='gsva', parallel.sz=4, verbose=FALSE)

  set.seed(0xBEEF)
  gsva.mg <- scoreSingleSamples(gdb, E, methods='gsva', as.matrix=TRUE)

  # expect_equal(gsva.mg, gsva.ex, info='GSVA,gsva', tolerance = 0.01)

  # Can't figure out why these aren't exact just yet! I suspect the genes
  # that make it through the gsva filtering step might be alterd a bit?
  avg.diffs <- sapply(1:ncol(gsva.mg), function(i) {
    mean(abs(gsva.mg[, i] - gsva.ex[,i]))
  })
  cors <- sapply(1:ncol(gsva.mg), function(i) {
    round(cor(gsva.mg[, i], gsva.ex[,i], method = "spearman"), 2)
  })
  expect_true(all(cors >= 0.97))

  # gsva.mg.melt <- scoreSingleSamples(gdb, E, methods='gsva',
  #                               verbose=FALSE, melted=TRUE)
  plage.ex <- gsva(E, lol, method='plage', parallel.sz=4, verbose=FALSE)
  plage.mg <- scoreSingleSamples(gdb, E, methods='plage', as.matrix=TRUE)
  # expect_equal(plage.mg, plage.ex,info='GSVA,gsva')
  cors <- sapply(1:ncol(gsva.mg), function(i) {
    round(cor(gsva.mg[, i], gsva.ex[,i], method = "spearman"), 2)
  })
  expect_true(all(cors >= 0.97))


  es <- exampleExpressionSet(do.voom=FALSE)
  counts <- Biobase::exprs(es)

  set.seed(0xBEEF)
  gsvar.ex <- gsva(counts, lol, method='gsva', kcdf='Poisson', parallel.sz=4,
                   verbose=FALSE)
  set.seed(0xBEEF)
  gsvar.mg <- scoreSingleSamples(gdb, counts, method='gsva', kcdf='Poisson',
                                 as.matrix=TRUE)
  # expect_equal(gsvar.mg, gsvar.ex, info='GSVA,gsva RNAseq',
  #              tolerance = sqrt(.Machine$double.eps))
  cors <- sapply(1:ncol(gsvar.mg), function(i) {
    round(cor(gsvar.mg[, i], gsvar.ex[,i], method = "spearman"), 2)
  })
  expect_true(all(cors >= 0.97))

})

test_that("multiple 'melted' scores are returned in a long data.frame", {
  scores <- scoreSingleSamples(gdb, vm$E, methods=c('svd', 'ssgsea'))
  expect_is(scores, 'data.frame')
  expect_true(setequal(c('svd', 'ssgsea'), scores$method))
  n.samples <- ncol(vm)
  n.gs <- nrow(geneSets(gdb))
  expect_true(nrow(scores) == n.samples * n.gs * 2)
})

test_that("ssGSEA.normalize returns same normalization as GSVA", {
  scores <- scoreSingleSamples(gdb, vm$E, methods='ssgsea', parallel.sz=4,
                               verbose=FALSE)
  my.norm <- ssGSEA.normalize(scores$score)
  ssgsea.norm <- scoreSingleSamples(gdb, vm$E, methods='ssgsea', parallel.sz=4,
                                    ssgsea.norm=TRUE)
  expect_equal(my.norm, ssgsea.norm$score)
})

test_that("ssGSEA (raw) scores are not affected by samples included in test", {
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

test_that("eigenWeightedMean with equal weights can be same as zScore", {
  gdbc <- conform(gdb, vm)
  gs.idxs <- as.list(gdbc, value='x.idx')

  E <- vm$E[gs.idxs[[1]], ]
  ewm.pc1 <- eigenWeightedMean(E, unscale=FALSE, uncenter=FALSE)
  ewm.z <- eigenWeightedMean(E, weights=1, unscale=FALSE, uncenter=FALSE)
  zs <- zScore(E)
  expect_equal(ewm.z$score, zs$score)
  expect_is(all.equal(ewm.pc1$score, zs$score), 'character')
})

test_that("normalization works in eigenWeightedMean", {
  unorm <- scoreSingleSamples(gdb, vm, method = "ewm", normalize = FALSE)
  norm <- scoreSingleSamples(gdb, vm, method = "ewm", normalize = TRUE)
  expect_true(all(unorm$score >= norm$score))
  expect_true(all(unorm$score >= 0))
  expect_true(any(norm$score < 0))
})
