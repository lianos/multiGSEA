context("calculateIndividualLogFC")

tt2dt <- function(x) {
  if (is(x, 'TopTags')) {
    onames <- c('PValue', 'FDR')
  } else {
    onames <- c('P.Value', 'adj.P.Val')
  }
  x <- as.data.frame(x)
  x$feature_id <- rownames(x)
  data.table::setnames(x, onames, c('pval', 'padj'))
}

## TODO: Test that logFC's are calculated correctly when using contrast
##       vectors
test_that("logFC's calculated from contrast vectors are correct", {
  y <- exampleExpressionSet('tumor-subtype', do.voom=FALSE)
  d0 <- y$design
  di <- model.matrix(~ PAM50subtype, data = y$samples)
  colnames(di) <- sub('PAM50subtype', '', colnames(di))

  ## limma ----------------------------------------------------------
  ## with intercept
  vmi <- limma::voom(y, di)
  fit <- limma::lmFit(vmi, vmi$design)
  e <- limma::eBayes(fit)
  tt.lumA <- tt2dt(limma::topTable(e, 'LumA', number=Inf, sort='none'))
  tt.her2 <- tt2dt(limma::topTable(e, 'Her2', number=Inf, sort='none'))

  ## no intercept
  cm <- limma::makeContrasts(her2.vs.basal=Her2 - Basal,
                             lumA.vs.basal=LumA - Basal,
                             levels=d0)
  vm0 <- limma::voom(y, d0)
  fit0 <- limma::lmFit(vm0, vm0$design)
  e0 <- limma::eBayes(limma::contrasts.fit(fit0, cm))

  tt0.lumA <- tt2dt(limma::topTable(e0, 'lumA.vs.basal', number=Inf, sort='none'))
  tt0.her2 <- tt2dt(limma::topTable(e0, 'her2.vs.basal', number=Inf, sort='none'))

  ## Checks results from analysis w/ and w/o intercepts
  expect_equal(tt.lumA, tt0.lumA)
  expect_equal(tt.her2, tt0.her2)

  ## logFC via multiGSEA codepath ----------------------------------------------
  ## 1. Using coef from design matrix
  my.tt.lumA <- calculateIndividualLogFC(vmi, vmi$design, 'LumA')
  my.tt.her2 <- calculateIndividualLogFC(vmi, vmi$design, 'Her2')
  expect_equal(tt.lumA, my.tt.lumA[, names(tt.lumA)], check.attributes=FALSE)
  expect_equal(tt.her2, my.tt.her2[, names(tt.her2)], check.attributes=FALSE)

  ## 2. Using a contrast vector
  my.tt0.lumA <- calculateIndividualLogFC(vm0, d0, cm[, 'lumA.vs.basal'])
  my.tt0.her2 <- calculateIndividualLogFC(vm0, d0, cm[, 'her2.vs.basal'])
  expect_equal(tt.lumA, my.tt0.lumA[, names(tt.lumA)], check.attributes=FALSE)
  expect_equal(tt.her2, my.tt0.her2[, names(tt.her2)], check.attributes=FALSE)
})

test_that("treat pvalues are legit", {
  lfc <- log2(1.25)
  y <- exampleExpressionSet(do.voom = FALSE)
  y <- edgeR::estimateDisp(y, y$design, robust=TRUE)

  vm <- limma::voom(y, y$design)

  ## limma/voom ----------------------------------------------------------------
  fit <- limma::lmFit(vm, vm$design)
  e <- limma::treat(fit, lfc=lfc)
  tt <- tt2dt(limma::topTreat(e, 'tumor', number=Inf, sort='none'))

  xx <- calculateIndividualLogFC(vm, vm$design, 'tumor', treat.lfc=lfc)

  expect_equal(tt$feature_id, xx$feature_id, info="voom")
  expect_equal(xx$logFC, tt$logFC, info="voom")
  expect_equal(xx$pval, tt$pval, info="voom")

  ## edgeR
  yfit <- edgeR::glmQLFit(y, y$design, robust = TRUE)
  res <- edgeR::glmTreat(yfit, coef='tumor', lfc = lfc)

  et <- tt2dt(edgeR::topTags(res, Inf, sort.by='none'))

  yy <- calculateIndividualLogFC(y, y$design, 'tumor', treat.lfc = lfc)

  expect_equal(et$feature_id, yy$feature_id, info="edgeR")
  expect_equal(yy$logFC, et$logFC, info="edgeR")
  expect_equal(yy$pval, et$pval, info="edgeR")
})

test_that("edgeR's glmLRT or QLF are used when asked", {
  y <- exampleExpressionSet(do.voom = FALSE)
  y <- edgeR::estimateDisp(y, y$design, robust = TRUE)
  gdb <- exampleGeneSetDb()

  ex.qlf <- glmQLFit(y, y$design, robust = TRUE) %>%
    glmQLFTest %>%
    topTags(n = Inf, sort.by = "none") %>%
    as.data.frame
  mgq <- multiGSEA(gdb, y, y$design, use.qlf = TRUE)
  expect_equal(logFC(mgq)$pval, ex.qlf$PValue, info = "QLF")

  ex.lrt <- glmFit(y, y$design) %>%
    glmLRT %>%
    topTags(n = Inf, sort.by = "none") %>%
    as.data.frame
  mgl <- multiGSEA(gdb, y, y$design, use.qlf = FALSE)
  expect_equal(logFC(mgl)$pval, ex.lrt$PValue, info = "LRT")

  # Pvalues from QLF and LRT should not be the same
  expect_false(isTRUE(all.equal(logFC(mgq)$pval, logFC(mgl)$pval)))
})

test_that("calculateIndividualLogFC supports basic ANOVA", {
  y <- exampleExpressionSet('tumor-subtype', do.voom=FALSE)
  d0 <- y$design
  di <- model.matrix(~ PAM50subtype, data = y$samples)
  vm <- voom(y, di)
  anova.res <- lmFit(vm, vm$design) %>%
    eBayes() %>%
    topTable(coef = 2:3, n = Inf)

  lfc <- calculateIndividualLogFC(vm, vm$design, 2:3)

  cmp <- merge(
    anova.res[, c("symbol", "entrez_id", "F", "P.Value")],
    lfc[, c("entrez_id", "F", "pval")], by = "entrez_id")
  expect_equal(cmp$F.x, cmp$F.y)
  expect_equal(cmp$P.Value, cmp$pval)
})

test_that("use explicit observation weights", {
  vm <- exampleExpressionSet(do.voom = TRUE)

  # we have already tested that this codepath produces the "standard"
  # limma::voom result
  vm.res <- calculateIndividualLogFC(vm, vm$design, "tumor")

  # explicit matrix and weights
  w.res <- calculateIndividualLogFC(vm$E, vm$design, "tumor",
                                    weights = vm$weights)
  expect_equal(w.res[["feature_id"]], vm.res[["feature_id"]])
  expect_equal(w.res[["logFC"]], vm.res[["logFC"]])
  expect_equal(w.res[["pval"]], vm.res[["pval"]])

  # standard voom w/ custom weights. Order of output is the same, but
  # everything (logFC, t-stats, etc.) are different
  set.seed(0xFEED)
  W <- runif(length(vm$E), min = min(vm$weights), max = max(vm$weights))
  W <- matrix(W, nrow = nrow(vm))
  x.res <- calculateIndividualLogFC(vm, vm$design, "tumor", weights = W)

  expect_equal(x.res[["feature_id"]], vm.res[["feature_id"]])
  mismatch.logFC <- sign(x.res[["logFC"]]) != sign(vm.res[["logFC"]])
  expect_true(mean(mismatch.logFC) > 0 & mean(mismatch.logFC) < 0.10)
  expect_false(isTRUE(all.equal(x.res[["pval"]], vm.res[["pval"]])))

  # custom weights through limma
  ex.res <- lmFit(vm$E, vm$design, weights = W) %>%
    eBayes() %>%
    topTable(coef = "tumor", n = Inf, sort.by = "none")
  expect_equal(x.res[["feature_id"]], rownames(ex.res))
  expect_equal(x.res[["logFC"]], ex.res[["logFC"]])
  expect_equal(x.res[["pval"]], ex.res[["P.Value"]])
})
