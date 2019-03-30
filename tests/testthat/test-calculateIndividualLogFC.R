context("calculateIndividualLogFC")

tt2dt <- function(x) {
  if (is(x, 'TopTags')) {
    onames <- c('PValue', 'FDR')
  } else {
    onames <- c('P.Value', 'adj.P.Val')
  }
  x <- as.data.frame(x)
  x$featureId <- rownames(x)
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

  expect_equal(tt$featureId, xx$featureId, info="voom")
  expect_equal(xx$logFC, tt$logFC, info="voom")
  expect_equal(xx$pval, tt$pval, info="voom")

  ## edgeR
  yfit <- edgeR::glmQLFit(y, y$design, robust = TRUE)
  res <- edgeR::glmTreat(yfit, coef='tumor', lfc = lfc)

  et <- tt2dt(edgeR::topTags(res, Inf, sort.by='none'))

  yy <- calculateIndividualLogFC(y, y$design, 'tumor', treat.lfc = lfc)

  expect_equal(et$featureId, yy$featureId, info="edgeR")
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
