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
  data.table::setkeyv(data.table::setDT(x), 'featureId')
  multiGSEA:::ret.df(x)
}

## TODO: Test that logFC's are calculated correctly when using contrast
##       vectors
test_that("logFC's calculated from contrast vectors are correct", {
  es <- exampleExpressionSet('tumor-subtype', do.voom=FALSE)
  d0 <- es@design
  di <- model.matrix(~ PAM50subtype, data=pData(es))
  colnames(di) <- sub('PAM50subtype', '', colnames(di))

  ## limma ----------------------------------------------------------
  ## with intercept
  vmi <- limma::voom(es, di)
  fit <- limma::lmFit(vmi, vmi$design)
  e <- limma::eBayes(fit)
  tt.lumA <- tt2dt(limma::topTable(e, 'LumA', number=Inf, sort='none'))
  tt.her2 <- tt2dt(limma::topTable(e, 'Her2', number=Inf, sort='none'))

  ## no intercept
  cm <- limma::makeContrasts(her2.vs.basal=Her2 - Basal,
                             lumA.vs.basal=LumA - Basal,
                             levels=d0)
  vm0 <- limma::voom(es, d0)
  fit0 <- limma::lmFit(vm0, vm0$design)
  e0 <- limma::eBayes(limma::contrasts.fit(fit0, cm))

  tt0.lumA <- tt2dt(limma::topTable(e0, 'lumA.vs.basal', number=Inf, sort='none'))
  tt0.her2 <- tt2dt(limma::topTable(e0, 'her2.vs.basal', number=Inf, sort='none'))

  ## Checks results from analysis w/ and w/o intercepts
  expect_equal(tt.lumA, tt0.lumA)
  expect_equal(tt.her2, tt0.her2)

  ## logFC via multiGSEA codepath ----------------------------------------------
  ## 1. Using coef from design matrix
  my.tt.lumA <- calculateIndividualLogFC(vmi, vmi$design, 'LumA',
                                         use.treat=FALSE)
  my.tt.her2 <- calculateIndividualLogFC(vmi, vmi$design, 'Her2',
                                         use.treat=FALSE)
  expect_equal(tt.lumA, my.tt.lumA[, names(tt.lumA)])
  expect_equal(tt.her2, my.tt.her2[, names(tt.her2)])

  ## 2. Using a contrast vector
  my.tt0.lumA <- calculateIndividualLogFC(vm0, d0, cm[, 'lumA.vs.basal'],
                                          use.treat=FALSE)
  my.tt0.her2 <- calculateIndividualLogFC(vm0, d0, cm[, 'her2.vs.basal'],
                                          use.treat=FALSE)
  expect_equal(tt.lumA, my.tt0.lumA[, names(tt.lumA)])
  expect_equal(tt.her2, my.tt0.her2[, names(tt.her2)])
})

test_that("treat pvalues are legit", {
  lfc <- log2(1.25)
  es <- exampleExpressionSet(do.voom=FALSE)
  d <- es@design

  vm <- limma::voom(es, d)
  y <- edgeR::DGEList(exprs(es), group=es$Cancer_Status, genes=fData(es))
  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, d, robust=TRUE)

  ## limma/voom ----------------------------------------------------------------
  fit <- limma::lmFit(vm, vm$design)
  e <- limma::treat(fit, lfc=lfc)
  tt <- tt2dt(limma::topTreat(e, 'tumor', number=Inf, sort='none'))

  xx <- calculateIndividualLogFC(vm, d, 'tumor', use.treat=TRUE, treat.lfc=lfc)

  expect_equal(tt$featureId, xx$featureId, info="voom")
  expect_equal(xx$logFC, tt$logFC, info="voom")
  expect_equal(xx$pval, tt$pval, info="voom")

  ## edgeR
  yfit <- edgeR::glmQLFit(y, d, robust=TRUE)
  res <- edgeR::glmTreat(yfit, coef='tumor', lfc=lfc)
  ## TODO: Fix testing error in line above:
  ## calling library(edgeR) fixes this here as well as the embedded call to
  ## this within multiGSEA::calculateIndividualLogFC
  ## Error:
  ##   1: edgeR::glmTreat(yfit, coef = "tumor", lfc = lfc)
  ##   (subscript) logical subscript too long
  ##   2: glmfit[i, ]
  ##   3: `[.DGEGLM`(glmfit, i, )
  ##   4: subsetListOfArrays(object, i, j, IJ = IJ, IX = IX, I = I, JX = JX)

  et <- tt2dt(edgeR::topTags(res, Inf, sort.by='none'))

  yy <- calculateIndividualLogFC(y, d, 'tumor', use.treat=TRUE, treat.lfc=lfc)

  expect_equal(et$featureId, yy$featureId, info="edgeR")
  expect_equal(yy$logFC, et$logFC, info="edgeR")
  expect_equal(yy$pval, et$pval, info="edgeR")
})
