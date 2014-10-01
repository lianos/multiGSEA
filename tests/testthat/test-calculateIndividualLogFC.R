context("Calculating individual logFC")

## TODO: Test that logFC's are calculated correctly when using contrast
##       vectors
test_that("logFC's calculated from contrast vectors are correct", {
  vm <- exampleExpressionSet('tumor-subtype', do.voom=TRUE)
  d <- model.matrix(~ PAM50subtype, vm$targets)
  colnames(d) <- sub('PAM50subtype', '', colnames(d))

  fit <- lmFit(vm, d)
  e <- eBayes(fit)
  tt.lumA <- topTable(e, 'LumA', number=Inf, sort='none')
  tt.her2 <- topTable(e, 'Her2', number=Inf, sort='none')

  d0 <- model.matrix(~ 0 + PAM50subtype, vm$targets)
  colnames(d0) <- sub('PAM50subtype', '', colnames(d0))
  cm <- makeContrasts(her2.vs.basal=Her2 - Basal,
                      lumA.vs.basal=LumA - Basal,
                      levels=d0)
  fit0 <- lmFit(vm, d0)
  e0 <- eBayes(contrasts.fit(fit0, cm))
  tt0.lumA <- topTable(e0, 'lumA.vs.basal', number=Inf, sort='none')
  tt0.her2 <- topTable(e0, 'her2.vs.basal', number=Inf, sort='none')

  ## Check to make sure that I really understand limma ;-)
  expect_equal(tt.lumA, tt0.lumA)
  expect_equal(tt.her2, tt0.her2)

  my.tt.lumA <- calculateIndividualLogFC(vm, d, 'LumA', use.t=FALSE)
  my.tt.her2 <- calculateIndividualLogFC(vm, d, 'Her2', use.t=FALSE)

  my.tt0.lumA <- calculateIndividualLogFC(vm, d0, cm[, 'lumA.vs.basal'],
                                          use.t=FALSE)
  my.tt0.her2 <- calculateIndividualLogFC(vm, d0, cm[, 'her2.vs.basal'],
                                          use.t=FALSE)

  expect_equal(tt.lumA$logFC, unname(my.tt.lumA))
  expect_equal(tt.her2$logFC, unname(my.tt.her2))
  expect_equal(tt.lumA$logFC, unname(my.tt0.lumA))
  expect_equal(tt.her2$logFC, unname(my.tt0.her2))
})
