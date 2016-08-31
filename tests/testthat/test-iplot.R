context("Interactive Plotting")

test_that("rbokeh boxplot with outliers is slow(?)", {
  if (FALSE) {
    vm <- exampleExpressionSet(do.voom=TRUE)
    gsl <- exampleGeneSets()
    gsd <- conform(GeneSetDb(gsl), vm)
    x <- multiGSEA(gsd, vm, vm$design, methods='camera')
    y <- 'c2'
    j <- 'SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP'
  }

  expect_true(TRUE)
})
