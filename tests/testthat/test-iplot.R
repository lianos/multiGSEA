context("Interactive Plotting")

test_that("rbokeh boxplot with outliers is slow(?)", {
  x <- exampleMultiGSEAResult()
  y <- 'c2'
  j <- 'SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP'
  iplot(x, y, j, type='boxplot')
  iplot(x, y, j, type='boxplot', width=350, height=350)

})
