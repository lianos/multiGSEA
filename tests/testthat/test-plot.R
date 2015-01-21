context("Plotting")

test_that("Plotting individual gene sets works", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), vm)

  mg  <- multiGSEA(gsd, vm, vm$design, methods='camera')

  ## Smallest pvalue from camera
  plot(mg, 'c2', 'LOPEZ_MESOTELIOMA_SURVIVAL_TIME_UP')
  plot(mg, 'c2', 'LOPEZ_MESOTELIOMA_SURVIVAL_TIME_UP', 't')

  plot(mg, 'c2', 'LOPEZ_MESOTELIOMA_SURVIVAL_TIME_UP', type='barcode')
  plot(mg, 'c2', 'LOPEZ_MESOTELIOMA_SURVIVAL_TIME_UP', 't', 'barcode')
})

## TODO: Remove plots that were serialized to HD?
