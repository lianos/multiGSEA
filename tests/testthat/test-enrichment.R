context("limma::kegga wrapper ('enrichment')")

dfinput <- system.file("extdata", "testdata", "dataframe-input.csv",
                       package = "multiGSEA")
dfinput <- read.csv(dfinput, stringsAsFactors = FALSE)
gdb. <- getMSigGeneSetDb("h", "human", "ensembl")

test_that("biased enrichment is a go", {
  # control for lenght bias
  lbias <- enrichment(gdb., dfinput, selected = "selected",
                      feature.bias = "effective_length",
                      plot.bias = FALSE)
  # control for expression bias
  ebias <- enrichment(gdb., dfinput, selected = "selected",
                      feature.bias = "AveExpr",
                      plot.bias = FALSE)

  # no bias
  nbias <- enrichment(gdb., dfinput, selected = "selected",
                      feature.bias = NULL)


  expect_equal(lbias$Pathway, ebias$Pathway)
  expect_equal(nbias$Pathway, ebias$Pathway)
  expect_equal(lbias$N, ebias$N)
  expect_equal(lbias$all, ebias$all)
  expect_equal(lbias$up, ebias$up)
  expect_equal(lbias$down, ebias$down)

  # expect_true(all(nbias$P.all <= ebias$P.all))
  if (FALSE) {
    plot(-log10(lbias$P.all), -log10(ebias$P.all))
  }
})

test_that("enrichment(feature.bias = NULL) performs unbiased enrichment", {
})
