context("Errors submitted by users")

test_that("multiGSEA fails with a 1-geneset GeneSetDb", {
  # Submitted by Thomas Sandmann:
  # https://github.com/lianos/multiGSEA/issues/7
  vm <- exampleExpressionSet()
  gdb <- exampleGeneSetDb()
  gdb <- gdb[geneSets(gdb)$name == 'GOZGIT_ESR1_TARGETS_DN']

  test.methods <- c("camera", "roast", "fry", "geneSetTest")
  mg <- multiGSEA(gdb, vm, vm$design, 'tumor',
                  methods=test.methods,
                  ## customize camera parameter:
                  inter.gene.cor=0.04)

  # Calling the show,MultiGSEAResult method is how this was identified,
  # so let's make sure that works, first.
  printed <- tryCatch(capture.output(show(mg)), error = function(e) NULL)
  expect_is(printed, "character")

  # Ensure that each method has an FDR/padj column
  for (method in test.methods) {
    res <- result(mg, method)
    expect_equal(nrow(res), 1L, info = method)
    expect_is(res[["padj"]], "numeric", info = method)
  }
})

test_that("multiGSEA pipeine can handle EList without a genes data.frame", {
  # Thanks to @RussBainer for reporting
  gdb <- exampleGeneSetDb()
  vm <- exampleExpressionSet()

  vm.noG <- vm
  vm.noG$genes <- NULL

  mg <- multiGSEA(gdb, vm, vm$design, "tumor", "camera")

  mg.noG <- expect_warning(
    multiGSEA(gdb, vm.noG, vm.noG$design, "tumor", "camera"),
    "no.*genes", ignore.case = TRUE)

  r <- result(mg)
  r.noG <- result(mg.noG)
  expect_equal(r.noG, r)
})
