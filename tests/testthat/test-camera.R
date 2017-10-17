context("camera")

test_that('camera runs equivalently from do.camera vs direct call', {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsd <- conform(exampleGeneSetDb(), vm)
  gsi <- as.list(gsd, value='x.idx')

  photo <- limma::camera(vm, gsi, vm$design, ncol(vm$design))
  my <- multiGSEA:::do.camera(gsd, vm, vm$design, ncol(vm$design))
  expect_true(setequal(rownames(photo), rownames(my)))
  my <- my[rownames(photo),]

  expect_equal(my, photo, check.attributes=FALSE)
})

test_that("camera result() is decorated correctly and has correct stats", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsd <- conform(exampleGeneSetDb(), vm)
  gsi <- as.list(gsd, value='x.idx')
  mg <- multiGSEA(gsd, vm, vm$design, ncol(vm$design), methods='camera')

  photo <- limma::camera(vm, gsi, vm$design, ncol(vm$design))
  res <- result(mg, 'camera', as.dt=TRUE)
  check <- setDF(res[, list(NGenes=n, Direction, PValue=pval, FDR=padj)])
  rownames(check) <- encode_gskey(res)

  expect_true(setequal(rownames(check), rownames(photo)))
  check <- check[rownames(photo),]
  expect_equal(check, photo, check.attributes=TRUE)
})
