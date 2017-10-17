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

test_that("cameraPR pass through method works like direct call", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsd <- conform(exampleGeneSetDb(), vm)
  gsi <- as.list(gsd, value='x.idx')

  lfc <- logFC(multiGSEA(gsd, vm, vm$design, 'tumor'))

  expected <- cameraPR(setNames(lfc$t, lfc$entrez_id), gsi, inter.gene.cor=0.01)

  mg <- multiGSEA(gsd, vm, vm$design, 'tumor',
                  methods=c('camera', 'cameraPR'),
                  inter.gene.cor=0.01)

  res.pr <- result(mg, 'cameraPR')
  res.pr$key <- encode_gskey(res.pr)
  expect_true(setequal(res.pr$key, rownames(expected)))

  expected <- expected[res.pr$key,,drop=FALSE]
  expect_equal(res.pr$n, expected$NGenes)
  expect_equal(res.pr$pval, expected$PValue)
  expect_equal(res.pr$Direction, expected$Direction)
  expect_equal(res.pr$padj, expected$FDR)
})
