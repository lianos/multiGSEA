context("Expression Utilities")
library(edgeR)

test_that("as.DGEList conversion works", {
  es <- exampleExpressionSet(do.voom=FALSE)

  y.expect <- calcNormFactors(DGEList(exprs(es), genes=fData(es)))
  y.expect$samples <- cbind(
    y.expect$samples,
    pData(es)[, setdiff(names(pData(es)), names(y.expect$samples))])

  ## From ExpressionSet to DGEList
  y.es <- as.DGEList(es, 'exprs')
  expect_equal(y.es$counts, y.expect$counts)
  expect_equal(y.es$genes, y.expect$genes)
  expect_equal(y.es$samples, y.expect$samples)
})

test_that("CPM is like edgeR::cpm", {
  es <- exampleExpressionSet(do.voom=FALSE)
  y <- as.DGEList(es, 'exprs')

  c.expect <- cpm(y, log=TRUE, prior.count=5)
  c.y <- CPM(y, prior.count=5, regularized=FALSE)
  expect_equivalent(c.y, c.expect) ## CPM adds some attribs

  es$norm.factors <- calcNormFactors(exprs(es))
  c.es <- CPM(es, prior.count=5, regularized=FALSE)
  expect_equivalent(c.es, c.expect)
})
