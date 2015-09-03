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

  ## The multiplication of 1 million to bring counts back to "normal" scale
  ## screws the equivalancy in cpm vs CPM(!)
})

test_that("TPM works as expected", {
  es <- exampleExpressionSet(do.voom=FALSE)
  y <- as.DGEList(es, 'exprs')
  cpm <- cpm(y, log=FALSE, prior.count=5) * 1e6
  glen <- multiGSEA:::create.glength.vector(cpm)
  mycpm <- CPM(y, log=FALSE, prior.count=5, regularized=FALSE) * 1e6
  etpm <- apply(cpm, 2, function(x) {
    rate <- log(x) - log(glen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  })

  tpm <- TPM(y, prior.count=5, regularized=FALSE, log=FALSE)
  ## expect_equivalent(log2(tpm), log2(etpm))
})
