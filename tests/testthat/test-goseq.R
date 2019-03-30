context("goseq")

test_that("multiGSEA(method='goseq') requires valid feature.bias vector", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsd <- conform(gsd, vm)

  mg <- multiGSEA(gsd, vm, vm$design)
  lfc <- logFC(mg)
  selected <- subset(lfc, significant)$featureId
  universe <- rownames(vm)

  expect_error({
    suppressWarnings({
      multiGSEA(gsd, vm, vm$design, methods='goseq', split.updown=FALSE)
    })
  })
})


test_that("internal goseq mimics goseq package", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsd <- conform(gsd, vm)

  ## Identify differentially expressed genes
  mg <- multiGSEA(gsd, vm, vm$design)
  lfc <- logFC(mg)
  selected <- subset(lfc, significant)$featureId
  universe <- rownames(vm)
  mylens <- setNames(vm$genes$size, rownames(vm))
  degenes <- setNames(integer(length(universe)), universe)
  degenes[selected] <- 1L

  ## Run internal version of goseq
  my.res <- suppressWarnings({
    multiGSEA::goseq(gsd, selected, universe, mylens, method='Wallenius',
                     use_genes_without_cat=TRUE, .pipelined=TRUE)
  })
  ## pwf <- attr(my.res, 'pwf')

  ## run goseq::goseq
  g2c <- transform(as.data.frame(gsd),
                   category=encode_gskey(collection, name),
                   stringsAsFactors=FALSE)
  g2c <- g2c[, c('category', 'featureId')]
  pwf <- suppressWarnings(goseq::nullp(degenes, bias.data=mylens))
  goseq.res <- suppressWarnings({
    goseq::goseq(pwf, gene2cat=g2c, method='Wallenius',
                 use_genes_without_cat=TRUE)
  })

  ## Match up and compare
  expect_true(setequal(my.res$category, goseq.res$category))
  goseq.res <- goseq.res[match(my.res$category, goseq.res$category),]
  expect_equal(goseq.res, my.res, check.attributes=FALSE)

  ## Run goseq through multiGSEA to make sure it matches goseq.res
  mg <- multiGSEA(gsd, vm, vm$design, methods="goseq", feature.bias=mylens)
  expect_true(setequal(resultNames(mg), c("goseq", "goseq.up", "goseq.down")))
  my2 <- result(mg, 'goseq')
  my2$key <- encode_gskey(my2)
  expect_equal(my2$key, goseq.res$category)
  expect_equal(my2$pval, goseq.res$over_represented_pvalue)
  expect_equal(my2$pval.under, goseq.res$under_represented_pvalue)
  expect_equal(my2$n, goseq.res$numInCat)
  expect_equal(my2$n.sig, goseq.res$numDEInCat)
})



test_that("goseq hypergeometric pvals close to do.hyperGeometricTest", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsd <- conform(gsd, vm)
  mylens <- setNames(vm$genes$size, rownames(vm))

  mg <- multiGSEA(gsd, vm, vm$design)
  lfc <- logFC(mg)
  selected <- subset(lfc, significant)$featureId
  universe <- rownames(vm)

  mygoh <- suppressWarnings({
    multiGSEA::goseq(gsd, selected, universe, method='Hypergeometric',
                     feature.bias=mylens)
  })

  myhyp <- multiGSEA::hyperGeometricTest(gsd, selected, universe)
  if (FALSE) {
    plot(-log10(mygoh$over_represented_pvalue), -log10(myhyp$pval))
    abline(0,1,col='red')
  }

  expect_equal(mygoh$over_represented_pvalue, myhyp$pval, tolerance=0.025)
})
