context("goseq")

test_that("multiGSEA(method='goseq') requires valid feature.bias vector", {
  vm <- exampleExpressionSet(do.voom=TRUE)
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
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsd <- conform(gsd, vm)

  mg <- multiGSEA(gsd, vm, vm$design)
  lfc <- logFC(mg)
  selected <- subset(lfc, significant)$featureId
  universe <- rownames(vm)
  mylens <- setNames(vm$genes$size, rownames(vm))
  degenes <- setNames(integer(length(universe)), universe)
  degenes[selected] <- 1L

  my.res <- suppressWarnings({
    multiGSEA::goseq(gsd, selected, universe, mylens, method='Wallenius',
                     use_genes_without_cat=TRUE, .pipelined=TRUE)
  })
  ## pwf <- attr(my.res, 'pwf')

  ## expected
  g2c <- transform(as.data.frame(gsd), category=name)
  g2c <- g2c[, c('category', 'featureId')]
  pwf <- suppressWarnings(goseq::nullp(degenes, bias.data=mylens))
  goseq.res <- suppressWarnings({
    goseq::goseq(pwf, gene2cat=g2c, method='Wallenius',
                 use_genes_without_cat=TRUE)
  })

  ## Match up and compare
  goseq.res <- goseq.res[match(my.res$name, goseq.res$category),]
  my <- data.table::setnames(data.table::copy(my.res),
                             c("name", "n", "n.drawn"),
                             c("category", "numInCat", "numDEInCat"))
  my <- my[, names(goseq.res)]
  expect_equal(my, goseq.res, check.attributes=FALSE)

  mgs <- multiGSEA(gsd, vm, vm$design, methods='goseq', split.updown=FALSE,
                   feature.bias=mylens)

  r <- result(mgs, 'goseq')
  data.table::setnames(r,
                       c('n.sig', 'pval', 'pval.under'),
                       c('n.drawn', 'over_represented_pvalue',
                         'under_represented_pvalue'))
  r <- r[, names(my.res)]
  expect_equal(r, my.res, info="multiGSEA(method='goseq')", check.attributes=FALSE)
})

test_that("goseq,split.updown=TRUE is sane", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsd <- conform(gsd, vm)

  mylens <- setNames(vm$genes$size, rownames(vm))

  mgs <- multiGSEA(gsd, vm, vm$design, methods='goseq', split.updown=TRUE,
                   feature.bias=mylens)
  mg <- multiGSEA(gsd, vm, vm$design, methods='goseq', split.updown=FALSE,
                  feature.bias=mylens)

  expect_true(setequal(resultNames(mgs), c('goseq', 'goseq.up', 'goseq.down')))
  expect_equal(result(mgs, 'goseq'), result(mg, 'goseq'))
})

test_that("goseq,split.updown=TRUE plays well with other multiGSEA methods", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsd <- conform(gsd, vm)
  mylens <- setNames(vm$genes$size, rownames(vm))

  mgs <- multiGSEA(gsd, vm, vm$design, methods='goseq', split.updown=TRUE,
                   feature.bias=mylens)
  mga <- multiGSEA(gsd, vm, vm$design, methods=c('goseq', 'camera'),
                   split.updown=TRUE, feature.bias=mylens)
  expect_true(setequal(resultNames(mga),
                       c('camera', 'goseq', 'goseq.up', 'goseq.down')))

  both <- intersect(resultNames(mgs), resultNames(mga))
  for (wut in both) {
    expect_equal(result(mgs, wut), result(mga, wut), info=wut)
  }
})


test_that("goseq hypergeometric test is like do.hyperGeometricTest", {
  vm <- exampleExpressionSet(do.voom=TRUE)
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

  expect_equal(mygoh$pval_over, myhyp$pval, tolerance=0.015)
})
