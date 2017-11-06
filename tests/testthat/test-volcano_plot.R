# context("Interactive Volcano")
#
# expect_that("volcano plots with hexbins", {
#   x <- exampleMultiGSEAResult()
#
#   stats='dge'; xaxis='logFC'; yaxis='pval'
#   xtfrm <- identity
#   ytfrm <- function(vals) -log10(vals)
#   xhex <- 1
#   yhex <- 0.05
#
#   dat <- volcanoStatsTable(x, stats, xaxis, yaxis, idx, xtfrm, ytfrm)
#
#
# })
