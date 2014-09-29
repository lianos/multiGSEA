## The datasets references in this file are generated in:
##   tests/testdata/setup-testdata.R

##' @export
##' @import Biobase
exampleExpressionSet <- function(dataset=c('tumor-vs-normal', 'tumor-subtype'),
                                 do.voom=TRUE, voom.plot=FALSE) {
  dataset <- match.arg(dataset)
  es.all <- readRDS(system.file('extdata', 'testdata', 'TCGA-BRCA-some.es.rds',
                                package='multiGSEA'))

  if (dataset == 'tumor-vs-normal') {
    es <- es.all
    design <- model.matrix(~ Cancer_Status, pData(es))
    colnames(design) <- sub('Cancer_Status', '', colnames(design))
  } else {
    es <- es.all[, es.all$Cancer_Status == 'tumor']
    pData(es) <- droplevels(pData(es))
    design <- model.matrix(~ 0 + PAM50subtype, pData(es.sub))
    colnames(design) <- sub('PAM50subtype', '', colnames(design))
  }

  if (do.voom) voom(es, design, plot=voom.plot) else es
}

##' Get sample gene sets
##'
##' @export
##'
##' @param as Character vector to specify the type of returned geneset
##' list that you want.
##'
##' @return A list of lists of entrezIDs when \code{as == 'lol'}, or
##' a list of integers into the rows of \code{exampleExpressionSet}
##' for the genes in the given geneset.
exampleGeneSets <- function(as=c('lol', 'limma')) {
  as <- match.arg(as)
  fn <- switch(as,
               lol='genesets-multiGSEA-list-of-lists.rds',
               limma='genesets-limma-idxvectors.rds')
  readRDS(system.file('extdata', 'testdata', fn, package='multiGSEA'))
}
