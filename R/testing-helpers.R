##' Fetch a pre-canned expression dataset that can be used for testing and
##' examples.
##'
##' The samples in this dataset are a subset of samples from TCGA-BRCA
##'
##' @export
##' @import Biobase
##'
##' @param dataset Character vector indicating what samples wanted, either
##' \code{"tumor-vs-normal"} for a tumor vs normal dataset from TCGA, or just the tumor samples
##' from the same (there are different subtypes in there)
##' @param do.voom If TRUE, a voomed EList is returned, otherwise an
##' ExpressionSet of counts.
exampleExpressionSet <- function(dataset=c('tumor-vs-normal', 'tumor-subtype'),
                                 do.voom=TRUE) {
  dataset <- match.arg(dataset)
  es.all <- readRDS(system.file('extdata', 'testdata', 'TCGA-BRCA-some.es.rds',
                                package='multiGSEA'))
  pData(es.all) <- within(pData(es.all), {
    Cancer_Status <- factor(as.character(Cancer_Status))
    PAM50subtype <- factor(as.character(PAM50subtype))
  })

  if (dataset == 'tumor-vs-normal') {
    es <- es.all
    design <- model.matrix(~ Cancer_Status, pData(es))
    colnames(design) <- sub('Cancer_Status', '', colnames(design))
  } else {
    es <- es.all[, es.all$Cancer_Status == 'tumor']
    pData(es) <- droplevels(pData(es))
    design <- model.matrix(~ 0 + PAM50subtype, pData(es))
    colnames(design) <- sub('PAM50subtype', '', colnames(design))
  }

  if (do.voom) voom(es, design, plot=FALSE) else es
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

## The datasets references in this file are generated in:
##   tests/testdata/setup-testdata.R
## Let's document them here

