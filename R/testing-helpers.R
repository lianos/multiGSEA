## The datasets references in this file are generated in:
##   tests/testdata/setup-testdata.R

##' Functions that load data for use in examples and testing.
##'
##' @description
##' We provide examplar expression data (counts or voomed) as well as exemplar
##' gene sets in different forms.
##'
##' @section exampleExpressionSet:
##' The expression data is a subset of the TCGA BRCA indication. Calling
##' \code{exampleExpressionSet} will load either \code{limma} or \code{Biobase}
##' into the calling environment so that the user can manipulate the object
##' that is returned.
##'
##' @rdname examples
##' @aliases exampleExpressionSet
##'
##' @export
##' @import Biobase
##'
##' @param dataset Character vector indicating what samples wanted, either
##'   \code{"tumor-vs-normal"} for a tumor vs normal dataset from TCGA, or
##'   just the tumor samples from the same annotated with subtype.
##' @param do.voom If TRUE, a voomed EList is returned, otherwise an
##'   ExpressionSet of counts.
exampleExpressionSet <- function(dataset=c('tumor-vs-normal', 'tumor-subtype'),
                                 do.voom=TRUE) {
  suppressPackageStartupMessages({
    ## voom is having some issues finding fData if I don't do this, and I don't
    ## want to waste time debuggin this furthe
    require('Biobase', character.only=TRUE)
  })
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

  out <- es
  attr(out, 'design') <- design
  if (do.voom) {
    ## require('limma', character.only=TRUE)
    out <- voom(es, design, plot=FALSE)
    out$genes <- fData(es)
  } else {
  }

  out
}


##' @section exampleGeneSets:
##' Returns gene sets as either a list of feature identifiers or integers
##' that index into a target expression object \code{x}.
##'
##' @rdname examples
##' @aliases exampleGeneSets
##'
##' @export
##'
##' @param x If provided, an expression/matrix object so that the genesets are
##'   returned as (integer) index vectors into the rows of x whose rownames
##'   match the ids in the geneset.
##' @return A list of lists of entrezIDs when \code{as == 'lol'}, or
##'   a list of integers into the rows of \code{exampleExpressionSet}
##'   for the genes in the given geneset.
exampleGeneSets <- function(x, unlist=!missing(x)) {
  gsl.fn <- system.file('extdata', 'testdata',
                        'genesets-multiGSEA-list-of-lists.rds',
                        package='multiGSEA')
  gsl <- readRDS(gsl.fn)
  if (!missing(x)) {
    gsl <- lapply(gsl, function(col) {
      lapply(col, function(gs) {
        out <- match(gs, rownames(x))
        out[!is.na(out)]
      })
    })
  }
  if (unlist) {
    gsl <- unlist(gsl, recursive=FALSE)
    names(gsl) <- sub('\\.', ';;', names(gsl))
  }
  gsl
}

##' @section exampleGeneSetDb:
##' Returns gene sets as a \code{GeneSetDb} object
##'
##' @rdname examples
##' @aliases exampleGeneSetDb
##' @export
exampleGeneSetDb <- function() {
  out <- GeneSetDb(exampleGeneSets())
  colls <- unique(collectionMetadata((out))$collection)
  fn <- function(x, y) {
    paste0('http://www.broadinstitute.org/gsea/msigdb/cards/', y, '.html')
  }
  for (col in colls) {
    geneSetCollectionURLfunction(out, col) <- fn
  }
  out
}

##' @section exampleGeneSetDF:
##' Returns a data.frame of gene set definitions. A data.frame of this form
##' can be passed into the \code{GeneSetDb} contructor.
##'
##' @export
##' @importFrom utils read.csv
##' @rdname examples
##' @aliases exampleGeneSetDF
exampleGeneSetDF <- function() {
  gs.df <- system.file('extdata', 'testdata', 'custom-sigs.csv',
                       package='multiGSEA')
  gs.df <- read.csv(gs.df, stringsAsFactors=FALSE)
  gs.df$featureId <- as.character(gs.df$featureId)
  gs.df
}


##' @export
##' @rdname examples
##' @aliases exampleGeneSetDF
exampleMultiGSEAResult <- function() {
  fn <- system.file('extdata', 'testdata', 'test-MultiGSEAResult.rds',
                    package='multiGSEA')
  readRDS(fn)
}
