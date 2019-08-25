# The datasets references in this file are generated in:
#   tests/testdata/setup-testdata.R

#' Functions that load data for use in examples and testing.
#'
#' We provide examplar expression data (counts or voomed) as well as exemplar
#' gene sets in different forms.
#'
#' @section exampleExpressionSet:
#' The expression data is a subset of the TCGA BRCA indication. Calling
#' `exampleExpressionSet(do.voom = TRUE)` will return a voomed `EList` version
#' of the data. When `do.voom = FALSE`, you will get a DGEList of the counts
#'
#' @rdname examples
#' @aliases exampleExpressionSet
#'
#' @export
#' @importFrom limma voom
#'
#' @param dataset Character vector indicating what samples wanted, either
#'   \code{"tumor-vs-normal"} for a tumor vs normal dataset from TCGA, or
#'   just the tumor samples from the same annotated with subtype.
#' @param do.voom If TRUE, a voomed EList is returned, otherwise an
#'   ExpressionSet of counts.
exampleExpressionSet <- function(dataset=c('tumor-vs-normal', 'tumor-subtype'),
                                 do.voom=TRUE) {
  dataset <- match.arg(dataset)
  fn <- system.file("extdata", "testdata", "TCGA-BRCA-some.DGEList.rds",
                    package = "multiGSEA", mustWork = TRUE)
  y.all <- readRDS(fn)

  # Two samples seem to be outliers:
  axe.samples <- c("TCGA-A2-A3XV-01A-21R-A239-07", "TCGA-A2-A3XU-01A-12R-A22U-07")
  y.all <- y.all[, !colnames(y.all) %in% axe.samples]

  if (dataset == 'tumor-vs-normal') {
    y <- y.all
    design <- model.matrix(~ Cancer_Status, y$samples)
    colnames(design) <- sub('Cancer_Status', '', colnames(design))
  } else {
    y <- y.all[, y.all$samples$Cancer_Status == 'tumor']
    y$samples <- droplevels(y$samples)
    design <- model.matrix(~ 0 + PAM50subtype, y$samples)
    colnames(design) <- sub('PAM50subtype', '', colnames(design))
  }

  y$design <- design

  if (do.voom) voom(y, y$design, plot = FALSE) else y
}


#' @section exampleGeneSets:
#' Returns gene sets as either a list of feature identifiers or integers
#' that index into a target expression object `x`.
#'
#' @rdname examples
#' @aliases exampleGeneSets
#'
#' @export
#'
#' @param x If provided, an expression/matrix object so that the genesets are
#'   returned as (integer) index vectors into the rows of x whose rownames
#'   match the ids in the geneset.
#' @return A list of lists of entrezIDs when \code{as == 'lol'}, or
#'   a list of integers into the rows of \code{exampleExpressionSet}
#'   for the genes in the given geneset.
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

#' @section exampleGeneSetDb:
#' Returns gene sets as a `GeneSetDb` object
#'
#' @rdname examples
#' @aliases exampleGeneSetDb
#' @export
exampleGeneSetDb <- function() {
  out <- GeneSetDb(exampleGeneSets())
  colls <- unique(collectionMetadata(out, as.dt = TRUE)$collection)
  fn <- function(x, y) {
    paste0('http://www.broadinstitute.org/gsea/msigdb/cards/', y, '.html')
  }
  # browser()
  for (col in colls) {
    geneSetCollectionURLfunction(out, col) <- fn
  }
  out
}

#' @section exampleGeneSetDF:
#' Returns a data.frame of gene set definitions. A data.frame of this form
#' can be passed into the `GeneSetDb()` contructor.
#'
#' @export
#' @importFrom utils read.csv
#' @rdname examples
#' @aliases exampleGeneSetDF
exampleGeneSetDF <- function() {
  gs.df <- system.file('extdata', 'testdata', 'custom-sigs.csv',
                       package='multiGSEA')
  read.csv(gs.df, stringsAsFactors=FALSE, colClasses='character')
}


#' @export
#' @rdname examples
#' @aliases exampleGeneSetDF
exampleMultiGSEAResult <- function(cached=TRUE) {
  if (cached) {
    fn <- system.file('extdata', 'testdata', 'test-MultiGSEAResult.rds',
                      package='multiGSEA')
    out <- readRDS(fn)
  } else {
    vm <- exampleExpressionSet()
    gdb <- exampleGeneSetDb()
    out <- multiGSEA(gdb, vm, vm$design, "tumor", c('camera', 'fry'))
  }
  out
}

#' `exampleDgeResult` returns a data.frame of differential expression results.
#' Currently they are for human ensembl genes. Setting the `induce.bias`
#' parameter to `"effective_length"` or `"AveExpr"` will munge the returned
#' result such that larger "bias" values will be associated to lower pvalues,
#' so we can more easily test biased enrichment approaches like [enrichtest()]
#' and [goseq()].
#'
#' @export
#' @rdname examples
exampleDgeResult <- function(species = "human", id.type = "ensembl",
                             induce.bias = NULL) {
  # we only have human/ensembl for now
  species <- match.arg(species, "human")
  id.type <- match.arg(id.type, "ensembl")
  dge.fn <- system.file("extdata", "testdata", "dataframe-input.csv",
                        package = "multiGSEA")
  out <- read.csv(dge.fn, stringsAsFactors = FALSE)
  if (is.character(induce.bias)) {
    bias <- match.arg(induce.bias, c("effective_length", "AveExpr"))
    o <- order(out[["pval"]])
    out[[bias]][o] <- sort(out[[bias]], decreasing = TRUE)
  }
  out
}
