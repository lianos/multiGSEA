## This code really shouldn't be in this package
##
## Utility functions for working with expression container objects, ie. with
## the following classes:
##  * ExpressionSet
##  * EList
##  * DGEList
##  * SummarizedExperiment

## Expression Values -----------------------------------------------------------

##' Provides a "regularized" version of edgeR::cpm
##'
##' A modified version of edgeR::cpm which uses interprets \code{prior.count}
##' to parameterize a Poisson distribution to sample read counts from for
##' each observation
##'
##' Note that when called with \code{regularized == FALSE}, this
##' function should behave almost identically to \code{\link[edgeR]{cpm}}.
##'
##' @export
##'
##' @param x A \code{DGEList} or something that can be converted to a count
##'   matrix via \code{as.matrix(x)}
##' @param normalized.lib.sizes \code{logical}: use normalized lib sizes?
##' @param log \code{log} return expression in log2? (defaults to \code{TRUE})
##' @param prior.count The lambda/mean of the Poisson distribution used to add
##'   pseudo counts to if \code{regularized == TRUE}. Otherwise a constant
##'   count added to each observation. Note that a A minimum of
##'   \code{prior.smooth} is always added to the poisson noise to avoid taking
##'   logs of 0
##' @param prior.smooth A small number that is added to the Poisson random var
##'   drawn by \code{prior.count} in order to avoid taking logs of 0. This will
##'   always be at least 0.25
##' @param regularized If TRUE, prior.count parameterize a Poission distribution
##'   to draw from, otherwise prior.count is constant. When \code{FALSE}, this
##'   function should behave almost identically to \code{\link[edgeR]{cpm}}
##' @param norm.factors precalculated norm.factors (numeric vector as long as
##'   \code{ncol(x)})
##' @param lib.size precalculate lib.size vector (numeric vector as long as
##'   \code{ncol(x)})
##' @param norm.factors.col A column to fetch from \code{.pData(x)} that has
##'   precalculated norm.factors
##' @param lib.size.col A column to fetch from \code{.pData(x)} that has
##'   precalculated lib.sizes
##' @param x.element Additional param to help get counts out of \code{x}: sent
##'   to \code{.exprs,element} param
##' @param ... not sure why I have this here
##' @return A \code{matrix,numeric} of the normalized/smoothed expression values
##'   from \code{x}.
CPM <- function(x, normalized.lib.sizes=TRUE, log=TRUE,
                prior.count=5, prior.smooth=0.25,
                regularized=prior.count != prior.smooth,
                norm.factors=NULL, lib.size=NULL,
                norm.factors.col='norm.factors',
                lib.size.col='lib.size',
                x.element='counts', ...) {
  e <- .exprs(x, x.element)
  if (is.null(lib.size) && !is.matrix(x)) {
    lib.size <- .pData(x)[[lib.size.col]]
  }
  if (is.null(lib.size)) {
    lib.size <- colSums(e)
  }
  if (is.null(norm.factors) && !is.matrix(x)) {
    norm.factors <- .pData(x)[[norm.factors.col]]
  }
  if (is.null(norm.factors)) {
    norm.factors <- rep(1, ncol(x))
  }
  if (normalized.lib.sizes) {
    lib.size <- lib.size * norm.factors
  }

  ## prior.count.scaled <- lib.size/mean(lib.size) * prior.count
  ## lib.size <- lib.size + 2 * prior.count.scaled
  prior.smooth <- min(0.25, prior.smooth)
  prior.count <- max(0L, if (regularized) ceiling(prior.count) else prior.count)
  norm.lib.size <- lib.size / mean(lib.size)
  prior.count.scaled <- norm.lib.size * prior.count
  prior.smooth.scaled <- norm.lib.size * prior.smooth

  e <- sapply(1:ncol(e), function(i) {
    if (regularized) {
      pcount <- rpois(nrow(e), prior.count) * norm.lib.size[i]
    } else {
      pcount <- prior.count.scaled[i]
    }
    out <- e[, i] + pcount
    replace(out, out == 0, prior.smooth.scaled[i])
  })
  colnames(e) <- colnames(x)

  ## lib.size <- lib.size + 2 * prior.count.min.scaled
  lib.size.offset <- if (prior.count == 0) {
    prior.smooth.scaled
  } else {
    prior.count.scaled
  }
  lib.size <- lib.size + 2 * lib.size.offset
  lib.size <- 1e-06 * lib.size

  out <- t(t(e) / lib.size)
  if (log) {
    out <- log2(out)
  }

  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)
  attr(out, 'lib.size') <- lib.size
  out
}

##' Calculates between-sample normalized transcripts per million
##'
##' There is some consensus that
##' \href{https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/}{TPM} is a good measure of transcript/gene abundance
##'
##' @export
##' @seealso \code{\link{CPM2TPM}}, \code{\link{CPM}}, \code{\link{RPKM}}
##' @inheritParams CPM
##' @inheritParams create.glength.vector
##' @return matrix of transcripts per million
TPM <- function(x, gene.info=NULL, normalized.lib.sizes=TRUE, log=TRUE,
                prior.count=5, prior.smooth=0.25,
                regularized=prior.count != prior.smooth,
                norm.factors=NULL, lib.size=NULL,
                norm.factors.col='norm.factors',
                lib.size.col='lib.size',
                x.element='counts', ...) {
  xx <- CPM(x, normalized.lib.sizes=normalized.lib.sizes,
            log=FALSE, prior.count=prior.count, prior.smooth=prior.smooth,
            regularized=regularized,
            norm.factors=norm.factors, lib.size=lib.size,
            norm.factors.col=norm.factors.col,
            lib.size.col=lib.size.col, x.element=x.element, ...)
  out <- CPM2TPM(xx, gene.info, ...)
  if (log) log2(out) else out
}

##' Convert CPM matrix to TPM
##'
##' @export
##'
##' @param x An expression matrix of counts per million (on the natural scale)
##' @param gene.info A data.frame with entrez to gene length?
##' @param id.col The colum in \code{gene.info} that has the ids of the genes,
##'   which is used to match against \code{rownames(x)}.
##' @param size.col The column in \code{gene.info} that has the gene lengths
##' @param impute.missing.lengths If \code{TRUE}, when the length for a gene
##'   is not found, we use the mean of all gene lengths.
##' @param ... Moar things
##' @inheritParams create.glength.vector
##' @return A matrix of TPMs
CPM2TPM <- function(x, gene.info=NULL, id.col='entrez.id', size.col='size',
                    impute.missing.lengths=TRUE, ...) {
  xlens <- create.glength.vector(x, gene.info, id.col, size.col,
                                 impute.missing.lengths)
  apply(x, 2, function(y) {
    yy <- y * 1e6 ## Bring it back to raw counts
    rate <- log(yy) - log(xlens)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  })
}

##' Modfied version of edgeR RPKM for easy access to gene size
##'
##' @importFrom stringr str_extract
##' @export
##'
##' @param x expression object
##' @param gene.info data.table that has the lengths for the genes in x.
##'   If \code{NULL} this will fetch a default gene size table that I extracted
##'   from the TCGA::COAD datasets in ExpressionPlot
##' @param regularized passed down to CPM, but set to \code{FALSE} by default
##'   here.
##' @param impute.missing.lengths If \code{TRUE}, when the length for a gene
##'   is not found, we use the mean of all gene lengths.
##' @param id.col The name of the column in \code{gene.info} that has the IDs
##'   that match to \code{rownames(x)}
##' @param size.col The name of the column in \code{gene.info} that has the
##'   gene lengths to normalize by
##' @param ... Moar args
##' @return a matrix of RPKM values
RPKM <- function(x, gene.info=NULL, regularized=FALSE,
                 impute.missing.lengths=TRUE,
                 id.col='entrez.id', size.col='size',
                 ...) {

  xlens <- create.glength.vector(x, gene.info, id.col, size.col,
                                 impute.missing.lengths)
  ## ---------------------------------------------------------------------------
  ## calculate rpkms
  ##
  ## Note I original divided by the length by making a (repetative) matrix of
  ## the length and subtracted that to make sure we did element-by-element
  ## operations. I then tested that output against this where we just use a
  ## single vector on the RHS, and they are equivalent.
  out <- CPM(x, regularized=regularized, ...) - log2(xlens / 1000)
  attr(out, 'K') <- xlens
  out
}

create.glength.vector <- function(x, gene.info=NULL,
                                  id.col='entrez.id', size.col='size',
                                  impute.missing.lengths=TRUE) {
  if (is.null(gene.info)) {
    fn <- system.file('extdata', 'gene-info', 'Homo_sapiens.csv',
                      package='multiGSEA')
    gene.info <- as.data.frame(fread(fn))
  }

  ## ---------------------------------------------------------------------------
  ## Process gene.info and match lenghts to rows of x
  stopifnot(is.data.frame(gene.info))
  stopifnot(id.col %in% names(gene.info))
  stopifnot(size.col %in% names(gene.info))
  ## Match lengths to genes/rows in x
  lens <- setNames(gene.info[[size.col]], gene.info[[id.col]])
  ## Work around GeneID: prefix bit
  x.is.gid <- substr(sample(rownames(x), 1), 1, 7) == 'GeneID:'
  if (x.is.gid) {
    names(lens) <- paste0('GeneID:', str_extract(names(lens), '(\\d)+'))
  }

  xref <- match(rownames(x), names(lens))
  xlens <- lens[xref]
  no.match <- is.na(xlens)
  if (any(no.match)) {
    if (!impute.missing.lengths) {
      stop("Can't match rownames(x) to ids in gene.info")
    }
    if (mean(no.match) > 0.30) {
      stop("Will not impute missing lengths with < 50% genes matched")
    }
    warning(sprintf("Imputing %.2f%% of gene lengths", mean(no.match)),
            immediate.=TRUE)
    xlens[no.match] <- median(xlens[!no.match])
  }

  xlens
}

## Accessor Functions ----------------------------------------------------------

##' Fetch the sample information \code{data.frame}
.pData <- function(x) {
  if (is(x, 'eSet')) {
    out <- pData(x)
  } else if (is(x, 'DGEList')) {
    out <- x$samples
  } else if (is(x, 'EList')) {
    out <- x$targets
  } else if (is(x, 'SummarizedExperiment')) {
    require('GenomicRanges', quietly=TRUE)
    out <- as.data.frame(colData(x))
  } else {
    stop("Don't know how to get pData for class:", class(x))
  }
  rownames(out) <- colnames(x)
  out
}

##' Fetch the row/gene information \code{data.frame}
.fData <- function(x) {
  if (is(x, 'eSet')) {
    out <- fData(x)
  } else if (is(x, 'DGEList') || is(x, 'EList')) {
    out <- x$genes
  } else if (is(x, 'SummarizedExperiment')) {
    if (is(rowData(x, 'GRanges'))) {
      out <- as.data.frame(mcols(x))
    } else {
      stopifnot(is(x, 'GRangesList'))
      out <- as.data.frame(mcols(rowData(x)))
    }
  } else {
    stop("Don't know how to get fData for class: ", class(x))
  }
  rownames(out) <- rownames(x)
  out
}

##' Fetch expression values out of expression container \code{x}
.exprs <- function(x, element='counts') {
  if (is(x, 'eSet')) {
    enames <- assayDataElementNames(x)
    if (length(enames) == 1L) {
      if (element != enames) {
        warning("Using ", enames, " for assayDataElement", immediate.=TRUE)
        element <- enames
      } else {
        stopifnot(element %in% enames)
      }
    }
    out <- assayDataElement(x, element)
  } else if (is(x, 'DGEList')) {
    out <- x$counts
  } else if (is(x, 'EList')) {
    out <- x$E
  } else if (is(x, 'SummarizedExperiment')) {
    enames <- assayNames(x)
    if (length(enames) == 1L) {
      if (element != enames) {
        warning("Using ", enames, " for assayDataElement", immediate.=TRUE)
        element <- enames
      } else {
        stopifnot(element %in% enames)
      }
    }
    out <- assay(x, element)
  } else if (is.matrix(x)) {
    out <- x
  } else {
    stop("Don't know how to get exprs from class:", class(x))
  }
  if (!is.matrix(out)) {
    stop("Failed to get expression out of object of class ", class(x),
         " using element: ", element)
  }
  rownames(out) <- rownames(x)
  colnames(out) <- colnames(x)
  out
}

## Conversion Functions --------------------------------------------------------

##' Utility function to convert ExpressionSet to DGEList
##'
##' @export as.DGEList
##' importFrom edgeR DGEList calcNormFactors
##'
##' @param x ExpressionSet to convert
##' @param element The name of the \code{assayDataElement(x)} to use for the
##'   count matrix from \code{x}
##' @param lib.size The total lib.size to use, will take \code{colSums} of
##'   count matrix if not defined: this is a \code{numeric} vector as long
##'   as the numbe of samples
##' @param norm.factors A numeric vector of norm.factors to use. If this is
##'   NULL, then they will be recalculated
##' @param group Either a string that is a column of \code{.pData(x)} to use
##'   as grouping information, or a factor for grouping information. If not
##'   used, the group for all samples will be set to \code{factor(1)}.
##' @param ... moar args
##' @return A DGEList with norm factors
as.DGEList <- function(x, element='counts', lib.size=NULL,
                       norm.factors=.pData(x)[['norm.factors']],
                       group=factor(rep(1L, ncol(x))), ...) {
  pd <- .pData(x)
  fd <- .fData(x)
  if (is.character(group) && length(group) == 1) {
    group <- pd[[group]]
  }
  if (!is.atomic(group)) {
    stop("group needs to be atomic")
  }
  if (length(group) != ncol(x)) {
    stop("Illegal group param")
  }
  y <- DGEList(.exprs(x, element), lib.size=lib.size, genes=.fData(x),
               group=group)
  if (is.null(norm.factors)) {
    y <- calcNormFactors(y, ...)
  } else {
    stopifnot(is.numeric(norm.factors) && length(norm.factors) == ncol(x))
    y$samples$norm.factors <- norm.factors
  }
  add.cols <- setdiff(names(pd), names(y$samples))
  y$samples <- cbind(y$samples, pd[, add.cols, drop=FALSE])
  y
}

##' Combines two (already concordant) DGEList objects
##'
##' @export combine.DGEList
##' @usage combine.DGEList(x, y)
##' @param x \code{DGEList}
##' @param y \code{DGEList}
##' @return combined \code{DGEList}
combine.DGEList <- function(x, y) {
  stopifnot(is(x, 'DGEList') && is(y, 'DGEList'))
  stopifnot(all.equal(rownames(x), rownames(y)))
  stopifnot(length(intersect(colnames(x), colnames(y))) == 0)
  x$counts <- cbind(x$counts, y$counts)
  x$samples <- rbind(x$samples, y$samples)
  x
}

## Melting ---------------------------------------------------------------------

##' Melts an Expression Container
##'
##' @export
##'
##' @param x A \code{ExpressionSet}, \code{DGEList}, etc.
##' @param pd The "phenodata" columns to keep in melted object. Either the
##'   \code{data.frame}, or the column names of the \code{data.frame} that
##'   will be extracted from the expression object.
##' @param fd The "feature data" columns to keep in melted object. Either the
##'   \code{data.frame}, or the column names of the \code{data.frame} that
##'   will be extracted from the expression object.
##' @param element The name of "the assay" that holds the expression data
##' @param exform A function to call on x to get the (possibly transformed)
##'   data
meltx <- function(x, pd=NULL, fd=NULL, element='counts', exform=NULL,
                  value.name='expr') {
  if (length(pd) == 0) pd <- NULL
  if (length(fd) == 0) fd <- NULL
  if (is.character(pd)) {
    if (is.matrix(x)) {
      stop("character vector not allowed for `pd` when `x` is a matrix")
    }
    .pd <- .pData(x)
    bad.pd <- setdiff(pd, names(.pd))
    if (length(bad.pd)) {
      stop("Illegal pheno data column names: ", paste(bad.pd, collapse=','))
    }
    pd <- .pd[, pd, drop=FALSE]
  }
  if (is.character(fd)) {
    if (is.matrix(x)) {
      stop("character vector not allowed for `fd` when `x` is a matrix")
    }
    .fd <- .fData(x)
    bad.fd <- setdiff(fd, names(.fd))
    if (length(bad.fd)) {
      stop("Illegal feature data column names: ", paste(bad.fd, collapse=','))
    }
    fd <- .fd[, fd, drop=FALSE]
  }
  stopifnot(is.null(pd) || (is.data.frame(pd) && nrow(pd) == ncol(x)))
  stopifnot(is.null(fd) || (is.data.frame(fd) && nrow(fd) == nrow(x)))
  if (is.function(exform)) {
    x <- exform(x)
  }
  if (is(x, 'eSet')) {
    meltx.eSet(x, pd, fd, element, value.name)
  } else if (is(x, 'DGEList')) {
    meltx.DGEList(x, pd, fd, value.name)
  } else if (is(x, 'matrix')) {
    meltx.matrix(x, pd, fd, value.name)
  } else {
    stop("Unknown expression object type for x: ", class(x))
  }
}

meltx.eSet <- function(x, pd, fd, element='counts', exform=NULL,
                       value.name='expr') {
  M <- assayDataElement(x, element)
  meltx.matrix(M, pd, fd, value.name)
}

meltx.DGEList <- function(x, pd, fd, value.name='expr') {
  meltx.matrix(x$counts, pd, fd, value.name)
}

meltx.matrix <- function(x, pd, fd, value.name='expr') {
  if (is.data.frame(pd)) {
    stopifnot(colnames(x) == rownames(pd))
  } else {
    stopifnot(is.null(pd))
  }
  if (is.data.frame(fd)) {
    stopifnot(rownames(x) == rownames(fd))
  } else {
    stopifnot(is.null(fd))
  }

  m <- melt(x)
  setnames(m, c('geneID', 'sample', value.name))
  m <- transform(m, geneID=as.character(geneID), sample=as.character(sample))

  if (!is.null(fd)) {
    fd.xref <- match(m$geneID, rownames(fd))
    m <- cbind(m, fd[fd.xref,,drop=FALSE])
  }
  if (!is.null(pd)) {
    pd.xref <- match(m$sample, rownames(pd))
    m <- cbind(m, pd[pd.xref,,drop=FALSE])
  }
  rownames(m) <- NULL
  m
}
