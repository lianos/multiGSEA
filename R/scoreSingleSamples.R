##' Generates single sample gene set scores across a datasets by many methods
##'
##' It is common to assess the activity of a gene set in a given sample. There
##' are many ways to do that, and this method is analogous to the
##' \code{\link{multiGSEA}} function in that it enables the user to run a
##' multitude of single-sample-gene-set-scoring algorithms over a target
##' expression matrix using a \code{\link{GeneSetDb}} object.
##'
##' Please refer to the "Generating Single Sample Gene Set Scores" of the
##' multiGSEA vignette for further exposition.
##'
##' The following \code{methods} are currenly provided.
##'
##' \describe{
##'   \item{ewm}{
##'     The \code{\link{eigenWeightedMean}} calculates the fraction each gene contributes
##'     to a pre-specified principal component. These contributions acts as
##'     weights over each gene, which are then used in a simple weighted mean
##'     calculation over all the genes in the geneset per sample. This is
##'     similar, in spirit, to the svd/gsdecon method ("gsd")}
##'   \item{gsd}{
##'     This method was first introduced by Jason Hackney in
##'     \href{https://doi.org/10.1038/ng.3520}{doi:10.1038/ng.3520}. Please
##'     refer to the \code{\link{gsdScore}} function for more information.}
##'   \item{ssgsea}{Using ssGSEA as implemented in the GSVA package}
##'   \item{zscore}{
##'     The features in the expression matrix are rowwise z transformed. The
##'     gene set level score is then calculated by adding up the zscores for
##'     the genes in the gene set, then dividing that number by either the
##'     the size (or its sqaure root (default)) of the gene set.}
##'   \item{mean}{
##'     Simply take the mean of the values from the expression matrix that are
##'     in teh gene set. Right or wrong, sometimes you just want the mean
##'     without transforming the data.
##'   }
##'   \item{gsva}{The gsva method of GSVA package}
##'   \item{plage}{Using "plage" as implemented in the GSVA package}
##' }
##'
##' @export
##' @importFrom matrixStats rowSds
##'
##' @param gdb A GeneSetDb
##' @param y An expression matrix to score genesets against
##' @param methods A character vector of methods to score samples by
##' @param as.matrix Return results as a list of matrices instead of a melted
##'   data.frame? Defaults to \code{FALSE}.
##' @param drop.sd Genes with a standard deviation across columns in \code{y}
##'   with a standard deviation less than this value will be dropped.
##' @param verbose make some noise? Defaults to \code{FALSE}.
##' @return A long data.frame with sample,method,score values per row. If
##'   \code{as.matrix=TRUE}, a matrix with as many rows as \code{geneSets(gdb)}
##'   and as many columns as \code{ncol(x)}
##'
##' @examples
##' library(reshape2)
##' gdb <- exampleGeneSetDb()
##' vm <- exampleExpressionSet()
##' scores <- scoreSingleSamples(gdb, vm, methods=c('ewm', 'ssgsea', 'zscore'),
##'                              uncenter=FALSE, unscale=FALSE,
##'                              ssgsea.norm=TRUE)
##' sw <- dcast(scores, name + sample ~ method, value.var='score')
##' corplot(sw[, c("ewm", "svd", "zscore")],
##'         title='Single Sample Score Comparison')
scoreSingleSamples <- function(gdb, y, methods='ewm', as.matrix=FALSE,
                               drop.sd=1e-4, verbose=FALSE, ..., as.dt=FALSE) {
  methods <- tolower(methods)
  bad.methods <- setdiff(methods, names(gs.score.map))
  if (length(bad.methods)) {
    stop("Uknown geneset scoring methods: ",
         paste(bad.methods, collapse=','),
         "\nValid methods are: ",
         paste(names(gs.score.map), collapse=','))
  }
  ## TODO: Enable dispatch on whatever `method`s user asks for
  stopifnot(is(gdb, 'GeneSetDb'))
  y <- as_matrix(y)

  ## Removing genes that have almost-zero std.dev across the dataset.
  ## sds <- apply(y, 1, sd, na.rm=TRUE)
  sds <- rowSds(y)
  sd0 <- sds < drop.sd
  y.all <- y
  y <- y.all[!sd0,,drop=FALSE]
  if (any(sd0)) {
    warning(sum(sd0), " row(s) removed from expression object (y) due to 0sd")
  }
  gdb <- conform(gdb, y, ...)

  if (is.null(colnames(y))) {
    colnames(y) <- if (ncol(y) == 1) 'score' else paste0('scores', seq(ncol(y)))
  }

  gs.names <- encode_gskey(geneSets(gdb))
  gs.idxs <- as.list(gdb, active.only=TRUE, value='x.idx')

  scores <- sapply(methods, function(method) {
    fn <- gs.score.map[[method]]
    out <- fn(gdb, y, method=method, as.matrix=as.matrix, verbose=verbose,
              gs.idxs=gs.idxs, ...)
    rownames(out) <- gs.names
    if (!as.matrix) {
      out <- melt.gs.scores(gdb, out)
      out$method <- method
    }
    out
  }, simplify=FALSE)

  if (length(scores) == 1L) {
    scores <- scores[[1L]]
  } else {
    if (!as.matrix) {
      scores <- rbindlist(scores)
    }
  }

  if (!as.matrix && !as.dt) setDF(scores)

  # I'm not a bad person, I just want to keep this S3 so end users can
  # use the data.frame results in dplyr chains.
  # if (is.matrix(scores)) {
  #   class(scores) <- c('sss_matrix', class(scores))
  # } else {
  #   class(scores) <- c('sss_frame', class(scores))
  # }

  scores
}

# The \code{singleSampleScores} function is the same as
# \code{scoreSingleSamples}, but switches the order of the first two arguments
# for easier piping.
#
# @rdname scoreSingleSamples
# @export
# @inheritParams scoreSingleSamples
# singleSampleScores <- function(y, gdb, methods='ssgsea', as.matrix=FALSE,
#                                drop.sd=1e-4, verbose=FALSE, ...) {
#   scoreSingleSamples(gdb, y, methods=methods, as.matrix=as.matrix,
#                      drop.sd=drop.sd, verbose=verbose, ...)
# }

##' Melts the geneset matrix scores from the do.scoreSingleSamples.* methods
##'
##' @param gdb \code{GeneSetDb} used for scoring
##' @param scores The \code{matrix} of geneset scores returned from the various
##'   \code{do.scoreSingleSamples.*} methods.
##' @param a melted \code{data.table} of scores
melt.gs.scores <- function(gdb, scores) {
  out <- cbind(geneSets(gdb, as.dt=TRUE)[, list(collection, name, n)], scores)
  out <- data.table::melt.data.table(out, c('collection', 'name', 'n'),
                                     variable.name='sample',
                                     value.name='score')
  out[, sample := as.character(sample)]
}

## Default to sqrt in denominator of zscores to stabilize the variance of
## the mean:
##
## Lee, E., et al. Inferring pathway activity toward precise disease
## classification. PLoS Comput. Biol. 4, e1000217 (2008).
##
##
do.scoreSingleSamples.zscore <- function(gdb, y, zsummary=c('mean', 'sqrt'),
                                         trim=0.10, gs.idxs=NULL, do.scale=TRUE,
                                         ...) {
  stopifnot(is.conformed(gdb, y))
  zsummary <- match.arg(zsummary)

  score.fn <- if (zsummary == 'mean') {
    function(vals) mean(vals, trim=trim, na.rm=TRUE)
  } else {
    function(vals) {
      keep <- !is.na(vals)
      sum(vals[keep]) / sqrt(sum(keep))
    }
  }

  if (is.null(gs.idxs)) gs.idxs <- as.list(gdb, active.only=TRUE, value='x.idx')
  if (do.scale) y <- t(scale(t(y)))

  scores <- sapply(1:ncol(y), function(y.col) {
    col.vals <- y[, y.col]
    sapply(seq(gs.idxs), function(gs.idx) {
      vidx <- gs.idxs[[gs.idx]]
      score.fn(col.vals[vidx])
    })
  })
  colnames(scores) <- colnames(y)

  scores
}

## Just take the average of the raw scores from the expression matrix
##
## Right or wrong, sometimes you want this (most often you want mean Z, though)
do.scoreSingleSamples.mean <- function(gdb, y, gs.idxs=NULL, ...) {
  do.scoreSingleSamples.zscore(gdb, y, zsummary='mean',
                               trim=0, gs.idxs=gs.idxs, do.scale=FALSE)
}

.xformGdbForGSVA <- function(gdb, y) {
  stopifnot(is.conformed(gdb, y))
  stopifnot(is.matrix(y))
  out <- lapply(as.list(gdb, value='x.idx'), function(i) {
    rownames(y)[i]
  })
  out
}

##' @importFrom GSVA gsva
do.scoreSingleSamples.gsva <- function(gdb, y, method, as.matrix=FALSE,
                                       parallel.sz=4, ssgsea.norm=FALSE,
                                       gs.idxs=NULL, ...) {
  # idxs <- .xformGdbForGSVA(gdb, y)
  if (is.null(gs.idxs)) {
    gs.idxs <- as.list(gdb, active.only=TRUE, value='x.idx')
  }
  idxs <- lapply(gs.idxs, function(i) rownames(y)[i])
  f <- formals(GSVA:::.gsva)
  args <- list(...)
  ## I want to explicity show that we are setting parallel.sz to 4 here, since
  ## it will defalut to "Infinity" (all your cores are belong to GSVA)
  args$parallel.sz <- parallel.sz
  args$ssgsea.norm <- ssgsea.norm
  take <- intersect(names(args), names(f))
  gargs <- list(expr=y, gset.idx.list=idxs, method=method)
  gargs <- c(gargs, args[take])

  gres <- do.call(gsva, gargs)
  if (is.list(gres)) {
    gres <- gres$es.obs
  }

  gres
}

##' Normalize a vector of ssGSEA scores in the ssGSEA way.
##'
##' ssGSEA normalization (as implemented in GSVA (ssgsea.norm)) normalizes the
##' individual scores based on ALL scores calculated across samples AND
##' genesets. It does NOTE normalize the scores within each geneset
##' independantly of the others.
##'
##' @param x a \code{numeric} vector of ssGSEA scores for a single signature
##' @param bounds the maximum and minimum scores obvserved used to normalize
##'   against.
##' @return normalized \code{numeric} vector of \code{x}
ssGSEA.normalize <- function(x, bounds=range(x)) {
  ## apply(es, 2, function(x, es) x / (range(es)[2] - range(es)[1]), es)
  stopifnot(length(bounds) == 2L)
  max.b <- max(bounds)
  min.b <- min(bounds)
  stopifnot(all(x <= max.b) && all(x >= min.b))
  x / (max.b - min.b)
}

## A no dependency call to GSDecon-like eigengene scoring
do.scoreSingleSamples.gsd <- function(gdb, y, as.matrix=FALSE, center=TRUE,
                                      scale=TRUE, uncenter=center,
                                      unscale=scale, gs.idxs=NULL, ...) {
  stopifnot(is.matrix(y))
  stopifnot(is.conformed(gdb, y))

  if (is.null(gs.idxs)) {
    gs.idxs <- as.list(gdb, active.only=TRUE, value='x.idx')
  }

  scores <- lapply(gs.idxs, function(idxs) {
    gsdScore(y[idxs,], center=center, scale=scale, uncenter=uncenter,
             unscale=unscale)
  })

  out <- t(sapply(scores, '[[', 'score'))
  rownames(out) <- names(gs.idxs)
  colnames(out) <- colnames(y)
  out
}

do.scoreSingleSamples.eigenWeightedMean <- function(gdb, y, eigengene=1L,
                                                    center=TRUE, scale=TRUE,
                                                    uncenter=center,
                                                    unscale=scale,
                                                    as.matrix=FALSE,
                                                    gs.idxs=NULL, ...) {
  stopifnot(is.matrix(y))
  stopifnot(is.conformed(gdb, y))

  if (is.null(gs.idxs)) {
    gs.idxs <- as.list(gdb, active.only=TRUE, value='x.idx')
  }

  scores <- sapply(gs.idxs, function(idxs) {
    eigenWeightedMean(y[idxs,], center=center, scale=scale,
                      uncenter=uncenter, unscale=unscale)$score
  })

  out <- t(scores)
  rownames(out) <- names(gs.idxs)
  colnames(out) <- colnames(y)
  out
}

gs.score.map <- list(
  zscore=do.scoreSingleSamples.zscore,
  gsva=do.scoreSingleSamples.gsva,
  plage=do.scoreSingleSamples.gsva,
  ssgsea=do.scoreSingleSamples.gsva,
  svd=do.scoreSingleSamples.gsd,
  gsd=do.scoreSingleSamples.gsd,
  ewm=do.scoreSingleSamples.eigenWeightedMean,
  mean=do.scoreSingleSamples.mean)
