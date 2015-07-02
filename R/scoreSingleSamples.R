##' Score each genes over the columns of an expression matrix.
##'
##' @export
##'
##' @param gdb A GeneSetDb
##' @param y An expression matrix to score genesets against
##' @param score.fn The function used to summarize the genes in a geneset.
##'   If left \code{NULL}, defaults to \code{mean}
##' @param trim The amount to trim in the trimmed mean of the individual
##'   geneset genes/features.
##' @return \code{matrix} with as many rows as \code{geneSets(gdb)} and
##'   as many columns as \code{ncol(x)}
scoreSingleSamples <- function(gdb, y, methods='plage', melted=FALSE, ...) {
  methods <- tolower(methods)
  bad.methods <- setdiff(methods, names(gs.score.map))
  if (length(bad.methods)) {
    stop("Uknown geneset scoring methods: ",
         paste(bad.methods, collapse=','))
  }
  ## TODO: Enable dispatch on whatever `method`s user asks for
  stopifnot(is(gdb, 'GeneSetDb'))
  if (is.vector(y)) {
    y <- t(t(y)) ## column vectorization that sets names to rownames
  }
  stopifnot(is.matrix(y) && is.numeric(y))
  stopifnot(is.conformed(gdb, y))
  if (is.null(colnames(y))) {
    colnames(y) <- if (ncol(y) == 1) 'score' else paste0('scores', seq(ncol(y)))
  }

  gs.names <- with(geneSets(gdb), paste(collection, name, sep=';;'))
  scores <- sapply(methods, function(method) {
    out <- gs.score.map[[method]](gdb, y, method=method, melted=melted, ...)
    rownames(out) <- gs.names
    if (melted) {
      out <- melt.gs.scores(gdb, out)
    }
    out
  }, simplify=FALSE)

  if (length(scores) == 1L) {
    scores <- scores[[1L]]
  }

  scores
}

##' Melts the geneset matrix scores from the do.scoreSingleSamples.* methods
##'
##' @param gdb \code{GeneSetDb} used for scoring
##' @param scores The \code{matrix} of geneset scores returned from the various
##'   \code{do.scoreSingleSamples.*} methods.
##' @param a melted \code{data.table} of scores
melt.gs.scores <- function(gdb, scores) {
  out <- cbind(geneSets(gdb)[, list(collection, name, n)],
               scores)
  data.table:::melt.data.table(out, c('collection', 'name', 'n'),
                               variable.name='sample',
                               value.name='score')
}

## Default to sqrt in denominator of zscores to stabilize the variance of
## the mean:
##
## Lee, E., et al. Inferring pathway activity toward precise disease
## classification. PLoS Comput. Biol. 4, e1000217 (2008).
##
##
do.scoreSingleSamples.zscore <- function(gdb, y, zsummary=c('sqrt', 'mean'),
                                    trim=0.10, melted=FALSE, ...) {
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

  gs <- geneSets(gdb)
  gs.idxs <- lapply(1:nrow(gs), function(i) {
    featureIds(gdb, gs$collection[i], gs$name[i], 'x.idx')
  })
  y <- t(scale(t(y)))
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

.xformGdbForGSVA <- function(gdb, y) {
  stopifnot(is.conformed(gdb, y))
  stopifnot(is.matrix(y))
  out <- lapply(as.list(gdb, value='x.idx'), function(i) {
    rownames(y)[i]
  })
  out
}

##' @importFrom GSVA gsva
do.scoreSingleSamples.gsva <- function(gdb, y, method, melted=FALSE,
                                  tweak.plage.sign=FALSE, ...) {
  idxs <- .xformGdbForGSVA(gdb, y)
  f <- formals(GSVA:::.gsva)
  args <- list(...)
  take <- intersect(names(args), names(f))
  gargs <- list(expr=y, gset.idx.list=idxs, method=method)
  gargs <- c(gargs, args[take])

  normalize.ssGSEA <- isTRUE(gargs$ssgsea.norm)

  gres <- do.call(gsva, gargs)
  if (is.list(gres)) {
    gres <- gres$es.obs
  }

  if (method == 'plage' && tweak.plage.sign) {
    ## The sign of the result can be flipped due to vagaries of SVD assigning
    ## the "correct" sign to either the right or left singula values, so let's
    ## put some duct tape on that and fix the sign
    zscores <- do.scoreSingleSamples.zscore(gdb, y)
    gres <- abs(gres) * sign(zscores)
  }

  gres
}

gs.score.map <- list(
  zscore=do.scoreSingleSamples.zscore,
  gsva=do.scoreSingleSamples.gsva,
  plage=do.scoreSingleSamples.gsva,
  ssgsea=do.scoreSingleSamples.gsva)
