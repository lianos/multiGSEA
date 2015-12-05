##' Score each genes over the columns of an expression matrix.
##'
##' These functions provide gene set scores at the sample level. No need for
##' designs and contrasts.
##'
##' Current scoring methods include:
##'
##' \enumerate{
##'   \item zscore: features in geneset rowwise z transformed, then normalized
##'         by size of geneset (or sqrt(size))
##'   \item gsva: GSVA package
##'   \item plage: from GSVA package
##'   \item ssgsea: from GSVA package
##'   \item gsdecon: Jason Hackeny's method (kind of like plage), but not
##'         really.
##'   \item svd: Simple score that (should) mimic Jason's SVD/eigengene score.
##'         This is included here mainly to facilitate use of this scoring
##'         method for people who would have a hard time installing GSDecon2
##' }
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
scoreSingleSamples <- function(gdb, y, methods='ssgsea', melted=FALSE, ...) {
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
  if (is.vector(y)) {
    y <- t(t(y)) ## column vectorization that sets names to rownames
  }
  stopifnot(is.matrix(y) && is.numeric(y))

  ## Removing genes that have almost-zero std.dev across the dataset.
  sds <- apply(y, 1, sd, na.rm=TRUE)
  sd0 <- sds < 1e-4
  y.all <- y
  y <- y.all[!sd0,]
  if (any(sd0)) {
    warning(sum(sd0), " rows from expression objet (y) due to 0sd")
  }
  gdb <- conform(gdb, y, ...)

  if (is.null(colnames(y))) {
    colnames(y) <- if (ncol(y) == 1) 'score' else paste0('scores', seq(ncol(y)))
  }

  gs.names <- with(geneSets(gdb, .external=FALSE), {
    paste(collection, name, sep=';;')
  })
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
  out <- cbind(geneSets(gdb, .external=FALSE)[, list(collection, name, n)],
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

  gs <- geneSets(gdb, .external=FALSE)
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

  ## ssGSEA normalization:
  ## apply(es, 2, function(x, es) x / (range(es)[2] - range(es)[1]), es)
  gres
}


##' Jason's method
do.scoreSingleSamples.gsdecon <- function(gdb, y, melted=FALSE, design=NULL,
                                          doPerm=FALSE, nPerm=249,
                                          pvalueCutoff=0.01, nComp=1, seed=NULL,
                                          ...) {
  if (!requireNamespace('GSDecon')) {
    stop("Jason Hackney's GSDecon package required")
  }
  if (is.null(design)) {
    design <- cbind(Intercept=rep(1L, ncol(y)))
  }
  stopifnot(is.matrix(design))
  stopifnot(nrow(design) == ncol(y))
  im <- incidenceMatrix(gdb)
  res <- GSDecon::decon(y, design, im, doPerm=doPerm, nPerm=nPerm,
                        pvalueCutoff=pvalueCutoff, nComp=nComp, seed=seed)
  out <- t(res@eigengenes)
  rownames(out) <- rownames(im)
  out
}

##' A no dependency call to GSDecon's eigengene scoring
##'
##' TODO: Need to debug this. For some reason the results aren't the same
##' as what GSDecon methos provides, even though the codepath should be
##' identical
do.scoreSingleSamples.svd <- function(gdb, y, melted=FALSE, uncenter=TRUE,
                                      unscale=TRUE, ...) {
  ## if (!requireNamespace('GSDecon')) {
  ##   stop("Jason Hackney's GSDecon package required")
  ## }
  stopifnot(is.matrix(y))
  stopifnot(is.conformed(gdb, y))

  y.scale <- t(scale(t(y)))
  cnt <- if (uncenter) {
    attributes(y.scale)$"scaled:center"
  } else {
    rep(0, nrow(y))
  }
  scl <- if (unscale) {
    attributes(y.scale)$"scaled:scale"
  } else {
    rep(1, nrow(y))
  }
  message("=== SVD ===")
  gs.idxs <- as.list(gdb, nested=FALSE, value='x.idx')
  out <- sapply(gs.idxs, function(idxs) {
    ## idxs <- setdiff(idxs, sd0.idx)
    SVD <- svd(y.scale[idxs,])
    ## GSDecon:::eigencomponent(SVD, n=1, center=cnt[idxs], scale=scl[idxs])
    n <- 1
    center <- cnt[idxs]
    scale <- scl[idxs]
    newD <- SVD$d
    newD[-n] <- 0
    sv <- SVD$u %*% diag(newD) %*% t(SVD$v)
    sv <- sweep(sweep(sv, 1, FUN = "*", scale), 1, FUN = "+", center)
    sv <- colMeans(sv)
    sv
  })
  browser()
  out <- t(out)
  rownames(out) <- names(gs.idxs)
  colnames(out) <- colnames(y)
  out
}

gs.score.map <- list(
  zscore=do.scoreSingleSamples.zscore,
  gsva=do.scoreSingleSamples.gsva,
  plage=do.scoreSingleSamples.gsva,
  ssgsea=do.scoreSingleSamples.gsva,
  gsdecon=do.scoreSingleSamples.gsdecon,
  svd=do.scoreSingleSamples.svd)
