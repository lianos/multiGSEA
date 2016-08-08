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
##'   \item svd: Simple score that (should) mimic Jason's SVD/eigengene score.
##'         This is included here mainly to facilitate use of this scoring
##'         method for people who would have a hard time installing GSDecon2
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
##' @return \code{matrix} with as many rows as \code{geneSets(gdb)} and
##'   as many columns as \code{ncol(x)}
scoreSingleSamples <- function(gdb, y, methods='ssgsea', as.matrix=FALSE,
                               drop.sd=1e-4, verbose=FALSE, ...) {
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
  ## sds <- apply(y, 1, sd, na.rm=TRUE)
  sds <- rowSds(y)
  sd0 <- sds < drop.sd
  y.all <- y
  y <- y.all[!sd0,]
  if (any(sd0)) {
    warning(sum(sd0), " row(s) removed from expression object (y) due to 0sd")
  }
  gdb <- conform(gdb, y, ...)

  if (is.null(colnames(y))) {
    colnames(y) <- if (ncol(y) == 1) 'score' else paste0('scores', seq(ncol(y)))
  }

  gs.names <- with(geneSets(gdb, .external=FALSE), {
    paste(collection, name, sep=';;')
  })

  scores <- sapply(methods, function(method) {
    fn <- gs.score.map[[method]]
    out <- fn(gdb, y, method=method, as.matrix=as.matrix, verbose=verbose, ...)
    rownames(out) <- gs.names
    if (!as.matrix) {
      out <- ret.df(melt.gs.scores(gdb, out))
      out$method <- method
    }
    out
  }, simplify=FALSE)

  if (length(scores) == 1L) {
    scores <- scores[[1L]]
  } else {
    if (!as.matrix) {
      scores <- ret.df(rbindlist(scores))
    }
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
  out <- data.table:::melt.data.table(out, c('collection', 'name', 'n'),
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
do.scoreSingleSamples.zscore <- function(gdb, y, zsummary=c('sqrt', 'mean'),
                                         trim=0.10, as.matrix=FALSE, ...) {
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
do.scoreSingleSamples.gsva <- function(gdb, y, method, as.matrix=FALSE,
                                       parallel.sz=4, ssgsea.norm=FALSE, ...) {
  idxs <- .xformGdbForGSVA(gdb, y)
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
##' @export
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

# Jason's method
# do.scoreSingleSamples.gsdecon <- function(gdb, y, as.matrix=FALSE, design=NULL,
#                                           doPerm=FALSE, nPerm=249,
#                                           pvalueCutoff=0.01, nComp=1, seed=NULL,
#                                           ...) {
#   if (!requireNamespace('GSDecon')) {
#     stop("Jason Hackney's GSDecon package required")
#   }
#   if (is.null(design)) {
#     design <- cbind(Intercept=rep(1L, ncol(y)))
#   }
#   stopifnot(is.matrix(design))
#   stopifnot(nrow(design) == ncol(y))
#   im <- incidenceMatrix(gdb)
#   res <- GSDecon::decon(y, design, im, doPerm=doPerm, nPerm=nPerm,
#                         pvalueCutoff=pvalueCutoff, nComp=nComp, seed=seed)
#   out <- t(res@eigengenes)
#   rownames(out) <- rownames(im)
#   out
# }

##' A no dependency call to GSDecon's eigengene scoring
##'
##' TODO: Need to debug this. For some reason the results aren't the same
##' as what GSDecon methos provides, even though the codepath should be
##' identical
do.scoreSingleSamples.svd <- function(gdb, y, as.matrix=FALSE, center=TRUE,
                                      scale=TRUE, uncenter=center,
                                      unscale=scale, ...) {
  stopifnot(is.matrix(y))
  stopifnot(is.conformed(gdb, y))

  gs.idxs <- as.list(gdb, nested=FALSE, value='x.idx')
  scores <- lapply(gs.idxs, function(idxs) {
    svdScore(y[idxs,], center=center, scale=scale, uncenter=uncenter,
             unscale=unscale)
  })

  out <- t(sapply(scores, '[[', 'score'))
  rownames(out) <- names(gs.idxs)
  colnames(out) <- colnames(y)
  out
}

gs.score.map <- list(
  zscore=do.scoreSingleSamples.zscore,
  gsva=do.scoreSingleSamples.gsva,
  plage=do.scoreSingleSamples.gsva,
  ssgsea=do.scoreSingleSamples.gsva,
  # gsdecon=do.scoreSingleSamples.gsdecon,
  svd=do.scoreSingleSamples.svd)
