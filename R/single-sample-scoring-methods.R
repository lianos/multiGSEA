##' Single sample gene set score by a weighted average of the genes in geneset
##'
##' Weights for the genes in \code{x} are calculated by the percent of which
##' they contribute to the principal component indicated by \code{eigengene}.
##'
##' You will generally want the rows of the gene x sample matrix \code{x} to
##' be z-transformed. If it is not already, ensure that \code{center} and
##' \code{scale} are set to \code{TRUE}.
##'
##' When uncenter and/or unscale are \code{FALSE}, it means that the scores
##' should be applied on the centered or scaled values, respectively.
##'
##' @export
##' @importFrom stats weighted.mean
##' @inheritParams gsdScore
##' @seealso scoreSingleSamples
##' @param eigengene the PC used to extract the gene weights from
##' @param weights a user can pass in a prespecified set of waits using a named
##'   numeric vector. The names must be a superset of \code{rownames(x)}. If
##'   this is \code{NULL}, we calculate the "eigenweights".
##' @return A list of useful transformation information. The caller is likely
##'   most interested in the \code{$score} vector, but other bits related to
##'   the SVD/PCA decomposition are included for the ride.
##' @examples
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- conform(getMSigGeneSetDb('h'), vm)
##' features <- featureIds(gdb, 'h', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
##'                        value='x.idx')
##' scores <- eigenWeightedMean(vm[features,])$score
##'
##' ## Use scoreSingleSamples to facilitate scoring of all gene sets
##' scores.all <- scoreSingleSamples(gdb, vm, 'ewm')
##' s2 <- with(subset(scores.all, name == 'HALLMARK_INTERFERON_GAMMA_RESPONSE'),
##'            setNames(score, sample))
##' all.equal(s2, scores) ## should be TRUE
eigenWeightedMean <- function(x, eigengene=1L, center=TRUE, scale=TRUE,
                              uncenter=center, unscale=scale, retx=FALSE,
                              weights=NULL, ...) {
  x <- as_matrix(x)

  if (is.numeric(weights)) {
    if (length(weights) == 1) {
      weights <- setNames(rep(weights, nrow(x)), rownames(x))
    }
    stopifnot(all(rownames(x) %in% names(weights)))
    res <- list(weights=weights[rownames(x)])
  } else {
    pc <- paste0('PC', eigengene)
    res <- gsdScore(x, eigengene, center, scale, uncenter=uncenter,
                    unscale=unscale, retx=FALSE)
    res[['weights']] <- setNames(res$factor.contrib[[pc]], rownames(x))
  }
  if (!uncenter || !unscale) {
    xx <- t(scale(t(x), center=!uncenter, scale=!unscale))
  } else {
    xx <- x
  }
  res[['score']] <- apply(xx, 2, weighted.mean, res[['weights']])
  res
}

##' Calculate single sample geneset score by average z-score method
##'
##' @export
##' @param x gene x sample matrix
##' @param summary sqrt or mean
##' @param trim calculate trimmed mean?
##' @examples
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- conform(getMSigGeneSetDb('h'), vm)
##' features <- featureIds(gdb, 'h', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
##'                        value='x.idx')
##' scores <- zScore(vm[features,])$score
##'
##' ## Use scoreSingleSamples to facilitate scoring of all gene sets
##' scores.all <- scoreSingleSamples(gdb, vm, 'zscore')
##' s2 <- with(subset(scores.all, name == 'HALLMARK_INTERFERON_GAMMA_RESPONSE'),
##'            setNames(score, sample))
##' all.equal(s2, scores) ## should be TRUE
zScore <- function(x, summary=c('mean', 'sqrt'), trim=0, ...) {
  x <- as_matrix(x)
  summary <- match.arg(summary)
  score.fn <- if (summary == 'mean') {
    function(vals) mean(vals, trim=trim, na.rm=TRUE)
  } else {
    function(vals) {
      keep <- !is.na(vals)
      sum(vals[keep]) / sqrt(sum(keep))
    }
  }

  xs <- t(scale(t(x)))
  scores <- apply(xs, 2L, score.fn)

  out <- list(
    score=setNames(as.vector(scores), colnames(x)),
    center=attributes(xs)$"scaled:center",
    scale=attributes(xs)$"scaled:scale")
  out
}

##' Single sample geneset score using SVD based eigengene value per sample.
##'
##' @description
##' This method was developed by Jason Hackney and first introduced in the
##' following paper \href{https://doi.org/10.1038/ng.3520}{doi:10.1038/ng.3520}.
##' It produces a single sample gene set score in values that are in
##' "expression space," the innards of which mimic something quite similar
##' to an eigengene based score.
##'
##' To easily use this method to score a number of gene setes across an
##' experiment, you'll want to have the \code{\link{scoreSingleSamples}} method
##' drive this function via specifying \code{"svd"} as one of the
##' \code{methods}.
##'
##' @details
##' The difference between this method vs the eigengene score is that the SVD is
##' used to calculate the eigengene. The vector of eigengenes (one score per
##' sample) is then multiplied through by the SVD's left matrix. This produces a
##' matrix which we then take the colSums of to get back to a single sample
##' score for the geneset.
##'
##' Why do all of that? You get data that is back "in expression space" and we
##' also run around the problem of sign of the eigenvector. The scores you get
##' are very similar to average zscores of the genes per sample, where the
##' average is weighted by the degree to which each gene contributes to the
##' principal component chosen by \code{eigengene}, as implemented in the
##' \code{\link{eigenWeightedMean}} function.
##'
##' \emph{The core functionality provided here is taken from the soon to be
##' released GSDecon package by Jason Hackney}
##'
##' @export
##' @param x An expression matrix of genes x samples. When using this to score
##'   geneset activity, you want to reduce the rows of \code{x} to be only the
##'   genes from the given gene set.
##' @param eigengene the "eigengene" you want to get the score for. only accepts
##'   a single value for now.
##' @param center,scale center and/or scale data before scoring?
##' @param uncenter,unscale uncenter and unscale the data data on the way out?
##'   Defaults to the respective values of \code{center} and \code{scale}
##' @param retx Works the same as `retx` from \code{\link[stats]{prcomp}}. If
##'   \code{TRUE}, will return a \code{ret$pca$x} matrix that has the rotated
##'   variables.
##' @return A list of useful transformation information. The caller is likely
##'   most interested in the \code{$score} vector, but other bits related to
##'   the SVD/PCA decomposition are included for the ride.
##'
##' @examples
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- conform(getMSigGeneSetDb('h'), vm)
##' features <- featureIds(gdb, 'h', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
##'                        value='x.idx')
##' scores <- gsdScore(vm[features,])$score
##'
##' ## Use scoreSingleSamples to facilitate scoring of all gene sets
##' scores.all <- scoreSingleSamples(gdb, vm, 'gsd')
##' s2 <- with(subset(scores.all, name == 'HALLMARK_INTERFERON_GAMMA_RESPONSE'),
##'            setNames(score, sample))
##' all.equal(s2, scores) ## should be TRUE
gsdScore <- function(x, eigengene=1L, center=TRUE, scale=TRUE,
                     uncenter=center, unscale=scale, retx=FALSE, ...) {
  eigengene <- as.integer(eigengene)
  stopifnot(!is.na(eigengene) && length(eigengene) == 1L)
  x <- as_matrix(x)
  xs <- t(scale(t(x), center=center, scale=scale))

  cnt <- attributes(xs)$"scaled:center"
  scl <- attributes(xs)$"scaled:scale"

  s <- svd(xs)
  newD <- s$d
  newD[-eigengene] <- 0
  s$D <- diag(newD)

  egene <- s$u %*% s$D %*% t(s$v)

  if (unscale) {
    egene <- sweep(egene, 1, FUN="*", scl)
  }
  if (uncenter) {
    egene <- sweep(egene, 1, FUN="+", cnt)
  }

  score <- colMeans(egene)
  names(score) <- colnames(x)

  ## Reconstruct some PCA output here just because we've already done the heavy
  ## lifting calculations
  pca.d <- s$d / sqrt(max(1L, nrow(x) - 1L))
  pca.v <- s$v
  dimnames(pca.v) <- list(colnames(x), paste0('PC', seq_len(ncol(s$v))))

  ## Make this look ilke the output of `prcomp`
  ## s$v is the matrix with the eigenvectors. these entries should be correlated
  ## to our svd scores.
  xproj <- xs %*% s$v
  ## Show something like the percent contribution of each gene to the PCs
  ## This code was inspired from the biplot function, where the vectors of
  ## a PCA biplot are drawn as a function of the projected data (ie. pca$x)
  axproj <- abs(xproj)
  ctrb <- sweep(axproj, 2, colSums(axproj), '/')
  colnames(ctrb) <- paste0('PC', seq(ncol(ctrb)))

  rn <- if (!is.null(rownames(ctrb))) rownames(ctrb) else paste0('r', 1:nrow(ctrb))
  ctrb <- cbind(
    data.frame(featureId=rn, stringsAsFactors=FALSE),
    as.data.frame(ctrb))

  pca <- list(sdev=pca.d, rotation=pca.v,
              center=if (center) cnt else 0,
              scale=if (scale) scl else 1,
              x=if (retx) xproj else NULL,
              percentVar=pca.d^2 / sum(pca.d^2))
  class(pca) <- 'prcomp'

  list(score=score, egene=egene,
       svd=s, pca=pca, factor.contrib=ctrb)
}
