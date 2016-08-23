##' Calculates the (optinally) expression-rescaled eigengene for each sample.
##'
##' This function returns the innards of a PCA over "sample space"
##'
##' @export
##' @param x matrix of genes x samples
##' @param eigengene the "eigengene" you want to get the score for. only accepts
##'   a single value for now.
##' @param center center the data first?
##' @param scale scale the data first?
##' @param uncenter uncenter the data on the way out?
##' @param unscale the data on the way out?
##' @param retx Works the same as `retx` from \code{\link[stats]{prcomp}}. If
##'   \code{TRUE}, will return a \code{ret$pca$x} matrix that has the rotated
##'   variables.
##' @return list of useful transformation information
svdScore <- function(x, eigengene=1L, center=TRUE, scale=FALSE,
                     uncenter=FALSE, unscale=FALSE, retx=FALSE) {
  eigengene <- as.integer(eigengene)
  stopifnot(!is.na(eigengene) && length(eigengene) == 1L)
  x <- as.matrix(x)
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

##' Calculate geneset score by average z-score method
##'
##' @export
##' @param x gene x sample matrix
##' @param summary sqrt or mean
##' @param trim calculate trimmed mean?
zScore <- function(x, summary=c('sqrt', 'mean'), trim=0) {
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
    scores=as.vector(scores),
    center=attributes(xs)$"scaled:center",
    scale=attributes(xs)$"scaled:scale")
  out
}
