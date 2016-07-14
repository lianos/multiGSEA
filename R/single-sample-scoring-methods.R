##' Calculates the (optinally) expression-rescaled eigengene for each sample.
##'
##' This function returns the innards of a PCA over "sample space"
##'
##' @export
##' @param x matrix of genes x samples
svdScore <- function(x, center=TRUE, scale=FALSE, uncenter=center,
                     unscale=scale, n=1L) {
  xs <- t(scale(t(x), center=center, scale=scale))

  cnt <- attributes(xs)$"scaled:center"
  scl <- attributes(xs)$"scaled:scale"

  s <- svd(xs)
  newD <- s$d
  newD[-n] <- 0
  s$D <- diag(newD)

  egene <- s$u %*% s$D %*% t(s$v)

  if (unscale) {
    egene <- sweep(egene, 1, FUN="*", scl)
  }
  if (uncenter) {
    egene <- sweep(egene, 1, FUN="+", cnt)
  }

  ## Reconstruct some PCA output here just because we've already done the heavy
  ## lifting calculations
  pca.d <- s$d / sqrt(max(1, nrow(x) - 1L))
  pca.v <- s$v
  dimnames(pca.v) <- list(colnames(x), paste0('PC', seq_len(ncol(s$v))))

  pca <- list(sdev=pca.d, rotation=s$v, center=if (center) cnt else NULL,
              scale=if (scale) scl else NULL,
              percentVar=pca.d^2 / sum(pca.d))

  list(score=colMeans(egene), egene=egene,
       svd=s, pca=pca,
       center=cnt, scale=scl)
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
