#' Single sample gene set score by a weighted average of the genes in geneset
#'
#' Weights for the genes in `x` are calculated by the percent of which
#' they contribute to the principal component indicated by `eigengene`.
#'
#' You will generally want the rows of the gene x sample matrix ``x` to
#' be z-transformed. If it is not already, ensure that `center` and
#' `scale` are set to `TRUE`.
#'
#' When uncenter and/or unscale are `FALSE`, it means that the scores
#' should be applied on the centered or scaled values, respectively.
#'
#' @section Normalization:
#' Scores can be normalized against a set of control genes. This results in
#' negative and postiive sample scores. Positive scores are ones where the
#' specific geneset score is higher than the aggregate control-geneset score.
#'
#' Genes used for the control set can either be randomly sampled from the
#' rows of the `all.x` expression matrix (when `normalize = TRUE`), or
#' explicitly specified by a row-identifier character vectore passed to the
#' `normalize` parameter. In both cases, the code prefers to select a
#' random-control geneset to be of equal size as `nrow(x)`. If that's not
#' possible, we use as many genes as we can get.
#'
#' Note that normalization requires an expression matrix to be passed into
#' the `all.x` parameter, whose columns match 1:1 to the columns in `x`.
#' Calling [scoreSingleSamples()] with `method = "ewm", normalize = TRUE`
#' handles this transparently.
#'
#' This idea to implement this method of normalizatition was inspried from
#' the `ctrl.score` normalization found in Seurat's
#' [Seurat::AddModuleScore()] function.
#'
#' @md
#' @export
#' @importFrom stats weighted.mean
#' @inheritParams gsdScore
#' @seealso scoreSingleSamples
#'
#' @param eigengene the PC used to extract the gene weights from
#' @param weights a user can pass in a prespecified set of waits using a named
#'   numeric vector. The names must be a superset of `rownames(x)`. If
#'   this is `NULL`, we calculate the "eigenweights".
#' @param normalize If `TRUE`, each score is normalized to a randomly
#'   selected geneset score. The size of the randomly selected geneset is
#'   the same as the corresponding geneset. This only works with the "ewm"
#'   method when unscale and uncenter are `TRUE`. By default, this is
#'   set to `FALSE`, and normalization does not happen. Instead of
#'   passing in `TRUE`, the user can pass in a vector of gene names
#'   (identifiers) to be considered for random geneset creation. If no
#'   genes are provided, then all genes in `y` are fair game.
#' @param all.x if the user is trying to normalize these scores, an expression
#'   matrix that has superset of the control genes needs to be provided, where
#'   the columns of `all.x` must correspond to this in `x`.
#' @return A list of useful transformation information. The caller is likely
#'   most interested in the `$score` vector, but other bits related to
#'   the SVD/PCA decomposition are included for the ride.
#' @examples
#' vm <- exampleExpressionSet(do.voom=TRUE)
#' gdb <- conform(getMSigGeneSetDb('h'), vm)
#' features <- featureIds(gdb, 'h', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
#'                        value='x.idx')
#' scores <- eigenWeightedMean(vm[features,])$score
#'
#' ## Use scoreSingleSamples to facilitate scoring of all gene sets
#' scores.all <- scoreSingleSamples(gdb, vm, 'ewm')
#' s2 <- with(subset(scores.all, name == 'HALLMARK_INTERFERON_GAMMA_RESPONSE'),
##'            setNames(score, sample))
#' all.equal(s2, scores) ## should be TRUE
eigenWeightedMean <- function(x, eigengene=1L, center=TRUE, scale=TRUE,
                              uncenter=center, unscale=scale, retx=FALSE,
                              weights=NULL, normalize = FALSE, all.x = NULL,
                              ..., .drop.sd = 1e-4) {
  x <- as_matrix(x)
  do.norm <- isTRUE(normalize) || (is.character(normalize) && length(normalize))

  if (is.numeric(weights)) {
    if (length(weights) == 1) {
      weights <- setNames(rep(weights, nrow(x)), rownames(x))
    }
    stopifnot(all(rownames(x) %in% names(weights)))
    res <- list(weights=weights[rownames(x)])
  } else {
    pc <- paste0('PC', eigengene)
    res <- gsdScore(x, eigengene, center, scale, uncenter=uncenter,
                    unscale=unscale, retx=FALSE, ..., .drop.sd = .drop.sd)
    res[['weights']] <- setNames(res$factor.contrib[[pc]], rownames(x))
  }
  if (!uncenter || !unscale) {
    xx <- t(scale(t(x), center=!uncenter, scale=!unscale))
  } else {
    xx <- x
  }
  res[['score']] <- apply(xx, 2, weighted.mean, res[['weights']])

  if (do.norm) {
    all.x <- try(as_matrix(all.x), silent = TRUE)
    if (!is.matrix(all.x)) {
      stop("`all.x` needs to be an expression matrix to normalize scores")
    }
    if (!isTRUE(all.equal(colnames(x), colnames(all.x)))) {
      stop("Must have 1:1 match for columns between `x` and `all.x`")
    }
    if (!(uncenter && unscale)) {
      warning("Normalization only works when uncenter & unscale are TRUE, ",
              "skipping normalization ...")
    } else {
      if (!is.character(normalize)) {
        normalize <- rownames(all.x)
      }
      normalize <- intersect(normalize, rownames(all.x))
      if (length(normalize) < nrow(x)) {
        warning(length(normalize) - nrow(x), " too few genes in all.x to use ",
                "for normalization: only using ", length(normalize), " genes")
      } else {
        normalize <- sample(normalize, nrow(x))
      }
      norm.x <- all.x[normalize,,drop = FALSE]
      norm.x <- colMeans(norm.x)
      res[["score"]] <- res[["score"]] - norm.x
    }
  }

  res
}

#' Just like eigenWeightedMean but direct from PCA (no gsdscore)
#'
#' @export
pcWeightedMean <- function(x, eigengene=1L, center=TRUE, scale=TRUE,
                           uncenter=center, unscale=scale, retx=FALSE,
                           weights=NULL, normalize = FALSE, all.x = NULL,
                           ..., .drop.sd = 1e-4) {
  stop("TODO: implement pcWeightedMean")
}

#' Calculate single sample geneset score by average z-score method
#'
#' @export
#' @param x gene x sample matrix
#' @param summary sqrt or mean
#' @param trim calculate trimmed mean?
#' @examples
#' vm <- exampleExpressionSet(do.voom=TRUE)
#' gdb <- conform(getMSigGeneSetDb('h'), vm)
#' features <- featureIds(gdb, 'h', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
#'                        value='x.idx')
#' scores <- zScore(vm[features,])$score
#'
#' ## Use scoreSingleSamples to facilitate scoring of all gene sets
#' scores.all <- scoreSingleSamples(gdb, vm, 'zscore')
#' s2 <- with(subset(scores.all, name == 'HALLMARK_INTERFERON_GAMMA_RESPONSE'),
#'            setNames(score, sample))
#' all.equal(s2, scores) ## should be TRUE
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

#' Single sample geneset score using SVD based eigengene value per sample.
#'
#' @description
#' This method was developed by Jason Hackney and first introduced in the
#' following paper [doi:10.1038/ng.3520](https://doi.org/10.1038/ng.3520).
#' It produces a single sample gene set score in values that are in
#' "expression space," the innards of which mimic something quite similar
#' to an eigengene based score.
#'
#' To easily use this method to score a number of gene setes across an
#' experiment, you'll want to have the [scoreSingleSamples()] method
#' drive this function via specifying `"svd"` as one of the
#' `methods`.
#'
#' @details
#' The difference between this method vs the eigengene score is that the SVD is
#' used to calculate the eigengene. The vector of eigengenes (one score per
#' sample) is then multiplied through by the SVD's left matrix. This produces a
#' matrix which we then take the colSums of to get back to a single sample
#' score for the geneset.
#'
#' Why do all of that? You get data that is back "in expression space" and we
#' also run around the problem of sign of the eigenvector. The scores you get
#' are very similar to average zscores of the genes per sample, where the
#' average is weighted by the degree to which each gene contributes to the
#' principal component chosen by `eigengene`, as implemented in the
#' [eigenWeightedMean()] function.
#'
#' *The core functionality provided here is taken from the soon to be
#' released GSDecon package by Jason Hackney*
#'
#' @export
#' @importFrom irlba svdr
#' @importFrom DelayedMatrixStats rowSds
#' @param x An expression matrix of genes x samples. When using this to score
#'   geneset activity, you want to reduce the rows of \code{x} to be only the
#'   genes from the given gene set.
#' @param eigengene the "eigengene" you want to get the score for. only accepts
#'   a single value for now.
#' @param center,scale center and/or scale data before scoring?
#' @param uncenter,unscale uncenter and unscale the data data on the way out?
#'   Defaults to the respective values of \code{center} and \code{scale}
#' @param retx Works the same as `retx` from \code{\link[stats]{prcomp}}. If
#'   \code{TRUE}, will return a \code{ret$pca$x} matrix that has the rotated
#'   variables.
#' @param ... these aren't used in here
#' @param .drop.sd When zero-sd (non varying) features are scaled, their values
#'   are `NaN`. When the Features with rowSds < this threshold (default 1e-4) are
#'   identified, and their scaled values are set to 0.
#' @return A list of useful transformation information. The caller is likely
#'   most interested in the \code{$score} vector, but other bits related to
#'   the SVD/PCA decomposition are included for the ride.
#'
#' @examples
#' vm <- exampleExpressionSet(do.voom=TRUE)
#' gdb <- conform(getMSigGeneSetDb('h'), vm)
#' features <- featureIds(gdb, 'h', 'HALLMARK_INTERFERON_GAMMA_RESPONSE',
#'                        value='x.idx')
#' scores <- gsdScore(vm[features,])$score
#'
#' ## Use scoreSingleSamples to facilitate scoring of all gene sets
#' scores.all <- scoreSingleSamples(gdb, vm, 'gsd')
#' s2 <- with(subset(scores.all, name == 'HALLMARK_INTERFERON_GAMMA_RESPONSE'),
#'            setNames(score, sample))
#' all.equal(s2, scores) ## should be TRUE
gsdScore <- function(x, eigengene = 1L, center = TRUE, scale = TRUE,
                     uncenter = center, unscale = scale, retx = FALSE, ...,
                     .use_irlba = FALSE, .drop.sd = 1e-4) {
  eigengene <- as.integer(eigengene)
  stopifnot(!is.na(eigengene) && length(eigengene) == 1L)
  x <- as_matrix(x)
  xs <- t(scale(t(x), center=center, scale=scale))

  # 0sd rows (genes) will deliver NaN's in the scaled matrix. This issue is
  # handled when we call scoreSingleSamples by removing 0sd rows, but when
  # this method is called directly, it will throw an error. In this case, we
  # will set the scaled values to 0.
  isnan <- is.nan(xs)
  if (scale && any(isnan)) {
    # Just making sure that the NaN entires are coming from rows with 0sd genes
    sds <- DelayedMatrixStats::rowSds(x)
    sd0 <- sds < .drop.sd
    sd0.idx <- which(sd0)
    if (length(sd0.idx)) {
      warning("Found NaN's in scaled matrix, replacing scaled values for ",
              "the ", length(sd0.idx), " features with 0-sd",
              immediate. = TRUE)
      xs[sd0.idx,] <- 0
      more.nan <- is.nan(xs)
      if (any(more.nan)) {
        stop("NaN features are found unrealted to 0-sd ... intervene, please")
      }
      # I tested adding a small amount of random noise to the features with
      # 0-sd, figuring that their contributions to the score will be
      # downweighted to zero and have minimal impact, but setting the scaled
      # values to 0 and testing seemed to perform slightly better
    } else {
      stop("NaN's found in scaled matrix for reasons unrelated to the ",
           "existence of 0-sd features")
    }
  } else if (any(isnan) || any(is.na(xs))) {
    stop("NaN's or NAs found in expression matrix without scaling")
  }

  cnt <- attributes(xs)$"scaled:center"
  scl <- attributes(xs)$"scaled:scale"

  if (.use_irlba) {
    s <- svdr(xs, k = min(eigengene, nrow(x), ncol(x)))
  } else {
    s <- svd(xs)
  }

  newD <- s$d
  newD[-eigengene] <- 0

  if (length(s$d) == 1L) {
    s$D <- matrix(s$d, nrow = 1L, ncol = 1L)
  } else {
    s$D <- diag(newD)
  }

  egene <- s$u %*% s$D %*% t(s$v)

  if (unscale) {
    egene <- sweep(egene, 1, FUN="*", scl)
  }
  if (uncenter) {
    egene <- sweep(egene, 1, FUN="+", cnt)
  }

  score <- colMeans(egene)
  names(score) <- colnames(x)

  # Reconstruct some PCA output here just because we've already done the heavy
  # lifting calculations
  pca.d <- s$d / sqrt(max(1L, nrow(x) - 1L))
  pca.v <- s$v
  dimnames(pca.v) <- list(colnames(x), paste0('PC', seq_len(ncol(s$v))))

  # Make this look ilke the output of `prcomp`
  # s$v is the matrix with the eigenvectors. these entries should be correlated
  # to our svd scores.
  xproj <- xs %*% s$v
  # Show something like the percent contribution of each gene to the PCs
  # This code was inspired from the biplot function, where the vectors of
  # a PCA biplot are drawn as a function of the projected data (ie. pca$x)
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
