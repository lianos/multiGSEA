#' Creates a "geneset smart" ComplexHeatmap::Heatmap
#'
#' Encapsulates many common "moves" you'll make when trying to make a heatmap,
#' especially if you are trying to show geneset activity across a panel of
#' samples.
#'
#' @section Renaming Heatmap Rows:
#' You might want to rename the rows of the heatmap
#' "on the way out". For instance, `x` might be a `DGEList` whose
#' `rownames()` are gene identifieres, but you want the rownames of the
#' heatmap to gene symbols. You can perform this renaming using the
#' `rename.rows` parameter.
#'
#' * If `rename.rows` is `NULL`, then nothing is done.
#' * If `rename.rows` is a `string`, then we assume that `x` has an associated
#'   metadata `data.frame` over its rows and that `rename.rows` names one of
#'   its columns, ie. `DGEList$genes[[rename.rows]]` or
#'   `fData(ExpressionSet)[[rename.rows]]`. The values in that column will
#'   be swapped out for `x`'s rownames
#' * If `rename.rows` is a two-column data.frame, the first column is assumed
#'   to be `rownames(x)` and the second is what you want to rename it to.
#'
#' @md
#' @export
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom viridis viridis
#'
#' @param x the data matrix
#' @param gdb `GeneSetDb` object that holds the genesets to plot. Defaults to
#'   `NULL`, which will plot all rows in `x`.
#' @param col a colorRamp(2) function
#' @param aggregate.by the method used to generate single-sample geneset
#'   scores. Default is `none` which plots heatmap at the gene level
#' @param split introduce row-segmentation based on genesets or collections?
#'   Defaults is `TRUE` which will create split heatmaps based on
#'   collection if `aggregate.by != 'none'`, or based on gene sets
#'   if `aggregate.by == "none"`.
#' @param scores If `aggregate.by != "none"` you can pass in a precomupted
#'   [scoreSingleSamples()] result, otherwise one will be
#'   computed internally. Note that if this is a `data.frame` of
#'   pre-computed scores, the `gdb` is largely irrelevant (but still
#'   required).
#' @param rm.dups if `aggregate.by == 'none'`, do we remove genes that
#'   appear in more than one geneset? Defaults to `FALSE`
#' @param recenter do you want to mean center the rows of the heatmap matrix
#'   prior to calling [ComplexHeatmap::Heatmap()]?
#' @param rescale do you want to standardize the row variance to one on the
#'   values of the heatmap matrix prior to calling
#'   [ComplexHeatmap::Heatmap()]?
#' @param rename.rows defaults to `NULL`, which induces no action. A `string`
#'   can be specified which will do a 'smart lookup' on meta information on the
#'   rows of `x` to swap out that column for the rownames of `x` (ie. from a
#'   `DGEList$genes` or `fData(ExpressionSet)`. A two-column data.frame can
#'   also be provided where the first column is assumed to be `rownames(x)` and
#'   the second is the renamed stuff.
#' @param zlim A `length(zlim) == 2` numeric vector that defines the min and max
#'   values from `x` for the `colorRamp2` call. If the heatmap that is being
#'   drawn is "0-centered"-ish, then this defines the real values of the
#'   fenceposts. If not, then these define the quantiles to trim off the top
#'   or bottom.
#' @param transpose Flip display so that rows are columns. Default is `FALSE`.
#' @param ... parameters to send down to [scoreSingleSamples()] or
#'   [ComplexHeatmap::Heatmap()].
#' @return A `Heatmap` object.
#'
#' @examples
#'
#' vm <- exampleExpressionSet()
#' gdb <- exampleGeneSetDb()
#' col.anno <- ComplexHeatmap::HeatmapAnnotation(
#'   vm$targets[, c("Cancer_Status", "PAM50subtype")],
#'   col = list(
#'     Cancer_Status = c(normal = "grey", tumor = "red"),
#'     PAM50subtype = c(Basal = "purple", Her2 = "green", LumA = "orange")))
#' mgh <- mgheatmap(vm, gdb, aggregate.by='ewm', split=TRUE,
#'                  top_annotation = col.anno, show_column_names = FALSE,
#'                  column_title = "Gene Set Activity in BRCA subset")
mgheatmap <- function(x, gdb = NULL, col=NULL,
                      aggregate.by=c('none', 'ewm', 'zscore'),
                      split=TRUE, scores=NULL,
                      name=NULL, rm.collection.prefix=TRUE,
                      rm.dups=FALSE, recenter=TRUE, rescale=FALSE,
                      rename.rows = NULL, zlim = NULL, transpose = FALSE, ...) {
  X <- as_matrix(x)

  if (is.null(gdb)) {
    # make a one geneset GeneSetDb
    faux.gs <- list(allgenes = rownames(x))
    gdb <- GeneSetDb(faux.gs, collectionName = "faux")
  }
  stopifnot(is(gdb, "GeneSetDb"))

  # split.by <- match.arg(split.by)
  drop1.split <- missing(split)
  stopifnot(is.logical(split) && length(split) == 1L)
  if (!is.null(scores)) stopifnot(is.data.frame(scores))
  if (!missing(zlim) && !is.null(zlim)) {
    stopifnot(
      is.numeric(zlim),
      length(zlim) == 2L,
      zlim[1] < zlim[2])
  }

  stopifnot(
    ncol(X) > 1L,
    !any(is.na(X)))

  if (is.null(scores)) {
    aggregate.by <- match.arg(aggregate.by)
  } else {
    stopifnot(
      is.character(aggregate.by),
      length(aggregate.by) == 1L,
      aggregate.by %in% scores$method)
  }

  gdbc <- suppressWarnings(conform(gdb, X, ...))
  gdbc.df <- as.data.frame(gdbc) ## keep only genes that matched in gdb.df
  gdbc.df$key <- encode_gskey(gdbc.df)

  if (aggregate.by != 'none') {
    if (is.null(scores)) {
      X <- scoreSingleSamples(gdb, X, methods=aggregate.by, as.matrix=TRUE, ...)
    } else {
      xs <- scores[scores[['method']] == aggregate.by,,drop=FALSE]
      xs$key <- encode_gskey(xs)
      X <- acast(xs, key ~ sample, value.var="score")
    }
    ## If we want to split, it (only?) makes sense to split by collection
    split <- if (split) split_gskey(rownames(X))$collection else NULL
  }

  if (aggregate.by == 'none') {
    ridx <- if (rm.dups) unique(gdbc.df$featureId) else gdbc.df$featureId
    X <- X[ridx,,drop=FALSE]
    split <- if (split) gdbc.df$key else NULL
  }

  if (recenter || rescale) {
    X <- t(scale(t(X), center=recenter, scale=rescale))
    isna <- which(is.na(X), arr.ind = TRUE)
    if (nrow(isna) > 0L) {
      na.rows <- unique(isna[, "row"])
      if (length(na.rows) == nrow(X)) {
        stop("All rows removed after `scale`")
      }
      warning(length(na.rows), " features NA'd during `scale`, ",
              "these are removed", immediate. = TRUE)
      X <- X[-na.rows,,drop = FALSE]
      split <- split[-na.rows]
    }
  }

  # What kind of colorscale are we going to use?
  # If this is 0-centered ish, we use a red-white-blue scheme, otherwise
  # we use viridis.
  if (is.null(col)) {
    # Is 0 close to the center of the score distribution?
    mean.X <- mean(X)
    zero.center <- mean.X >= -0.2 && mean.X <= 0.2
    if (zero.center) {
      if (missing(zlim)) {
        fpost <- quantile(abs(X), 0.975)
        zlim <- c(-fpost, fpost)
      } else if (is.null(zlim)) {
        zlim <- c(min(X), max(X))
      } else {
        stopifnot(zlim[1L] < 0, zlim[2L] > 0)
      }
      col <- colorRamp2(
        c(zlim[1L], 0, zlim[2L]),
        c('#1F294E', 'white', '#6E0F11'))
    } else {
      if (missing(zlim)) {
        fpost <- quantile(X, c(0.025, 0.975))
      } else if (is.null(zlim)) {
        fpost <- c(min(X), max(X))
      } else {
        stopifnot(all(zlim >= 0), all(zlim <= 1))
        fpost <- quantile(X, zlim)
      }
      breaks <- quantile(X, seq(0, 1, by = 0.25))
      if (fpost[1L] > breaks[2L] || fpost[2L] < breaks[4L]) {
        stop("Illegal values for zlim")
      }
      breaks[1] <- fpost[1]
      breaks[5] <- fpost[2]
      col <- colorRamp2(breaks, viridis::viridis(5))
    }
  }
  stopifnot(is.function(col))

  if (drop1.split && !is.null(split) && length(unique(split)) == 1L) {
    split <- NULL
  }

  if (rm.collection.prefix) {
    if (aggregate.by != 'none') {
      rownames(X) <- split_gskey(rownames(X))$name
    } else {
      if (!is.null(split)) split <- split_gskey(split)$name
    }
  }

  ## Catch Heatmap arguments in `...` and build a list do do.call() them down
  ## into the function call.
  dot.args <- list(...)
  hm.args.default <- as.list(formals(Heatmap))

  if (is.null(name)) {
    name <- if (aggregate.by == 'none') 'value' else 'score'
  }
  hm.args <- dot.args[intersect(names(dot.args), names(hm.args.default))]
  hm.args[['matrix']] <- X
  hm.args[['col']] <- col
  hm.args[['split']] <- split
  hm.args[['name']] <- name

  H <- do.call(ComplexHeatmap::Heatmap, hm.args)

  if (is.data.frame(rename.rows)) {
    xref <- match(rownames(H@matrix), rename.rows[[1L]])
    vals <- rename.rows[[2L]][xref]
    isna <- is.na(vals) | nchar(vals) == 0L
    if (any(isna)) {
      msg <- sprintf("%d (%.2f%%) of rownames could not be renamed",
                     sum(isna), mean(isna))
      vals[isna] <- rownames(H@matrix)[isna]
    }
    rownames(H@matrix) <- vals
  }

  H
}

# mgheatmap <- function(x, ...) {
#   ## I'm not a bad person, I just want to keep this S3 so end users can
#   ## use the data.frame results in dplyr chains.
#   UseMethod("mgheatmap")
# }
#
# mgheatmap.sss_frame <- function(x, col=NULL, aggregate.by=x$method[1L],
#                                 split=TRUE, name=NULL,
#                                 rm.collection.prefix=TRUE, recenter=TRUE,
#                                 rescale=TRUE, ...) {
#   stopifnot(
#     is.character(aggregate.by),
#     length(aggregate.by) == 1L,
#     aggregate.by %in% x$method)
#
#   x$key <- encode_gskey(x)
#   xs <- subset(x, method == aggregate.by)
#   X <- acast(xs, key ~ sample, value.var="score")
#   mgheatmap(X, )
# }
#
# mgheatmap.default <- function()
