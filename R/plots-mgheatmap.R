#' Create a "geneset smart" heatmap.
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
#'
#' @param the `GeneSetDb` object that holds the genesets to plot
#' @param x the data matrix
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
#' @param ... parameters to send down to [scoreSingleSamples()] or
#'   [ComplexHeatmap::Heatmap()].
#' @return A `Heatmap` object.
#'
#' @examples
#'
#' vm <- exampleExpressionSet()
#' gdb <- exampleGeneSetDb()
#' mgh <- mgheatmap(gdb, vm, aggregate.by='ewm', split=TRUE)
mgheatmap <- function(gdb, x, col=NULL,
                      aggregate.by=c('none', 'ewm', 'zscore'),
                      split=TRUE, scores=NULL,
                      name=NULL, rm.collection.prefix=TRUE,
                      rm.dups=FALSE, recenter=TRUE, rescale=TRUE,
                      rename.rows = NULL, ...) {
  stopifnot(is(gdb, "GeneSetDb"))
  # split.by <- match.arg(split.by)
  drop1.split <- missing(split)
  stopifnot(is.logical(split) && length(split) == 1L)
  if (!is.null(scores)) stopifnot(is.data.frame(scores))

  X <- as_matrix(x)
  stopifnot(ncol(X) > 1L)

  if (!is.null(rename.rows)) {
    stopifnot(is.character(rename.rows), length(rename.rows) == 1L)
    if (is.character(rename.rows)) {
      if (is(x, "DGEList")) {
        mdf <- x$genes
      } else if (is(x, "eSet")) {
        mdf <- Biobase::fData(x)
      } else if (is(x, "SummarizedExperiment")) {
        mdf <- as.data.frame(SummarizedExperiment::mcols(x))
      }
      rename <- mdf[[rename.rows]]
      if (!is.character(rename)) {
        stop("rename.rows column (", rename.rows, ") is not a character ",
             "metadata for x")
      }
      rename.rows <- data.frame(rn = rownames(x), rename = rename,
                                stringsAsFactors = FALSE)
    }
    stopifnot(
      is.data.frame(rename.rows),
      ncol(rename.rows) == 2L,
      is.character(rename.rows[[2L]]))
  }

  if (is.null(col)) {
    col <- colorRamp2(c(-2, 0, 2), c('#1F294E', 'white', '#6E0F11'))
  }
  stopifnot(is.function(col))

  if (is.null(scores)) {
    aggregate.by <- match.arg(aggregate.by)
  } else {
    stopifnot(
      is.character(aggregate.by),
      length(aggregate.by) == 1L,
      aggregate.by %in% scores$method)
  }

  gdbc <- suppressWarnings(conform(gdb, X, min.gs.size=2L))
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

  if (recenter || rescale) {
    X <- t(scale(t(X), center=recenter, scale=rescale))
  }

  if (aggregate.by == 'none') {
    ridx <- if (rm.dups) unique(gdbc.df$featureId) else gdbc.df$featureId
    X <- X[ridx,,drop=FALSE]
    split <- if (split) gdbc.df$key else NULL
  }

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
