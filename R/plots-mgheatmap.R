##' Create a "geneset smart" heatmap.
##'
##' @export
##' @importFrom circlize colorRamp2
##' @importFrom ComplexHeatmap Heatmap
##'
##' @param the \code{GeneSetDb} object that holds the genesets to plot
##' @param x the data matrix
##' @param col a colorRamp(2) funciton
##' @param aggregate.by the method used to generate single-sample geneset
##'   scores. Default is \code{none} which plots heatmap at the gene level
##' @param split.by introduce row-segmentation based on genesets or collections?
##'   Defaults to \code{'none'}
##' @param rm.dups if \code{aggregate.by == 'none'}, do we remove genes that
##'   appear in more than one geneset? Defaults to \code{FALSE}
##' @param recenter do you want to mean center the heatmap matrix prior to
##'   display?
##' @param rescale do you want to rescale the heatmap matrix prior to display?
##' @param ... parameters to send down to \code{\link{scoreSingleSample}} or
##'   \code{\link[ComplexHeatmap]{Heatmap}}.
##' @return list(heatmap=ComplexHeatmap, matrix=X)
##'
##' @examples
##'
##' vm <- exampleExpressionSet()
##' gdb <- exampleGeneSetDb()
##' mgh <- mgheatmap(gdb, vm, aggregate.by='ewm', split=TRUE)
mgheatmap <- function(gdb, x, col=NULL,
                      aggregate.by=c('none', 'ewm', 'zscore'),
                      # split.by=c('none', 'geneset', 'collection'),
                      split=TRUE,
                      name=NULL, rm.collection.prefix=TRUE,
                      rm.dups=FALSE, recenter=TRUE, rescale=TRUE,
                      ...) {
  stopifnot(is(gdb, "GeneSetDb"))
  # split.by <- match.arg(split.by)
  drop1.split <- missing(split)
  stopifnot(is.logical(split) && length(split) == 1L)

  X <- as_matrix(x)
  stopifnot(ncol(X) > 1)
  if (is.null(col)) {
    col <- colorRamp2(c(-2, 0, 2), c('#1F294E', 'white', '#6E0F11'))
  }
  stopifnot(is.function(col))
  aggregate.by <- match.arg(aggregate.by)

  gdbc <- conform(gdb, X, min.gs.size=2L)
  gdbc.df <- as.data.frame(gdbc) ## keep only genes that matched in gdb.df
  gdbc.df$key <- encode_gskey(gdbc.df)

  if (aggregate.by != 'none') {
    X <- scoreSingleSamples(gdb, X, methods=aggregate.by, as.matrix=TRUE, ...)
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
  H
}
