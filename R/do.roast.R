##' @include validateInputs.R
NULL

validate.inputs.roast <- .validate.inputs.full.design
validate.x.roast <- validate.X

## do.roast <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
##                      use.cache=TRUE, set.statistic='mean', gene.weights=NULL,
##                      array.weights=NULL, weights=NULL, block=NULL, correlation,
##                      var.prior=NULL, df.prior=NULL, trend.var=FALSE, nrot=10000,
##                      approx.zscore=TRUE, ...) {
do.roast <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
                     use.cache=TRUE, ...) {
  stopifnot(is.conformed(gsd, x))
  args <- list(...)
  call.args <- as.list(formals(limma::mroast.default))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  ## This sucks
  extra.args <- c('set.statistic', 'nrot', 'approx.zscore')
  set.statistic <- call.args$set.statistic
  nrot <- call.args$nrot
  approx.zscore <- call.args$approx.zscore

  cache.fn <- cache.data.fn('roast', design, contrast, extra.args,
                            outdir=outdir, ext='rds')
  if (file.exists(cache.fn) && use.cache) {
    res <- readRDS(cache.fn)
    ## TODO: check the genesets returned from pvals to ensure that they match
    ##       the active genesets
  } else {
    gs.idxs <- as.expression.indexes(gsd, value='x.idx')
    call.args[['y']] <- x
    call.args[['index']] <- gs.idxs
    call.args[['design']] <- design
    call.args[['contrast']] <- contrast
    call.args[['sort']] <- 'none'
    call.args[['...']] <- NULL
    ## earlier versions of edgeR::mroast was double-passing var.prior
    if (packageVersion('edgeR') < '3.11') {
      call.args[['var.prior']] <- NULL
      call.args[['df.prior']] <- NULL
    }
    res <- do.call(limma::mroast, call.args)
    if (is.character(outdir) && isTRUE(file.exists(outdir))) {
      saveRDS(res, cache.fn)
    }
  }

  out <- cbind(geneSets(gsd)[, list(collection, name)], as.data.table(res))
  out[, NGenes := NULL]

  setnames(out,
           c('PValue', 'FDR', 'PValue.Mixed', 'FDR.Mixed'),
           c('pval', 'padj', 'pval.mixed', 'padj.mixed'))

  ## the `result` method adds the rank based on what user wants to rank by
  ## out[, rank := rank(pval)]
  ## out[, rank.up := Inf]
  ## out[, rank.down := Inf]
  ## out[Direction == 'Up', rank.up := rank(pval)]
  ## out[Direction == 'Down', rank.down := rank(pval)]

  out
}
