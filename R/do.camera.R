##' @include validateInputs.R
NULL

validate.inputs.camera <- .validate.inputs.full.design
validate.x.camera <- validate.X

##' Runs camera on the given experiment.
##'
##' @section Gene Set Enrichment Methods:
##'
##' camera tends to be very conservative, especially when trying to estimate
##' residual correlation of gene sets from experiments with small N. For this
##' reason,  Gordon has added the ability provide a prespecified
##' \code{inter.gene.correlation} value.
##' \href{https://support.bioconductor.org/p/70005/#70195}{
##' He suggests to try a small positive number (0.05)}
do.camera <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
                      use.cache=TRUE, weights=NULL, use.ranks=FALSE,
                      allow.neg.cor=TRUE, inter.gene.cor=NULL, trend.var=FALSE,
                      sort=FALSE,
                      ...) {
  stopifnot(is.conformed(gsd, x))
  if (!is.null(inter.gene.cor)) {
    ## Preset inter.gene.cor values were implemented in limma v3.24.14
    if (packageVersion('limma') < '3.24.14') {
      warning("inter.gene.cor values for camera require limma >- 3.24.14",
              immediate.=TRUE)
    }
  }
  extra.args <- c('use.ranks', 'allow.neg.cor', 'inter.gene.cor', 'trend.var',
                  'sort')
  cache.fn <- cache.data.fn('camera', design, contrast, extra.args,
                            outdir=outdir, ext='rds')
  gs.idxs <- as.expression.indexes(gsd, value='x.idx')

  if (file.exists(cache.fn) && use.cache) {
    res <- readRDS(cache.fn)
    ## TODO: check the genesets returned from pvals to ensure that they match
    ##       the active genesets
  } else {
    res <- camera(x, gs.idxs, design, contrast, weights=weights,
                  use.ranks=use.ranks, allow.neg.cor=allow.neg.cor,
                  inter.gene.cor=inter.gene.cor, trend.var=trend.var,
                  sort=FALSE, ...)
    if (is.character(outdir) && isTRUE(file.exists(outdir))) {
      saveRDS(res, cache.fn)
    }
  }

  out <- cbind(geneSets(gsd)[, list(collection, name)], as.data.table(res))
  out[, NGenes := NULL]
  setnames(out, c('PValue', 'FDR'), c('pval', 'padj'))
}
