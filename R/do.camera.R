##' Runs camera on the given experiment.
##'
do.camera <- function(x, gs.table, design, contrast, outdir=NULL,
                      use.cache=TRUE, weights=NULL, use.ranks=FALSE,
                      allow.neg.cor=TRUE, trend.var=FALSE, sort=FALSE, ...) {
  extra.args <- c('allow.neg.cor', 'trend.var', 'sort')
  gs.idxs <- setNames(gs.table@table$membership,
                      paste(gs.table@table$group, gs.table@table$id, sep='.'))

  cache.fn <- local({
    extra.args <- list(use.ranks=use.ranks, allow.neg.cor=allow.neg.cor,
                       trend.var=trend.var, sort=sort)
    cache.data.fn('camera', design, contrast, extra.args=extra.args,
                  outdir=outdir, ext='rds')
  })

  if (file.exists(cache.fn) && use.cache) {
    res <- readRDS(cache.fn)
  } else {
    res <- camera(x, gs.idxs, design, contrast, weights=weights,
                  use.ranks=use.ranks, allow.neg.cor=allow.neg.cor,
                  trend.var=trend.var, sort=sort, ...)
    if (!is.null(outdir)) {
      saveRDS(res, cache.fn)
    }
  }

  xref <- match(rownames(res), names(gs.idxs))
  out <- cbind(gs.table@table[xref], as.data.table(res))
  stopifnot(all.equal(out$n, out$NGenes))
  out[, NGenes := NULL]
  setnames(out, c('PValue', 'FDR'), c('pval', 'padj'))

  out <- out[order(pval),]
  out[, padj.by.group := p.adjust(pval, 'BH'), by='group']

  f.ids <- featureNames(gs.table)
  out[, feature.id := lapply(membership, function(m) f.ids[m])]
}
