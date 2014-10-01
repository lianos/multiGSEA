do.roast <- function(x, gs.table, design, contrast, outdir=NULL,
                     use.cache=TRUE, set.statistic='mean', gene.weight=NULL,
                     array.weights=NULL, weights=NULL, block=NULL, correlation,
                     var.prior=NULL, df.prior=NULL, trend.var=FALSE, nrot=10000,
                     approx.zscore=TRUE, ..., .seed=NULL) {
  extra.args <- c('set.statistic', 'gene.weight', 'array.weights', 'weights',
                  'block', 'correlation', 'var.prior', 'df.prior', 'trend.var',
                  'nrot', 'approx.zscore')
  gs.idxs <- setNames(gs.table@table$membership,
                      paste(gs.table@table$group, gs.table@table$id, sep='.'))

  if (is.numeric(.seed)) {
    set.seed(.seed)
  }

  cache.fn <- local({
    extra.args <- list(set.statistic=set.statistic, nrot=nrot,
                       approx.zscore=approx.zscore)
    cache.data.fn('roast', design, contrast, extra.args=extra.args,
                  outdir=outdir, ext='rds')
  })

  if (file.exists(cache.fn) && use.cache) {
    res <- readRDS(cache.fn)
  } else {
    res <- mroast(x, gs.idxs, design, contrast, set.statistic=set.statistic,
                  gene.weight=gene.weight, array.weights=array.weights,
                  weights=weights, block=block, correlation=correlation,
                  var.prior=var.prior, df.prior=df.prior, trend.var=trend.var,
                  nrot=nrot, approx.zscore=approx.zscore, ...)
    if (!is.null(outdir)) {
      saveRDS(res, cache.fn)
    }
  }

  xref <- match(rownames(res), names(gs.idxs))
  out <- cbind(gs.table@table[xref], as.data.table(res))
  stopifnot(all.equal(out$n, out$NGenes))
  out[, NGenes := NULL]

  setnames(out,
           c('PValue', 'FDR', 'PValue.Mixed', 'FDR.Mixed'),
           c('pval', 'padj', 'pval.mixed', 'padj.mixed'))

  out <- out[order(pval),]
  out[, padj.by.group := p.adjust(pval, 'BH'), by='group']
  out[, rank := 1:.N]
  out[, rank.up := NA_integer_]
  out[, rank.down := NA_integer_]
  out[Direction == 'Up', rank.up := seq(.N)]
  out[Direction == 'Down', rank.down := seq(.N)]

  f.ids <- featureNames(gs.table)
  out[, feature.id := lapply(membership, function(m) f.ids[m])]
}

