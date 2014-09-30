##' Runs camera on the given experiment.
##'
do.camera <- function(x, gs.table, design, contrast, outdir=NULL,
                      use.cache=TRUE, ...) {
  gs.idxs <- setNames(gs.table@table$membership,
                      paste(gs.table@table$group, gs.table@table$id, sep='.'))

  res <- camera(x, gs.idxs, design, contrast, ...)
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
