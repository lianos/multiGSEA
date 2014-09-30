do.roast <- function(x, gs.table, design, contrast, outdir=NULL,
                     use.cache=TRUE, nrot=10000, ..., .seed=NULL) {
  gs.idxs <- setNames(gs.table@table$membership,
                      paste(gs.table@table$group, gs.table@table$id, sep='.'))

  if (is.numeric(.seed)) {
    set.seed(.seed)
  }

  res <- mroast(x, gs.idxs, design, contrast, nrot=nrot, sort="none", ...)
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

