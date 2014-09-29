do.roast <- function(x, gs.table, design, contrast, ...) {
  lookup <- gs.table@feature.lookup
  lookup <- lookup[order(x.index)]
  f.ids <- rep(NA_character_, length(gs.table@table$membership[[1]]))
  if (length(f.ids) > max(lookup$x.index)) {
    stop("We have a problem with the lookup table")
  }
  f.ids[lookup$x.index] <- lookup$x.id

  res <- mroast(x, gs.table@table$membership, design, contrast, ...)
#
#   out <- cbind(gs.table@table[as.integer(rownames(res))],
#                as.data.table(res))
#   stopifnot(all.equal(out$n, out$NGenes))
#   setnames(out, c('PValue', 'FDR'), c('pval', 'padj'))
#
#   out <- out[order(pval),]
#   out[, feature.id := lapply(membership, function(m) f.ids[m])]
#   out[, padj.by.group := p.adjust(pval, 'BH'), by='group']
#   out
}
