## Signature for geneSetTest is
##   index, statistics, alternative="mixed" ,type="auto",
##   ranks.only=TRUE, nsim=9999

do.gst <- function(x, gs.table, design, contrast,
                   alternative="mixed", type="auto", ranks.only=TRUE,
                   nsim=9999, robust.fit=FALSE, robust.eBayes=FALSE,
                   score.by=c('logFC', 't'), ...) {
  score.by <- match.arg(score.by)
  if (ncol(x) > 1) {
    stats <- calculateIndividualLogFC(x, design, contrast,
                                      robust.fit=robust.fit,
                                      robust.eBayes=robust.eBayes, ...)
    stats <- stats[[score.by]]
  } else {
    ## These are already precomputed logFC's (or whatever they are)
    stats <- x[, 1L]
  }

  results <- mclapply(gs.table@table$membership, function(idx) {
    geneSetTest(idx, stats, alternative, type, ranks.only, nsim)
  })

  out <- copy(gs.table@table)
  out[, pval := unlist(results)]
  out[, padj := p.adjust(pval, 'BH')]
  out[, padj.by.group := p.adjust(pval, 'BH'), by='group']

  f.ids <- featureNames(gs.table)
  out[, feature.id := lapply(membership, function(m) f.ids[m])]
}
