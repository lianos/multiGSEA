## Signature for geneSetTest is
##   index, statistics, alternative="mixed" ,type="auto",
##   ranks.only=TRUE, nsim=9999

do.gst <- function(x, gs.table, design, contrast,
                   alternative="mixed", type="auto", ranks.only=TRUE,
                   nsim=9999, robust.fit=FALSE, robust.eBayes=FALSE, ...) {
  if (ncol(x) > 1) {
    ## Do limma fits on this thing and let it rip
    fit <- lmFit(x, design, method=if (robust.fit) 'robust' else 'ls', ...)
    fit  <- eBayes(fit, robust=robust.eBayes)
    stats <- topTable(fit, contrast, number=Inf, sort.by='none')$t
  } else {
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
