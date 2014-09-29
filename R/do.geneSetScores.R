do.geneSetScores <- function(x, gs.table, design, contrast,
                             robust.fit=FALSE, robust.eBayes=FALSE, ...) {
  if (ncol(x) > 1) {
    ## Do limma fits on this thing and let it rip
    fit <- lmFit(x, design, method=if (robust.fit) 'robust' else 'ls', ...)
    fit  <- eBayes(fit, robust=robust.eBayes)
    stats <- topTable(fit, contrast, number=Inf, sort.by='none')$t
  } else {
    stats <- x[, 1L]
  }

  scores <- lapply(gs.table@table$membership, function(idx) {
    data.table(mean.stat=mean(stats[idx], na.rm=TRUE),
               mean.stat.t=mean(stats[idx], na.rm=TRUE, trim=0.10))
  })
  scores <- rbindlist(scores)

  cbind(gs.table@table[, list(group, id)], scores)
}
