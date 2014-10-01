calculateIndividualLogFC <- function(x, design, contrast, robust.fit=FALSE,
                                     robust.eBayes=FALSE, ...) {
  if (ncol(x) > 1) {
    ## Do limma fits on this thing and let it rip
    fit <- lmFit(x, design, method=if (robust.fit) 'robust' else 'ls', ...)
    fit  <- eBayes(fit, robust=robust.eBayes)
    stats <- topTable(fit, contrast, number=Inf, sort.by='none')$t
    names(stats) <- rownames(x)
  } else {
    stats <- x[, 1L]
    names(stats) <- rownames(x)
  }
  stats
}
