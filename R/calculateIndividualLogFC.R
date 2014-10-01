calculateIndividualLogFC <- function(x, design, contrast, robust.fit=FALSE,
                                     robust.eBayes=FALSE, use.t=TRUE, ...) {
  if (ncol(x) > 1) {
    ## Do limma fits on this thing and let it rip
    fit <- lmFit(x, design, method=if (robust.fit) 'robust' else 'ls', ...)
    if (length(contrast) > 1) {
      if (length(contrast) != ncol(design)) {
        stop("Invalid contrast vector, must be as long as columns in design")
      }
      fit <- contrasts.fit(fit, contrast)
      contrast <- 1
    }
    fit  <- eBayes(fit, robust=robust.eBayes)
    tt <- topTable(fit, contrast, number=Inf, sort.by='none')
    stats <- if (use.t) tt$t else tt$logFC
    names(stats) <- rownames(x)
  } else {
    stats <- x[, 1L]
    names(stats) <- rownames(x)
  }
  stats
}
