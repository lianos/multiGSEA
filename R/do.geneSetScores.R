do.geneSetScores <- function(x, design, contrast, gs.table,
                             logFC.stats=NULL,
                             robust.fit=FALSE, robust.eBayes=FALSE, ...) {
  if (is.null(logFC.stats)) {
    logFC.stats <- calculateIndividualLogFC(x, design, contrast, robust.fit,
                                            robust.eBayes, ...)
  } else {
    if (!is.numeric(logFC.stats) && length(stats) != nrow(x)) {
      stop("Illegal passed in value for logFC.stats")
    }
    if (!all(names(logFC.stats) == rownames(x))) {
      stop("names on passed-in logFC.stats do not match rownames(x")
    }
  }

  scores <- lapply(gs.table@table$membership, function(idx) {
    data.table(mean.stat=mean(logFC.stats[idx], na.rm=TRUE),
               mean.stat.t=mean(logFC.stats[idx], na.rm=TRUE, trim=0.10))
  })
  scores <- rbindlist(scores)

  cbind(gs.table@table[, list(group, id)], scores)
}

