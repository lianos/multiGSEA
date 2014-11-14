do.geneSetScores <- function(x, gs.table, design, contrast, logFC.stats=NULL,
                             robust.fit=FALSE, robust.eBayes=FALSE,
                             trim=0.10, ...) {
  if (is.null(logFC.stats)) {
    logFC.stats <- calculateIndividualLogFC(x, design, contrast, robust.fit,
                                            robust.eBayes, provide='table', ...)
  } else if (!is.data.frame(logFC.stats)) {
    stop("data.frame required for logFC.stats parameter")
  }
  missed.cols <- setdiff(c('logFC', 't'), names(logFC.stats))
  if (length(missed.cols)) {
    stop("Missing required columns in logFC.stats: ",
         paste(missed.cols, collasep=','))
  }
  lfc <- as.data.table(logFC.stats)
  lfc[, fid := rownames(logFC.stats)]

  if (!is.numeric(lfc$logFC) || !is.numeric(lfc$t)) {
    stop("logFC and t must be numeric columns")
  }
  if (nrow(lfc) != nrow(x)) {
    stop("nrow(logFC.stats) != nrow(x)")
  }

  scores <- lapply(gs.table@table$membership, function(idx) {
    data.table(mean.logFC=mean(lfc$logFC[idx], na.rm=TRUE),
               tmean.logFC=mean(lfc$logFC[idx], na.rm=TRUE, trim=trim),
               mean.t=mean(lfc$t[idx], na.rm=TRUE),
               tmean.t=mean(lfc$t[idx], na.rm=TRUE, trim=trim))
  })
  scores <- rbindlist(scores)

  cbind(gs.table@table[, list(group, id)], scores)
}

