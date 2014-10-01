##' Generates GSEA plots for the genesets over a given experiment.
##'
##' @param x The expression object used in the multiGSEA call
##' @param design The desing matrix for this GSEA
##' @param contrast The contrast being analyzed in the GSEA
##' @param result The result data.table returned (or being generated) in the
##' \code{multiGSEA} call.
##' @param outdir The path to the parent directory that holds the output for
##' this run.
##' @param use.cache If true, the function checks if the image already exists
##' before replotting it.
##' @param logFC pre-computed logFCs for the contrast we are running GSEA
##' against.
##' @param padj.threshold plot generation to only include genesets that
##' have some value in an FDR column that is less than or equal to this
##' threshold.
##'
##' @return A data.table that has the group, id, and path columns. path
##' is the path to the image for the given geneset.
generate.GSEA.plots <- function(x, design, contrast, result, outdir,
                                use.cache=TRUE, logFC=NULL,
                                padj.threshold=1) {
  if (is.null(logFC)) {
    logFC <- calculateIndividualLogFC(x, design, contrast, ...)
  }
  if (!all(names(logFC) == rownames(x))) {
    stop("The logFC values do not match the rownames if x")
  }

  d.summary <- sprintf('logFC(%s)', design.params.name(design, contrast))
  fn.base <- paste0('%s-', d.summary, '.png')

  if (padj.threshold < 1) {
    padj.cols <- grep('^padj\\.', colnames(result))
    padj.m <- as.matrix(result[, padj.cols, with=FALSE])
    keep <- rowSums(padj.m <= padj.threshold) > 0
    if (!any(keep)) {
      warning("No plots generated due to FDR cutoff", immediate.=TRUE)
      return(data.table(group=character(), id=character(), gs.name=character(),
                        fn=character()))
    }
    result <- result[keep,]
  }

  if (nrow(result) == 0) {
    return(NULL)
  }

  bg.dens <- density(logFC, na.rm=TRUE)
  xrange <- c(min(logFC, na.rm=TRUE) - 0.5, max(logFC, na.rm=TRUE) + 0.5)
  f.ids <- names(logFC)

  fns <- lapply(1:nrow(result), function(i) {
    gs.group <- result$group[i]
    gs.id <- result$id[i]
    gs.name <- paste(gs.group, gs.id, sep='.')
    fn.out <- file.path(outdir, 'images', sprintf(fn.base, gs.name))

    if (!file.exists(fn.out) || !use.cache) {
      gs.features <- result$feature.id[[i]]
      gs.logFC <- logFC[gs.features]
      gs.logFC <- gs.logFC[!is.na(gs.logFC)]
      if (length(gs.logFC) < 3) {
        gs.dens <- NULL
        ymax <- base::max(bg.dens$y)
      } else {
        gs.dens <- density(gs.logFC, na.rm=TRUE)
        ymax <- base::max(bg.dens$y, gs.dens$y)
      }

      xstats <- sprintf('(%d genes)', length(gs.logFC))

      png(fn.out, 800, 800, res=150)
      plot(bg.dens, lwd=2, xlab="logFC", main=gs.name, sub=xstats,
           ylim=c(0, ymax), xlim=xrange)
      if (!is.null(gs.dens)) {
        lines(gs.dens, col='red', lwd=3)
      }
      rug(gs.logFC, col='#FF000033', lwd=3)
      legend('topright', legend=c("all logFC", "geneset logFC"),
             text.col=c('black', 'red'))
      dev.off()
    }
    data.table(group=gs.group, id=gs.id, gs.name=gs.name, fn=fn.out)
  })
  all.fns <- rbindlist(fns)
  all.fns
}
