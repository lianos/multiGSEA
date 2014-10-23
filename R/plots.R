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
##' @return A data.table that has the group, id, gene set name, image page
##'   and whether or not the image was plotted columns for each row in \code{x}
generate.GSEA.plots <- function(x, design, contrast, result, outdir,
                                use.cache=TRUE, logFC=NULL,
                                padj.threshold=1,
                                score.by=c('logFC', 't'), ...) {
  score.by <- match.arg(score.by)
  if (is.null(logFC)) {
    logFC <- calculateIndividualLogFC(x, design, contrast, provide=score.by,
                                      ...)
  }
  if (!all(names(logFC) == rownames(x))) {
    stop("The logFC values do not match the rownames if x")
  }

  x.axis <- switch(score.by, logFC='logFC', t='t-statistic')
  d.summary <- sprintf('%s(%s)', score.by, design.params.name(design, contrast))
  fn.base <- paste0('%s-', d.summary, '.png')

  bg.dens <- density(logFC, na.rm=TRUE)
  xrange <- c(min(logFC, na.rm=TRUE) - 0.5, max(logFC, na.rm=TRUE) + 0.5)
  f.ids <- names(logFC)

  should.plot <- significantGeneSets(result, 'adjusted', padj.threshold)
  work.me <- result[, list(group, id)]
  work.me[, gs.name := paste(group, id, sep='.')]
  work.me[, fn := file.path(outdir, 'images', sprintf(fn.base, gs.name))]
  work.me[, fn.exists.pre := file.exists(fn)]
  work.me[, do.plot := should.plot & (!fn.exists.pre | !use.cache)]

  for (i in which(work.me$do.plot)) {
    gs.features <- result$feature.id[[i]]
    gs.logFC <- logFC[gs.features]
    gs.logFC <- gs.logFC[!is.na(gs.logFC)]

    if (length(gs.logFC) < 3) {
      ## You can't plot a density of this
      gs.dens <- NULL
      ymax <- base::max(bg.dens$y)
    } else {
      gs.dens <- density(gs.logFC, na.rm=TRUE)
      ymax <- base::max(bg.dens$y, gs.dens$y)
    }

    xstats <- sprintf('(%d genes)', length(gs.logFC))

    png(work.me$fn[i], 800, 800, res=150)
    plot(bg.dens, main=work.me$gs.name[i], sub=xstats, xlab=x.axis,
         ylim=c(0, ymax), xlim=xrange, lwd=2, )
    if (!is.null(gs.dens)) {
      lines(gs.dens, col='red', lwd=3)
    }
    rug(gs.logFC, col='#FF000033', lwd=3)
    legend('topright', legend=c("all", "geneset"),
           text.col=c('black', 'red'))
    dev.off()
  }

  work.me[, fn.exists.post := file.exists(fn)]
  work.me
}
