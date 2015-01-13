setMethod("plot", 'MultiGSEAResult',
function(x, y, j, value=c('logFC', 't'), type=c('density', 'barcode', 'mini'),
         main=NULL, ...) {
  value <- match.arg(value)
  type <- match.arg(type)
  if (!is.character(y) && !is.character(j)) {
    stop("y,j should be characters to identifiy a geneset by collection,name")
  }
  if (missing(main)) {
    main <- paste(y, j, sep=':')
  }

  fn <- getFunction(sprintf('plot.gsea.%s', type))
  fn(x, y, j, value, main, ...)
})

plot.gsea.density <- function(x, y, j, value, main, bg.density=NULL,
                              xlim=base::range(c(bg.density$x, gs.density$x)),
                              ylim=c(0, base::max(c(gs.density$y, bg.density$y)) + 0.1),
                              jitter.sig.higher=FALSE, ...) {
  lfc <- logFC(x)
  if (!is(bg.density, 'density')) {
    bg.density <- density(lfc[[value]], na.rm=TRUE)
  }

  gs.fids <- featureIds(x@gsd, y, j, 'x.id')
  gs.lfc <- lfc[gs.fids]
  gs.lfc  <- gs.lfc[!is.na(featureId)]
  gs.n <- nrow(gs.lfc)
  gs.lfc$xval <- gs.lfc[[value]] ##jitter(gs.lfc[[value]])

  if (jitter.sig.higher) {
    gs.lfc$yval <- ifelse(gs.lfc$significant, rep(0.05, gs.n), rep(0.02, gs.n))
  } else {
    gs.lfc$yval <- rep(0.02, gs.n)
  }
  gs.lfc[, yval := pmax(jitter(yval, amount=0.02), 0)]

  if (gs.n >= 3) {
    gs.density <- density(gs.lfc[[value]])
  } else {
    gs.density <- list(x=numeric(), y=numeric())
  }

  plot(bg.density, main=main, sub=sprintf('(%d features)', gs.n), xlab=value,
       ylim=ylim, xlim=xlim, lwd=2)
  if (is(gs.density, 'density')) {
    lines(gs.density, col='red', lwd=3)
  }
  with(subset(gs.lfc, !significant), {
    points(xval, yval, pch=16, col='#FF0000AA', cex=0.5)
  })
  with(subset(gs.lfc, significant), {
    points(xval, yval, pch=16, col='#FF0000AA', cex=0.8)
    points(xval, yval, pch=1L, col='#000000AA', cex=0.8, lwd=2)
  })
  legend('topright', legend=c("all", "geneset"),
         text.col=c('black', 'red'))

  attr(gs.lfc, 'plot.attrs') <- list(xlim=xlim, ylim=ylim)
  invisible(gs.lfc)
}

plot.gsea.barcode <- function(x, y, j, value, main, ...) {
  lfc <- logFC(x)

  gs.fids <- featureIds(x@gsd, y, j, 'x.id')
  idx <- lfc[gs.fids, which=TRUE]
  idx <- idx[!is.na(idx)]
  gs.n <- length(idx)

  main <- sprintf("%s (%s values)", main, value)
  plot.barcode(lfc[[value]], idx, weights.label=value, main=main)
}

plot.gsea.mini <- function(x, y, j, value, main, ...) {
  stop("plot.gsea.mini not yet implemented")
}


##' limma's barcodeplot of one or two genesets (with enrichment)
##'
##'
##' This was imported from the barcode plot in limma_3.22.1 (the function was
##' annotated with a last update date of 20 October 2008)
plot.barcode <- function(statistics, index=NULL, index2=NULL, gene.weights=NULL,
                         weights.label="Weight", labels=c("Up", "Down"),
                         quantiles=c(-1, 1), col.bars=NULL, worm=TRUE,
                         span.worm=0.45, ...) {
  ##	Check statistics
  if (!is.vector(statistics, mode="numeric")) {
    stop("statistics should be a numeric vector")
  }
  nstat <- length(statistics)

  ##	Check index
  if (is.null(index)) {
    if (is.null(gene.weights)) {
      stop("Must specify at least one of index or gene.weights")
    } else {
      if (length(gene.weights) == nstat) {
        index <- rep_len(TRUE, nstat)
        index2 <- NULL
      } else {
        stop("No index and length(gene.weights) doesn't equal length(statistics)")
      }
    }
  } else {
    if (any(is.na(index)))
      stop("Need to provide index without NAs")
    if (is.logical(index))
      if(length(index) != nstat)
        stop("Length of index disagrees with statistics")
    if (length(index) > nstat)
      stop("Length of index disagrees with statistics")
  }

  ##	Check index2
  if (!is.null(index2)) {
    if (!is.null(gene.weights))
      warning("gene.weights ignored")
    gene.weights <- statistics
    gene.weights[] <- 0
    gene.weights[index] <- 1
    gene.weights[index2] <- -1
    index <- rep_len(TRUE, nstat)
    index2 <- NULL
  }

  ##	Check gene.weights
  if (!is.null(gene.weights)){
    if (!is.vector(gene.weights, mode = "numeric"))
      stop("gene.weights should be a numeric vector")
    if (any(is.na(gene.weights)))
      stop("Need to provide gene.weights without NAs")
    if (all(gene.weights == 0))
      stop("gene.weights equal to zero: no selected genes to plot")
    if (length(gene.weights) != length(statistics[index]))
      stop("Length of gene.weights disagrees with size of set")

    one <- all(gene.weights >= 0) | all(gene.weights <= 0)

    if (one) {
      index2 <- NULL
      gene.weights1 <- rep_len(0, nstat)
      names(gene.weights1) <- names(statistics)
      gene.weights1[index] <- gene.weights

      index <- rep_len(FALSE, nstat)
      names(index) <- names(statistics)
      index[gene.weights1 != 0] <- TRUE

      gene.weights1 <- gene.weights1[gene.weights1 != 0]
      gene.weights <- gene.weights1
    } else {
      gene.weights12 <- rep_len(0, nstat)
      names(gene.weights12) <- names(statistics)
      gene.weights12[index] <- gene.weights

      index <- index2 <- rep_len(FALSE, nstat)
      names(index) <- names(index2) <- names(statistics)
      index[gene.weights12 > 0] <- TRUE
      index2[gene.weights12 < 0] <- TRUE

      gene.weights1 <- gene.weights12[gene.weights12 > 0]
      gene.weights2 <- gene.weights12[gene.weights12 < 0]

      gene.weights <- gene.weights1
    }
  }

  ##	Are there up and down sets?
  TWO <- !is.null(index2)

  ##	Convert indexes to logical and add gene.weights
  set1 <- set2 <- data.frame(idx = rep.int(FALSE,nstat), weight = NA, wt = NA)
  rownames(set1) <- rownames(set2) <- names(statistics)
  set1$idx[index] <- TRUE
  if (TWO) set2$idx[index2] <- TRUE

  if (length(gene.weights)) {
    set1$weight <- 0
    set1$weight[index]<- gene.weights
    set1$wt <- abs(set1$weight)/sum(abs(set1$weight))
    if (TWO) {
      set2$weight <- 0
      set2$weight[index2] <- gene.weights2

      SUM <- sum(abs(set1$weight), abs(set2$weight))
      set1$wt <- abs(set1$weight)/SUM
      set2$wt <- abs(set2$weight)/SUM

    }
  }

  ## TODO: enable sorting of statistics from low to high
  ## Sort statistics and indexes
  ostat <- order(statistics, na.last = TRUE, decreasing=TRUE)
  statistics <- statistics[ostat]
  set1 <- set1[ostat,]

  if (TWO) set2 <- set2[ostat,]

  ##	Trim off missing values
  n <- sum(!is.na(statistics))
  if (n==0L) {
    message("No valid statistics")
    return(invisible())
  }
  if (n < nstat) {
    statistics <- statistics[1:n]
    set1 <- set1[1:n,]
    if (TWO) set2 <- set2[1:n,]
  }

  ##	Convert indexes to integer
  r <- which(set1$idx)

  if (TWO) {
    r2 <- which(set2$idx)
    if(!length(r2)) TWO <- FALSE
  }

  #	Check there is something to plot
  if (!length(r)) {
    if(TWO) {
      r <- r2
      set1 <- set2
      TWO <- FALSE

    } else {
      message("No selected genes to plot")
      return(invisible())
    }
  }

  ##	Are there unequal weights?
  WTS <- FALSE
  wt1 <- set1$wt[r]
  len.up <- 1

  if (all(!is.na(wt1))) {
    len.up <- set1$weight[r]/max(abs(set1$weight[r]))

    anydifferent <- function(x) {
      if(length(x) < 2) return(FALSE)
      r <- range(x)
      (r[2] > r[1])
    }

    if (!TWO) if (anydifferent(wt1)) WTS <- TRUE

    if (TWO) {
      wt12 <- c(set1$wt[r], abs(set2$wt[r2]))
      if(anydifferent(wt12)) WTS <- TRUE

      max.wt <- max(set1$wt[r], set2$wt[r2])
      len.up <- set1$wt[r]/max.wt
      len.down <- set2$wt[r2]/max.wt
    }
  }

  pos.dir <- all(len.up > 0)

  #	Plot setting
  if (WTS) shift <- 0.1 else shift <- 0

  #	Check other arguments
  quantiles <- sort(quantiles)

  #	color settting
  ALP <- 0.4

  if (is.null(col.bars)) {
    if (TWO) {
      col.bars <- c("red", "blue")
      if(WTS) col.bars.alpha <- c(rgb(1,0,0,alpha=ALP), rgb(0,0,1,alpha=ALP))
      else col.bars.alpha <- col.bars

    } else {
      col.bars <- "black"
      if(WTS) col.bars.alpha <- rgb(0,0,0,alpha=ALP)
      else col.bars.alpha <- col.bars
    }
  } else {
    if (TWO) {
      if(length(col.bars) == 1) col.bars <- rep(col.bars, 2)
      RGB <- col2rgb(col.bars)/255
      red <- RGB[1,1]
      green <- RGB[2,1]
      blue <- RGB[3,1]
      red2 <- RGB[1,2]
      green2 <- RGB[2,2]
      blue2 <- RGB[3,2]
      if (WTS) {
        col.bars.alpha <- c(rgb(red, green, blue, alpha=ALP),
                            rgb(red2, green2, blue2, alpha=ALP))
      } else {
        col.bars.alpha <- col.bars
      }
    } else {
      RGB <- col2rgb(col.bars)/255
      red <- RGB[1,1]
      green <- RGB[2,1]
      blue <- RGB[3,1]
      if(WTS) col.bars.alpha <- rgb(red, green, blue, alpha=ALP)
      else col.bars.alpha <- col.bars
    }
  }

  ## TODO: Enable different labeling of of statistics
  ## worm plot setting
  ylim.worm <- ylim <- c(-1, 1)
  ylab.worm <- ""
  xlab.worm <- "statistics"

  if (!TWO) ylim.worm <- c(0, 1)

  if (worm) {
    ylim.worm <- c(-2.1, 2.1)
    if (!TWO) ylim.worm <- c(0, 2.1)
  }

  ylim[2] <- ylim[2] + 0.5
  if (TWO) ylim[1] <- ylim[1] - 0.5

  if (TWO) {
    plot(1:n, xlim=c(0,n), ylim=c(ylim.worm[1]-shift, ylim.worm[2]+shift),
         type="n", axes=FALSE, xlab=xlab.worm, ylab=ylab.worm, ...)
  }
  if (!TWO) {
    plot(1:n, xlim=c(0,n),
         ylim=c(ylim.worm[1]-shift*(!pos.dir), ylim.worm[2]+shift*pos.dir),
         type="n", axes=FALSE, xlab=xlab.worm,ylab=ylab.worm, ...)
  }

  npos <- sum(statistics > quantiles[2])
  nneg <- sum(statistics < quantiles[1])

  lwd <- 50/length(r)
  lwd <- min(1.9, lwd)
  lwd <- max(0.2, lwd)

  if (TWO){
    lwd2 <- 50/length(r2)
    lwd2 <- min(1.9, lwd2)
    lwd2 <- max(0.2, lwd2)

    lwd <- lwd2 <- min(lwd, lwd2)
  }

  barlim <- ylim[2] - c(1.5, 0.5)

  if (!pos.dir) {
    rect.yb <- 0.5
    rect.yt <- 1
    rect(npos + 0.5, rect.yb, n - nneg + 0.5, rect.yt,
         col = "lightgray", border = NA)
    if (npos)
      rect(0.5, rect.yb, npos + 0.5, rect.yt, col = "pink",
           border = NA)
    if (nneg)
      rect(n - nneg + 0.5, rect.yb, n + 0.5, rect.yt, col = "lightblue",
           border = NA)

    segments(r, barlim[2]/2, r, barlim[2], lwd = lwd, col = col.bars.alpha[1])
    segments(r, barlim[2]/2-shift, r, barlim[2]/2*(1+len.up)-shift, lwd = lwd,
             col = col.bars[1])
  }

  if (pos.dir) {
    rect.yb <- -0.5
    if (!TWO) rect.yb <- 0
    rect.yt <- 0.5

    rect(npos + 0.5, rect.yb, n - nneg + 0.5, rect.yt, col = "lightgray",
         border = NA)
    if (npos)
      rect(0.5, rect.yb, npos + 0.5, rect.yt, col = "pink", border = NA)
    if (nneg)
      rect(n - nneg + 0.5, rect.yb, n + 0.5, rect.yt, col = "lightblue",
           border = NA)

    segments(r, barlim[1], r, barlim[2]/2, lwd = lwd, col = col.bars.alpha[1])
    segments(r, barlim[2]/2+shift, r, barlim[2]/2*(1+len.up)+shift, lwd = lwd,
             col = col.bars[1])
  }

  if (TWO) {
    barlim2 <- ylim[1] + c(0.5, 1.5)
    segments(r2, barlim2[1]/2, r2, barlim2[2], lwd = lwd2,
             col = col.bars.alpha[2])
    segments(r2, barlim2[1]/2*(1+len.down)-shift, r2, barlim2[1]/2-shift,
             lwd = lwd2, col = col.bars[2])
  }

  lab.at <- 0
  if (!TWO) lab.at <- 0.5
  axis(side = 2, at = lab.at, padj = 3.8, cex.axis = 0.85, labels = labels[1],
       tick = FALSE)
  axis(side = 4, at = lab.at, padj = -3.8, cex.axis = 0.85, labels = labels[2],
       tick = FALSE)

  ## label statistics on x axis
  prob <- (10:0)/10
  axis(at = seq(1,n,len=11), side = 1, cex.axis = 0.7, las = 2,
       labels = format(quantile(statistics, p = prob), digits = 1))

  ## create worm
  if (worm) {
    # rescale x to new range
    rescale <- function(x, newrange, oldrange=range(x)) {
      newrange[1] + (x-oldrange[1]) / (oldrange[2]-oldrange[1]) * (newrange[2] - newrange[1])
    }

    # calculate enrichment
    if (!WTS) {
      ave.enrich1 <- length(r)/n
      worm1 <- tricubeMovingAverage(set1$idx, span = span.worm)/ave.enrich1
      if(TWO) {
        ave.enrich2 <- length(r2)/n
        worm2 <- tricubeMovingAverage(-set2$idx, span = span.worm)/ave.enrich2
      }
    }

    ## calculate weighted enrichment
    if (WTS) {
      ave.enrich1 <- mean(set1$wt)
      worm1 <- tricubeMovingAverage(set1$wt, span = span.worm)/ave.enrich1

      if (TWO) {
        ave.enrich2 <- mean(set2$wt)
        worm2 <- tricubeMovingAverage(-set2$wt, span = span.worm)/ave.enrich2
      }
    }

    # rescale worm
    max.worm1 <- max(worm1)
    r.worm1 <- c(0,max.worm1)
    worm1.scale <- rescale(worm1,
                           newrange=c(1.1+shift*pos.dir,2.1+shift*pos.dir),
                           oldrange=r.worm1)

    if (TWO) {
      min.worm2 <- min(worm2)
      r.worm2 <- c(min.worm2,0)
      worm2.scale <- rescale(worm2,
                             newrange=c(-2.1-shift,-1.1-shift),
                             oldrange=r.worm2)
    }

    ## plot worms
    if (!TWO) {
      lines(x = 1:n, y = worm1.scale, col = col.bars[1], lwd = 2)
      abline(h = rescale(1,newrange=c(1.1+shift*pos.dir,2.1+shift*pos.dir),oldrange=r.worm1),
             lty=2)
      axis(side = 2,
           at = c(1.1+shift*pos.dir, 2.1+shift*pos.dir),
           cex.axis = 0.8,
           labels = c(0, format(max.worm1, digits = 2)))
      axis(side = 2,
           labels = "Enrichment",
           at = 1.6+shift*pos.dir, padj = -0.6,
           tick = FALSE,
           cex.axis = 0.8)
    }

    if (TWO) {
      lines(x = 1:n, y = worm1.scale, col = col.bars[1], lwd = 2)
      abline(h = rescale(1,newrange=c(1.1+shift,2.1+shift),oldrange=r.worm1),
             lty=2)
      lines(x = 1:n, y = worm2.scale, col = col.bars[2], lwd = 2)
      abline(h = rescale(-1,newrange=c(-2.1-shift,-1.1-shift),oldrange=r.worm2),
             lty=2)

      axis(side = 2,
           at = c(1.1+shift, 2.1+shift),
           cex.axis = 0.7,
           labels = c(0, format(max.worm1, digits = 2)))
      axis(side = 2,
           at = c(-1.1-shift, -2.1-shift),
           cex.axis = 0.7,
           labels = c(0, format(-min.worm2, digits = 2)))
      axis(side = 2,
           labels = "Enrichment",
           at = 1.6+shift,
           tick = FALSE,
           padj = -0.6,
           cex.axis = 0.7)
      axis(side = 2,
           labels = "Enrichment",
           at = -1.6-shift,
           tick = FALSE,
           padj = -0.6,
           cex.axis = 0.7)
    }
  }

  ## add gene.weights axis
  if (WTS){
    if (!TWO){
      if (pos.dir) {
        axis(side = 2,
             at = c(0.5+shift, 1+shift),
             cex.axis = 0.48,
             padj = 1.6,
             labels = c(0, format(max(set1$weight[r]), digits = 2)))
        axis(side = 2,
             labels = weights.label[1],
             at = 0.75+shift,
             padj = 1,
             tick = FALSE,
             cex.axis = 0.5)
      }

      if (!pos.dir){
        axis(side = 2,
             at = c(0-shift, 0.5-shift),
             cex.axis = 0.48,
             padj = 1.6,
             labels = c(format(min(set1$weight[r]), digits = 2), 0))
        axis(side = 2,
             labels = weights.label[1],
             at = 0.25-shift,
             padj = 1,
             tick = FALSE,
             cex.axis = 0.5)
      }
    }

    if (TWO){
      max.weight <- max(set1$weight[r], abs(set2$weight[r2]))
      axis(side = 2,
           at = c(0.5+shift, 1+shift),
           cex.axis = 0.43,
           labels = c(0, format(max.weight, digits = 2, scientific = FALSE)),
           padj = 1.6)
      axis(side = 2,
           labels = weights.label[1],
           at = 0.75+shift,
           padj = 1,
           tick = FALSE,
           cex.axis = 0.46)
      axis(side = 2,
           at = c(-0.5-shift, -1-shift),
           cex.axis = 0.43,
           labels = c(0, format(-max.weight, digits = 2, scientific = FALSE)),
           padj = 1.6)
      axis(side = 2,
           labels = weights.label[1],
           at = -0.75-shift,
           padj = 1,
           tick = FALSE,
           cex.axis = 0.46)
    }
  }

  invisible()
}

## -----------------------------------------------------------------------------
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
##' @return A data.table that has the collection, id, gene set name, image page
##'   and whether or not the image was plotted columns for each row in \code{x}
generate.GSEA.plots <- function(x, design, contrast, result, outdir,
                                use.cache=TRUE, logFC.stats=NULL,
                                padj.threshold=1,
                                score.by=c('logFC', 't'), ...) {
  score.by <- match.arg(score.by)
  score.by <- 'logFC'
  if (is.null(logFC.stats)) {
    logFC.stats <- calculateIndividualLogFC(x, design, contrast,
                                            provide='table', ...)
  }
  if (!all(rownames(logFC.stats) == rownames(x))) {
    stop("The logFC values do not match the rownames if x")
  }

  x.axis <- switch(score.by, logFC='logFC', t='t-statistic')
  d.summary <- sprintf('%s(%s)', score.by, design.params.name(design, contrast))
  fn.base <- paste0('%s-', d.summary, '.png')

  scores <- setNames(logFC.stats[[score.by]], rownames(logFC.stats))

  bg.dens <- density(scores, na.rm=TRUE)
  xrange <- c(min(scores, na.rm=TRUE) - 0.5, max(scores, na.rm=TRUE) + 0.5)
  f.ids <- names(scores)

  should.plot <- significantGeneSets(result, 'adjusted', padj.threshold)
  work.me <- result[, list(collection, id)]
  work.me[, gs.name := paste(collection, id, sep='.')]
  work.me[, fn := file.path(outdir, 'images', sprintf(fn.base, gs.name))]
  work.me[, fn.exists.pre := file.exists(fn)]
  work.me[, do.plot := should.plot & (!fn.exists.pre | !use.cache)]

  for (i in which(work.me$do.plot)) {
    gs.features <- result$feature.id[[i]]
    gs.logFC <- scores[gs.features]
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
