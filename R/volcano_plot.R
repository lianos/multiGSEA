.volcano.sources <- c(
  'MultiGSEAResultContainer',
  'MultiGSEAResult',
  'data.frame')

##' Create an interactive volcano plot
##'
##' @export
##'
##' @param highlight_genes A vector of featureIds to highlight, or a GeneSetDb
##'   that we can extract the featureIds from for this purpose.
volcano_plot <- function(x, stats='dge', xaxis='logFC', yaxis='pval', idx,
                         highlight_genes=NULL,
                         xtfrm=identity, ytfrm=function(vals) -log10(vals),
                         horiz_lines=c('padj'=0.10),
                         xgran=NULL, ygran=NULL) {
  if (FALSE) {
    x <- readRDS('~/tmp/schmidt/multiGSEA-EP-uber_hWT_tKO-hWT_tWT.rds')
    stats='dge'; xaxis='logFC'; yaxis='pval'; idx='idx';
    xtfrm=identity; ytfrm=function(vals) -log10(vals)
    xgran=NULL; ygran=NULL;
    # xgran=1.5; ygran=0.05;
    horiz_lines=c('padj'=0.10);
    highlight_genes=x;
  }
  ## NOTE: I should use S3 or S4 here, but I'm lazy right now.
  dat <- volcano.stats.table(x, stats, xaxis, yaxis, idx, xtfrm, ytfrm)

  yvals <- dat[['yaxis']]
  yrange <- range(yvals)
  ypad <- 0.05 * diff(yrange)
  ylim <- c(yrange[1] - ypad, yrange[2] + ypad)

  xvals <- dat[['xaxis']]
  xrange <- range(xvals)
  xpad <- 0.05 * diff(xrange)
  xlim <- c(xrange[1] - xpad, xrange[2] + xpad)

  ## Check if we want to hexbin some of the data: values less than
  ## abs(input$xgran) on the xaxis or less than input$ygran on the y axis
  ## should be hexbinized
  do.hex <- is.numeric(xgran) && is.numeric(ygran)
  if (do.hex) {
    xthresh <- xtfrm(xgran)
    ythresh <- ytfrm(ygran)
    hex.me <- abs(dat[['xaxis']]) <= xthresh & dat[['yaxis']] <= ythresh
    hex <- dat[hex.me,,drop=FALSE]
    pts <- dat[!hex.me,,drop=FALSE]
  } else {
    pts <- dat
    hex <- dat[FALSE,]
  }

  ## Build the figure
  ## Setup the initial figure
  tools <- c('box_select', 'box_zoom', 'reset', 'save')
  p <- figure(xlab=xaxis, ylab=sprintf('-log10(%s)', yaxis),
              xlim=xlim, ylim=ylim, webgl=TRUE, tools=tools)
  if (do.hex) {
    p <- ly_hexbin(p, 'xaxis', 'yaxis', data=hex, xbins=30)
  }
  if (nrow(pts) > 0) {
    col <- rep("#404040", nrow(pts))
    if (!is.null(highlight_genes)) {
      fids <- extract.genes(highlight_genes)
      col[pts$featureId %in% fids] <-"#EE4000"
    }

    p <- ly_points(p, 'xaxis', 'yaxis', data=pts, lname='points', size=4,
                   color=col, legend=FALSE,
                   hover=list(symbol,
                              logFC=sprintf('%.3f', xaxis),
                              pval=sprintf('%.3f', pval),
                              qval=sprintf('%.3f', padj)))
  }

  ## Add horizontal lines to indicate where padj of 0.10 lands
  if (is.numeric(horiz_lines)) {
    # q.thresh <- c(0.05, 0.10, 0.20)
    q.thresh <- c(horiz_lines)
    horiz.unit <- names(horiz_lines)[1L]
    if (yaxis == horiz.unit) {
      ypos <- horiz_lines
    } else if (yaxis == 'pval' && horiz.unit == 'padj') {
      ## Find the pvals that are closest to qvals without going over
      ypos <- sapply(q.thresh, function(val) {
        pdiffs <- abs(dat$padj - val)
        idx <- which.min(pdiffs)
        idiff <- pdiffs[idx]
        if (idiff < 1e-2) {
          # message("Closest padj: ", dat$padj[idx])
          return(dat$pval[idx])
        } else {
          # message("No pval with qval close to: ", sprintf('%.3f', val))
          return(NA)
        }
      })
    } else {
      warning("Not drawing horiz_lines in volcano", immediate.=TRUE)
      ypos <- rep(NA_real_, length(q.thresh))
    }
    names(ypos) <- NULL
    ypos <- ypos[!is.na(ypos)]
    if (length(ypos) > 0) {
      ypos <- ytfrm(ypos)
      ablabel <- if (horiz.unit == 'padj') 'q-value' else horiz.unit
      p <- ly_abline(p, h=ypos, color='red', type=2)
      p <- ly_text(p, xrange[1], ypos[1],
                   sprintf('%s: %.2f', ablabel, horiz_lines[1L]),
                   color='red', font_size='9pt')
    }
  }

  p$data <- list(pts=pts, hex=hex)
  p
}

##' @export
extract.genes <- function(x, ...) {
  stopifnot(is.character(x) || is(x, 'GeneSetDb') || is(x, 'MultiGSEAResult'))
  if (is.character(x)) {
    return(x)
  }
  featureIds(x)
}

##' @export
volcano.source.type <- function(x) {
  is.valid <- sapply(.volcano.sources, function(src) is(x, src))
  if (sum(is.valid) != 1L) {
    stop("Illegal object type for source of volcano: ", class(x)[1L])
  }
  .volcano.sources[is.valid]
}

##' @export
volcano.stats.table <- function(x, stats='dge', xaxis='logFC', yaxis='padj',
                                idx='idx',
                                xtfrm=identity,
                                ytfrm=function(vals) -log10(vals)) {
  type <- volcano.source.type(x)
  stopifnot(is.function(xtfrm), is.function(ytfrm))
  if (is(x, 'MultiGSEAResultContainer')) {
    x <- x()$mg
  }
  if (is(x, 'MultiGSEAResult')) {
    stats <- match.arg(stats, c('dge', resultNames(x)))
    if (stats == 'dge') {
      x <- logFC(x)
      idx <- 'featureId'
    } else {
      x <- result(x, stats)
      idx <- 'idx'
      x[[idx]] <- paste(x$collection, x$name, sep=';;')
    }
  }

  if (!is.data.frame(x)) {
    stop('x should have been converted into a data.frame by now')
  }
  x <- as.data.frame(x)
  missed.cols <- setdiff(c(xaxis, yaxis, idx), names(x))
  if (length(missed.cols)) {
    stop("Missing columns from stats data.frame from volcano object:\n  ",
         paste(missed.cols, collapse=','))
  }
  x[['xaxis']] <- xtfrm(x[[xaxis]])
  x[['yaxis']] <- ytfrm(x[[yaxis]])
  x[['idx']] <- x[[idx]]
  x
}
