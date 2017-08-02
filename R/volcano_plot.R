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
                         xtfrm=base::identity,
                         ytfrm=function(vals) -log10(vals),
                         xlab=xaxis, ylab=sprintf('-log10(%s)', yaxis),
                         highlight_genes=NULL,
                         horiz_lines=c('padj'=0.10),
                         xhex=NULL, yhex=NULL,
                         point.size=5,
                         tools=c('box_select', 'reset', 'save'),
                         width=NULL, height=NULL,
                         shiny_source='mgvolcano', ggtheme=theme_bw(), ...) {
  ## NOTE: I should use S3 or S4 here, but I'm lazy right now.
  dat <- volcano.stats.table(x, stats, xaxis, yaxis, idx, xtfrm, ytfrm)

  yvals <- dat[['.yvt']]
  yrange <- range(yvals)
  ypad <- 0.05 * diff(yrange)
  ylim <- c(yrange[1] - ypad, yrange[2] + ypad)

  xvals <- dat[['.xvt']]
  xrange <- range(xvals)
  xpad <- 0.05 * diff(xrange)
  xlim <- c(xrange[1] - xpad, xrange[2] + xpad)

  ## Check if we want to hexbin some of the data: values less than
  ## abs(input$xhex) on the xaxis or less than input$yhex on the y axis
  ## should be hexbinized
  do.hex <- is.numeric(xhex) && is.numeric(yhex)
  if (do.hex) {
    hex.me <- abs(dat[['.xv']]) <= xhex | dat[['.yvt']] <= yhex
  } else {
    hex.me <- rep(FALSE, nrow(dat))
  }
  hex <- dat[hex.me,,drop=FALSE]
  pts <- dat[!hex.me,,drop=FALSE]

  # silence R CMD check NOTEs
  .xvt <- .yvt <- .xv <- .yv <- NULL
  gg <- ggplot(dat, aes(.xvt, .yvt))
  if (nrow(pts)) {
    if ('symbol' %in% names(dat)) {
      gg <- gg + suppressWarnings({
        geom_point(aes(key=featureId,
                       text=paste0('Symbol: ', symbol, '<br>',
                                   xaxis, ': ', sprintf('%.3f', .xv), '<br>',
                                   yaxis, ': ', sprintf('%.3f', .yv))),
                   data=pts)
      })
    } else {
      gg <- gg +
        geom_point(aes(key=featureId,
                       text=paste0('featureId: ', featureId, '<br>',
                                   xaxis, ': ', sprintf('%.3f', .xv), '<br>',
                                   yaxis, ': ', sprintf('%.3f', .yv))),
                   data=pts)
    }
  }
  if (nrow(hex)) {
    gg <- gg + geom_hex(data=hex, bins=50)
  }

  ## Add horizontal lines to indicate where padj of 0.10 lands
  lpos <- NULL
  if (is.numeric(horiz_lines) && !any(is.na(horiz_lines))) {
    # q.thresh <- c(0.05, 0.10, 0.20)
    q.thresh <- c(horiz_lines)
    horiz.unit <- names(horiz_lines)[1L]
    if (yaxis == horiz.unit) {
      ypos <- horiz_lines
    } else if (yaxis == 'pval' && horiz.unit == 'padj') {
      ## Find the pvals that are closest to qvals without going over
      ypos <- sapply(q.thresh, function(val) {
        approx.target.from.transformed(val, dat$pval, dat$padj)
      })
    } else {
      warning("Not drawing horiz_lines in volcano", immediate.=TRUE)
      ypos <- rep(NA_real_, length(q.thresh))
    }
    names(ypos) <- NULL
    ypos <- ypos[!is.na(ypos)]
    if (length(ypos) > 0) {
      lpos <- ytfrm(ypos)
      ablabel <- if (horiz.unit == 'padj') 'q-value' else horiz.unit
      lbl <- sprintf('%s: %.2f', ablabel, horiz_lines[1L])
      # txtdat <- data.frame(x=xrange[1], y=lpos[1], lbl=lbl)
      # gg <- gg +
      #   geom_hline(yintercept=lpos, color='red', linetype='dashed') +
      #   geom_text(aes(x, y, label=lbl), data=txtdat, vjust=1.5, hjust=0.2, color='red')
      #
      # gg <- gg +
      #   geom_hline(aes(text=lbl), yintercept=lpos, color='red', linetype='dashed')
    }
  }

  if (is(ggtheme, 'theme')) {
    gg <- gg + ggtheme
  }

  ## ggplotly messages to use github: hadley/ggplot2
  p <- suppressMessages(ggplotly(gg, width=width, height=height, tooltip='text'))
  if (is.numeric(lpos)) {
    p <- add_lines(p, x=seq(xrange[1] - 1, xrange[2] + 1, length=100),
                   y=rep(lpos, 100),
                   inherit=FALSE,
                   text=lbl, color=I("red"),
                   line=list(dash='dash', width=2),
                   hoverinfo='text')
  }
  p <- layout(p, dragmode="select", autosize=FALSE, xaxis=list(title=xlab),
              yaxis=list(title=ylab))
  p <- plotly_build(p)

  ## Disable hovering on hexbin
  p$x$source <- shiny_source
  if (nrow(hex)) {
    hexidx <- 2:(length(p$x$data) -1L)
    for (idx in hexidx) {
      p$x$data[[idx]]$hoverinfo <- 'none'
      # len <- length(p$x$data[[idx]]$x)
      # p$x$data[[idx]]$text <- sprintf('count: %d', rep(len))
    }
  }

  config(p, collaborate=FALSE, displaylogo=FALSE)
}

##' Get the approximate nominal pvalue for a target qvalue given the
##' distribution of adjusted pvalues from the nominal ones.
##'
##' @param pval the adjust pvalue you want the nominal value for
##' @param pvals the distribution of nominal pvalues
##' @param padjs the adjusted pvalues from \code{pvals}
##' @param thresh how close padj has to be in padjs to get its nominal
##'   counterpart
##' @return numeric answer, or NA if can't find nominal pvalue within given
##'   threshold for padjs
approx.target.from.transformed <- function(target, orig, xformed,
                                           thresh=1e-2) {
  xdiffs <- abs(xformed - target)
  idx <- which.min(xdiffs)
  idiff <- xdiffs[idx]
  if (idiff < thresh) {
    # message("Closest padj: ", dat$padj[idx])
    out <- orig[idx]
  } else {
    # message("No pval with qval close to: ", sprintf('%.3f', val))
    out <- NA_real_
  }
  out
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
volcano.stats.table <- function(x, stats='dge', xaxis='logFC', yaxis='pval',
                                idx='idx',
                                xtfrm=identity,
                                ytfrm=function(vals) -log10(vals)) {
  stopifnot(is.function(xtfrm), is.function(ytfrm))
  type <- volcano.source.type(x)
  if (is(x, 'MultiGSEAResultContainer')) {
    x <- x$mg
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
  if (!'featureId' %in% names(x)) {
    ids <- rownames(x)
    if (is.null(ids)) {
      ids <- as.character(1:nrow(x))
    }
    x[['featureId']] <- ids
  }
  x[['.xv']] <- x[[xaxis]]
  x[['.yv']] <- x[[yaxis]]
  x[['.xvt']] <- xtfrm(x[[xaxis]])
  x[['.yvt']] <- ytfrm(x[[yaxis]])
  x
}
