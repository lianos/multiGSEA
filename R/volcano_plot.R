.volcano.sources <- c(
  'MultiGSEAResultContainer',
  'MultiGSEAResult',
  'data.frame')

##' Create an interactive volcano plot
##'
##' @export
##'
##' @param highlight A vector of featureIds to highlight, or a GeneSetDb
##'   that we can extract the featureIds from for this purpose.
volcano_plot <- function(x, stats='dge', xaxis='logFC', yaxis='pval', idx,
                         xtfrm=base::identity,
                         ytfrm=function(vals) -log10(vals),
                         xlab=xaxis, ylab=sprintf('-log10(%s)', yaxis),
                         highlight=NULL,
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
  if (is.character(highlight)) {
    hlite <- dat[dat[['featureId']] %in% highlight,,drop=FALSE]
  } else {
    hlite <- dat[FALSE,,drop=FALSE]
  }
  pts <- dat[!hex.me,,drop=FALSE]
  if (nrow(hlite)) {
    pts <- pts[!pts[['featureId']] %in% hlite[['featureId']],,drop=FALSE]
  }

  # silence R CMD check NOTEs
  .xvt <- .yvt <- .xv <- .yv <- NULL
  gg <- ggplot(dat, aes(.xvt, .yvt))
  if (nrow(hex)) {
    gg <- gg + geom_hex(data=hex, bins=50)
  }
  gg <- mg_add_points(gg, pts)
  gg <- mg_add_points(gg, hlite, color='red')

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
    }
  }

  if (is(ggtheme, 'theme')) {
    gg <- gg + ggtheme
  }

  ## ggplotly messages to use github: hadley/ggplot2
  # p <- suppressMessages(ggplotly(gg, width=width, height=height, tooltip='text'))
  p <- suppressMessages({
    ggplotly(gg, width=width, height=height, source=shiny_source,
             layerData=if (nrow(hex)) 2 else 1)
  })
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

  ## We don't need to do this when we don't ggplotly(..., tooltip='text')
  ## Disable hovering on hexbin
  # p$x$source <- shiny_source
  # if (nrow(hex)) {
  #   # hexidx <- 2:(length(p$x$data) -1L)
  #   end.idx <- length(p$x$data)
  #   if (nrow(pts)) end.idx <- end.idx - 1L
  #   if (nrow(hlite)) end.idx <- end.idx - 1L
  #   hexidx <- 1:end.idx
  #   for (idx in hexidx) {
  #     p$x$data[[idx]]$hoverinfo <- 'none'
  #     # len <- length(p$x$data[[idx]]$x)
  #     # p$x$data[[idx]]$text <- sprintf('count: %d', rep(len))
  #   }
  # }

  config(p, collaborate=FALSE, displaylogo=FALSE)
}

mg_add_points <- function(gg, dat, color='black') {
  if (is.null(dat) || nrow(dat) == 0) {
    return(gg)
  }
  if ('symbol' %in% names(dat)) {
    gg <- gg + suppressWarnings({
      geom_point(aes(key=featureId,
                     text=paste0('Symbol: ', symbol, '<br>',
                                 xaxis, ': ', sprintf('%.3f', .xv), '<br>',
                                 yaxis, ': ', sprintf('%.3f', .yv))),
                 color=color, data=dat)
    })
  } else {
    gg <- gg +
      geom_point(aes(key=featureId,
                     text=paste0('featureId: ', featureId, '<br>',
                                 xaxis, ': ', sprintf('%.3f', .xv), '<br>',
                                 yaxis, ': ', sprintf('%.3f', .yv))),
                 color=color, data=dat)
  }
  gg
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
      x <- logFC(x, as.dt=TRUE)
      idx <- 'featureId'
    } else {
      x <- result(x, stats, as.dt=TRUE)
      idx <- 'idx'
      x[[idx]] <- encode_gskey(x)
    }
  }

  if (!is.data.frame(x)) {
    stop('x should have been converted into a data.frame by now')
  }
  x <- setDF(data.table::copy(x))
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
