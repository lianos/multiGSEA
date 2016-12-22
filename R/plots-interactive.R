##' Visualize gene level behavior of genes within a geneset across a contrast.
##'
##' @description
##' It is informative to look at the individual log fold changes of the genes
##' within a gene set to explore the degree to which they (1) are coherent with
##' respect to each other; and (2) see how the compare to the background
##' distribution of log fold changes of the entire transcriptome.
##'
##' You can visualize this behavior via a \code{type='density'} plot, or a
##' \code{type='boxplot'}. It is also common to plot either the individual
##' log fold changes (\code{value="logFC"}) or t-statistics (\code{value="t"}).
##'
##' @rdname iplot
##' @export
##'
##' @param x A \code{MultiGSEAResult} object
##' @param y the name of the gene set collection
##' @param j the name of the gene set name
##' @param value plot individual log fold changes (default) or t-statistics
##'   (\code{"logFC"} or \code{"t"}, respectively)
##' @param type plot the distributions as a \code{"density"} plot or
##'   \code{"boxplot"}
##' @param tools the tools to display in the rbokeh plot
##' @param main A title to display. If not specified, the gene set name
##'   will be used, otherwise you can pass in a custom title, or \code{NULL}
##'   will disable the title altogether.
##' @param with.legend Draws a legend to map point color to meaning. There are
##'   three levels a point (gene level statistic) can be color as, "notsig",
##'   "psig", and "sig". "notsig" implies that the FDR \gt 10%%, "psig" means that
##'   FDR \lte 10%%, but the logFC is "unremarkable" (\lt < 1), and "sig" means
##'   that both the FDR \lte 10% and the logFC \gte 1
##' @param with.data if \code{TRUE}, the data used for the plot is added to
##'   the outgoing rbokeh plot object (list) as \code{$data}. Default is
##'   \code{FALSE}.
##' @return the rbokeh (list) plot object
iplot <- function(x, y, j, value=c('logFC', 't'),
                  type=c('density', 'boxplot'),
                  tools=c('wheel_zoom', 'box_select', 'reset', 'save'),
                  main=NULL, with.legend=TRUE, with.data=FALSE, ...) {
  if (FALSE) {
    x <- xmg; y <- 'h'; j <- 'HALLMARK_E2F_TARGETS'; value <- 'logFC';
    main <- NULL; type <- 'boxplot'; with.legend <- TRUE
  }

  stopifnot(is(x, 'MultiGSEAResult'))
  type <- match.arg(type)
  value <- match.arg(value)

  if (missing(main)) {
    main <- sprintf("%s (%s)", j, y)
  }

  dat <- local({
    lfc <- copy(logFC(x, .external=FALSE))[, group := 'bg']
    lfc[, val := lfc[[value]]]
    gs.stats <- geneSet(x, y, j, .external=FALSE)[, group := 'geneset']
    gs.stats[, val := gs.stats[[value]]]
    kcols <- intersect(names(gs.stats), names(lfc))
    out <- bind_rows(select_(lfc, .dots=kcols), select_(gs.stats, .dots=kcols))
    out[, significant := {
      tmp <- ifelse(significant, 'sig', 'notsig')
      ifelse(padj < 0.10 & tmp == 'notsig', 'psig', tmp)
    }]
  })

  if (type == 'density') {
    out <- iplot.density.rbokeh(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, tools=tools,
                                with.data=with.data, ...)
  } else if (type == 'boxplot') {
    out <- iplot.boxplot.rbokeh(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, tools=tools,
                                with.data=with.data, ...)
  } else if (type == 'volcano') {
    out <- iplot.volcano.rbokeh(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, tools=tools,
                                with.data=with.data, ...)
  }

  out
}

## rbokeh ----------------------------------------------------------------------
##' @rdname iplot
iplot.boxplot.rbokeh <- function(x, y, j, value, main, dat, with.legend=TRUE,
                                 tools=c('wheel_zoom', 'box_select', 'reset', 'save'),
                                 with.data=FALSE, ...) {
  gs <- subset(dat, group == 'geneset') %>%
   transform(jgrp=catjitter(group, 0.5), stringsAsFactors=FALSE)
  bg <- subset(dat, group != 'geneset')
  n.gs <- sum(dat$group == 'geneset')
  if (value == 't') {
    value <- 't-statistic'
  }
  p <- figure(xlab=sprintf("Gene Set Group (%d genes)", n.gs), ylab=value,
              tools=tools) %>%
    ly_boxplot(x="group", y="val", data=bg, fill_color='white',
               outlier_size=5) %>%
    ly_boxplot(x="group", y="val", data=gs,  fill_color='white',
               outlier_size=NA)
  if ('symbol' %in% names(gs)) {
    p <- ly_points(p, x="jgrp", y="val", data=gs, color="significant", size=5,
                   hover=list(symbol, logFC, padj), lname='points',
                   legend=with.legend)
  } else {
    p <- ly_points(p, x="jgrp", y="val", data=gs, color="significant", size=5,
                   hover=list(featureId, logFC, padj), lname='points',
                   legend=with.legend)
  }

  if (with.data) {
    p$data <- gs
  }

  p
}

##' @rdname iplot
iplot.density.rbokeh <- function(x, y, j, value, main, dat, with.legend=TRUE,
                                 tools=c('wheel_zoom', 'box_select', 'reset', 'save'),
                                 with.data=FALSE, ...) {
  stopifnot(is(x, 'MultiGSEAResult'))

  gs.dat <- subset(dat, group == 'geneset')
  cols <- c('bg'='black', 'geneset'='red',
            'notsig'='grey', 'psig'='lightblue', 'sig'='darkblue')

  if (value == 't') {
    value <- 't-statistic'
    # gs.dat$y <- 0.005
    gs.dat$y <- 0.0015 + runif(nrow(gs.dat), 0, 0.004)
    jitter <- 0.005
  } else {
    ## gs.dat$y <- c('notsig'=0.1, 'psig'=0.2, 'sig'=0.3)[gs.dat$significant]
    gs.dat$y <- 0.005 + runif(nrow(gs.dat), 0, 0.040)
    jitter <- 0.05
  }

  bg <- subset(dat, group == 'bg')

  p <- figure(xlab=sprintf("%s (%d genes)", value, nrow(gs.dat)),
              tools=tools) %>%
    ly_density(x="val", data=bg, color="black", width=3) %>%
    ly_density(x="val", data=gs.dat, color="red", width=3)
  if ('symbol' %in% names(gs.dat)) {
    p <- ly_points(p, x="val", y="y", data=gs.dat, color="significant", size=5,
                   hover=list(symbol, logFC, padj), legend=with.legend,
                   lname='points')
  } else {
    p <- ly_points(p, x="val", y="y", data=gs.dat, color="significant", size=5,
                   hover=list(featureId, logFC, padj), legend=with.legend,
                   lname='points')
  }

  if (with.data) {
    p$data <- gs.dat
  }

  p
}

if (FALSE) {
## ggplotly --------------------------------------------------------------------
iplot.boxplot.gg <- function(x, y, j, value, main, dat, ...) {
  cols <- c('notsig'='grey', 'psig'='lightblue', 'sig'='darkblue')
  bg <- subset(dat, group == 'bg')
  gs <- subset(dat, group == 'geneset')
  plabel <- 'symbol: %s<br>logFC: %.2f<br>FDR: %.2f'
  value <- if (value == 't') 't-statistic' else value
  gg <- ggplot(dat, aes(group, val)) +
    geom_boxplot(data=bg) +
    geom_boxplot(data=gs, outlier.shape=NA) +
    geom_point(aes(color=significant,
                   text=sprintf(plabel, symbol, logFC, padj)),
               data=gs,
               position=position_jitter(width=0.25)) +
    ylab(value) +
    xlab(sprintf("Gene Set Group<br>(%d genes)", sum(dat$group == 'geneset'))) +
    scale_color_manual(values=cols)
  if (!is.null(main)) {
    gg <- gg + ggtitle(main)
  }
  ggplotly(gg)
}

iplot.density.gg <- function(x, y, j, value, main, dat, with.legend=TRUE, ...) {
  stopifnot(is(x, 'MultiGSEAResult'))

  gs.dat <- subset(dat, group == 'geneset')
  cols <- c('bg'='black', 'geneset'='red',
            'notsig'='grey', 'psig'='lightblue', 'sig'='darkblue')

  if (value == 't') {
    value <- 't-statistic'
    gs.dat$y <- 0.005
    jitter <- 0.005
  } else {
    ## gs.dat$y <- c('notsig'=0.1, 'psig'=0.2, 'sig'=0.3)[gs.dat$significant]
    gs.dat$y <- 0.1
    jitter <- 0.05
  }

  plabel <- "symbol: %s<br>logFC: %.2f<br>FDR: %.2f"
  gg <- ggplot(dat, aes(val, color=group)) +
    stat_density(lwd=1, geom='line') +
    geom_point(aes(val, y, color=significant,
                   text=sprintf(plabel, symbol, logFC, padj)),
               data=gs.dat,
               position=position_jitter(width=0, height=jitter)) +
    xlab(sprintf("%s (%d genes)", value, sum(dat$group == 'geneset'))) +
    scale_color_manual(values=cols)
  if (!is.null(main)) {
    gg <- gg + ggtitle(main)
  }
  ggplotly(gg)
}

}
