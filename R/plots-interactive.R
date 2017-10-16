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
##'   "psig", and "sig". "notsig" implies that the FDR >= 10\%, "psig" means that
##'   FDR <= 10\%, but the logFC is "unremarkable" (< 1), and "sig" means
##'   that both the FDR <= 10\% and the logFC >= 1
##' @param with.data if \code{TRUE}, the data used for the plot is added to
##'   the outgoing rbokeh plot object (list) as \code{$data}. Default is
##'   \code{FALSE}.
##' @return the rbokeh (list) plot object
iplot <- function(x, y, j, value=c('logFC', 't'),
                  type=c('density', 'boxplot'),
                  tools=c('wheel_zoom', 'box_select', 'reset', 'save'),
                  main=NULL, with.legend=TRUE, with.data=FALSE,
                  shiny_source='mggenes', width=NULL, height=NULL,
                  ggtheme=theme_bw(), trim=0.005, ...) {
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

  # silence R CMD check NOTEs
  val <- NULL
  dat <- local({
    lfc <- copy(logFC(x, as.dt=TRUE))[, group := 'bg']
    lfc[, val := lfc[[value]]]
    gs.stats <- copy(geneSet(x, y, j, as.dt=TRUE))[, group := 'geneset']
    gs.stats[, val := gs.stats[[value]]]
    kcols <- intersect(names(gs.stats), names(lfc))

    out <- rbind(lfc[, kcols, with=FALSE], gs.stats[, kcols, with=FALSE])
    out[, significant := {
      tmp <- ifelse(significant, 'sig', 'notsig')
      ifelse(padj < 0.10 & tmp == 'notsig', 'psig', tmp)
    }]
  })

  if (type == 'density') {
    # out <- iplot.density.rbokeh(x, y, j, value, main, dat=dat,
    #                             with.legend=with.legend, tools=tools,
    #                             with.data=with.data, ...)
    out <- iplot.density.plotly(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, tools=tools,
                                with.data=with.data, shiny_source=shiny_source,
                                ggtheme=ggtheme, trim=trim, ...)
  } else if (type == 'boxplot') {
    # out <- iplot.boxplot.rbokeh(x, y, j, value, main, dat=dat,
    #                             with.legend=with.legend, tools=tools,
    #                             with.data=with.data, ...)
    out <- iplot.boxplot.plotly(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, tools=tools,
                                with.data=with.data, shiny_source=shiny_source,
                                width=width, height=height, ggtheme=ggtheme,
                                trim=trim, ...)
  } else if (type == 'volcano') {
    # out <- iplot.volcano.rbokeh(x, y, j, value, main, dat=dat,
    #                             with.legend=with.legend, tools=tools,
    #                             with.data=with.data, ...)
    # out <- iplot.volcano.plotly(x, y, j, value, main, dat=dat,
    #                             with.legend=with.legend, tools=tools,
    #                             with.data=with.data, width=width, height=height,
    #                             shiny_source=shiny_source, ggtheme=ggtheme, ...)
  }

  out
}

## plotly ======================================================================

##' @rdname iplot
iplot.density.plotly <- function(x, y, j, value, main, dat, with.legend=TRUE,
                                 with.points=TRUE,  with.data=FALSE,
                                 shiny_source='mggenes',
                                 legend.pos=c('inside', 'outside'),
                                 height=NULL, width=NULL, trim=0.02,
                                 square=TRUE, ...) {
  stopifnot(is(x, 'MultiGSEAResult'))
  legend.pos <- match.arg(legend.pos)
  gs.dat <- subset(dat, group == 'geneset') %>% setDF
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

  bgd <- density(bg$val)
  gsd <- density(gs.dat$val)
  lmeta <- list(width=3)

  if (is.numeric(trim) && trim != 0) {
    xrange <- quantile(dat$val, c(trim, 1 - trim))
    xrange[1] <- min(xrange[1], min(gs.dat$val), min(gsd$x))
    xrange[2] <- max(xrange[1], max(gs.dat$val), max(gsd$x))
  } else {
    xrange <- range(dat$val)
  }
  if (square) {
    extreme <- max(abs(xrange))
    xrange <- c(-extreme, extreme)
  }

  p <- plot_ly(source=shiny_source, width=width, height=height) %>%
    add_lines(x=bgd$x, y=bgd$y, name='All Genes', hoverinfo='none', line=lmeta) %>%
    add_lines(x=gsd$x, y=gsd$y, name='Geneset', hoverinfo='none', line=lmeta) %>%
    layout(xaxis=list(title="logFC", range=xrange),
           yaxis=list(title="Density"),
           title=main, dragmode="select")
  if ('symbol' %in% names(gs.dat) && with.points) {
    p <- add_markers(p, x=~val, y=~y, key=~featureId, data=gs.dat, name="Genes",
                     hoverinfo='text',
                     text=~paste0('Symbol: ', symbol, '<br>',
                                  'logFC: ', sprintf('%.3f', logFC), '<br>',
                                  'FDR: ', sprintf('%.3f', padj)))
  } else if (with.points) {
    p <- add_markers(p, x=~val, y=~y, key=~featureId, data=gs.dat, name="Genes",
                     text=~paste0('featureId: ', featureId, '<br>',
                                  'logFC: ', logFC, '<br>',
                                  'FDR: ', padj))
  }
  if (legend.pos == 'inside') {
    p <- layout(p, legend=list(x=0.75, y=1))
  }

  config(p, collaborate=FALSE, displaylogo=FALSE)
}

iplot.boxplot.plotly <- function(x, y, j, value, main, dat, with.legend=TRUE,
                                 with.points=TRUE, with.data=FALSE,
                                 shiny_source='mggenes', height=NULL,
                                 width=NULL, ggtheme=theme_bw(), trim=0.02,
                                 ...) {
  is.gs <- dat[['group']] == 'geneset'
  gs <- subset(dat, is.gs) %>% setDF
  bg <- subset(dat, !is.gs) %>% setDF
  n.gs <- sum(is.gs)
  if (value == 't') {
    value <- 't-statistic'
  }

  # silence R CMD check NOTEs
  val <- symbol <- NULL
  gg <- ggplot(bg, aes(group, val)) +
    geom_boxplot(data=bg) +
    geom_boxplot(outlier.shape=NA, data=gs)
  if ('symbol' %in% names(bg) && with.points) {
    gg <- gg +
      suppressWarnings({
        geom_jitter(aes(key=featureId,
                        text=paste0('Symbol: ', symbol, '<br>',
                                    'logFC: ', sprintf('%.3f', logFC), '<br>',
                                    'FDR: ', sprintf('%.3f', padj))),
                    data=gs, width=0.2)
      })

  } else if (with.points) {
    gg <- gg +
      suppressWarnings({
        geom_jitter(aes(key=featureId,
                        text=paste(paste0('featureId: ', featureId, '<br>',
                                          'logFC: ', logFC, '<br>',
                                          'FDR: ', padj))),
                    data=gs, width=0.2)
      })
  }

  if (is(ggtheme, 'theme')) {
    gg <- gg + ggtheme
  }

  ## ggplotly keeps suggesting to use the github/ggplot2
  p <- suppressMessages(ggplotly(gg, width=width, height=height, tooltip='text')) %>%
    layout(yaxis=list(title=value), dragmode="select") %>%
    plotly_build

  p$x$source <- shiny_source
  ## Hacks to hide hover events on boxplots and remove outliers
  ## https://plot.ly/ggplot2/box-plots/#outliers
  p$x$data[[1]]$hoverinfo <- 'none'
  p$x$data[[1]]$marker <- list(opacity=0)
  p$x$data[[2]]$hoverinfo <- 'none'
  p$x$data[[2]]$marker <- list(opacity=0)
  config(p, collaborate=FALSE, displaylogo=FALSE)
}

## rbokeh ----------------------------------------------------------------------

if (FALSE) {
## @rdname iplot
iplot.boxplot.rbokeh <- function(x, y, j, value, main, dat, with.legend=TRUE,
                                 tools=c('wheel_zoom', 'box_select', 'reset', 'save'),
                                 with.data=FALSE, ...) {
  gs <- subset(dat, group == 'geneset') %>%
   transform(jgrp=catjitter(group, 0.5), stringsAsFactors=FALSE)
  gs[['__index']] <- gs[['featureId']]
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

## @rdname iplot
iplot.density.rbokeh <- function(x, y, j, value, main, dat, with.legend=TRUE,
                                 tools=c('wheel_zoom', 'box_select', 'reset', 'save'),
                                 with.data=FALSE, ...) {
  stopifnot(is(x, 'MultiGSEAResult'))

  dat[['__index']] <- dat[['featureId']]
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
}
