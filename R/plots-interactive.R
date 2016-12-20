##' Create interactive GSEA plots
##'
##' @export
iplot <- function(x, y, j, value=c('logFC', 't'),
                  type=c('density', 'boxplot', 'volcano'),
                  tools=c('wheel_zoom', 'box_select', 'reset', 'save'),
                  main=NULL, with.legend=TRUE, with.data=FALSE, ...) {
  if (FALSE) {
    x <- xmg; y <- 'h'; j <- 'HALLMARK_E2F_TARGETS'; value <- 'logFC';
    main <- NULL; type <- 'boxplot'; with.legend <- TRUE
  }

  stopifnot(is(x, 'MultiGSEAResult'))
  type <- match.arg(type)
  value <- match.arg(value)

  if (type == 'volcano') {
    stop("volcano not yet implemented")
  }

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
