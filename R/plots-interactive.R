##' Create interactive GSEA plots
##'
##' @export
iplot <- function(x, y, j, value=c('logFC', 't'),
                  type=c('density', 'boxplot', 'volcano'),
                  main=NULL, with.legend=TRUE, ...) {
  if (FALSE) {
    x <- xmg; y <- 'h'; j <- 'HALLMARK_E2F_TARGETS'; value <- 'logFC';
    main <- NULL;
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
    lfc <- transform(logFC(x, .external=FALSE), group='bg')
    lfc$val <- lfc[[value]]
    gs.stats <- transform(geneSet(x, y, j, .external=FALSE),
                          group='geneset', collection=NULL, name=NULL)
    gs.stats$val <- gs.stats[[value]]
    rbindlist(list(lfc, gs.stats))
  })
  dat$significant <- ifelse(dat$significant, 'sig', 'notsig')
  dat$significant <- ifelse(dat$padj < 0.10 & dat$significant == 'notsig',
                            'psig', dat$significant)

  if (type == 'density') {
    out <- iplot.density.rbokeh(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, ...)
  } else if (type == 'boxplot') {
    out <- iplot.boxplot.rbokeh(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, ...)
  } else if (type == 'volcano') {
    out <- iplot.volcano.rbokeh(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, ...)
  }

  out
}

## rbokeh ----------------------------------------------------------------------
iplot.boxplot.rbokeh <- function(x, y, j, value, main, dat, with.legend=TRUE, ...) {
  gs <- subset(dat, group == 'geneset') %>%
   transform(jgrp=catjitter(group, 0.5), stringsAsFactors=FALSE)
  bg <- subset(dat, group != 'geneset')
  n.gs <- sum(dat$group == 'geneset')
  p <- figure(xlab=sprintf("Gene Set Group (%d genes)", n.gs)) %>%
    ly_boxplot(x="group", y="val", data=bg, fill_color='white',
               outlier_size=5) %>%
    ly_boxplot(x="group", y="val", data=gs,  fill_color='white',
               outlier_size=NA) %>%
    ly_points(x="jgrp", y="val", data=gs, color="significant", size=5,
              hover=list(symbol, logFC, padj), legend=with.legend)
  p
}

iplot.density.rbokeh <- function(x, y, j, value, main, dat, with.legend=TRUE,
                                 ...) {
  stopifnot(is(x, 'MultiGSEAResult'))

  gs.dat <- subset(dat, group == 'geneset')
  cols <- c('bg'='black', 'geneset'='red',
            'notsig'='grey', 'psig'='lightblue', 'sig'='darkblue')

  if (value == 't') {
    value <- 't-statistic'
    # gs.dat$y <- 0.005
    gs.dat$y <- 0.005 + runif(nrow(gs.dat), 0, 0.040)
    jitter <- 0.005
  } else {
    ## gs.dat$y <- c('notsig'=0.1, 'psig'=0.2, 'sig'=0.3)[gs.dat$significant]
    gs.dat$y <- 0.005 + runif(nrow(gs.dat), 0, 0.040)
    jitter <- 0.05
  }

  bg <- subset(dat, group == 'bg')

  p <- figure() %>%
    ly_density(x="val", data=bg, color="black", width=3) %>%
    ly_density(x="val", data=gs.dat, color="red", width=3) %>%
    ly_points(x="val", y="y", data=gs.dat, color="significant", size=5,
              hover=list(symbol, logFC, padj), legend=with.legend)
  p
}

## ggplotly --------------------------------------------------------------------
iplot.boxplot.gg <- function(x, y, j, value, main, dat, ...) {
  cols <- c('notsig'='grey', 'psig'='lightblue', 'sig'='darkblue')
  bg <- filter(dat, group == 'bg')
  gs <- filter(dat, group == 'geneset')
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
