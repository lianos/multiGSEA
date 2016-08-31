##' Create interactive GSEA plots
##'
##' @export
iplot <- function(x, y, j, value=c('logFC', 't'),
                  type=c('density', 'boxplot', 'volcano'),
                  main=NULL, ...) {
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

  if (type == 'volcano') {
    out <- iplot.volcano.rbokeh(x, y, j, value, main, dat=dat, ...)
  } else if (type == 'density') {
    out <- iplot.density.rbokeh(x, y, j, value, main, dat=dat, ...)
  } else if (type == 'boxplot') {
    # out <- iplot.boxplot.gg(x, y, j, value, main, dat=dat, ...)
    out <- iplot.boxplot.rbokeh(x, y, j, value, main, dat=dat, ...)
    # out <- iplot.boxplot.gg(x, y, j, value, main, dat=dat, ...)
  }

  out
}

## boxplot ---------------------------------------------------------------------
iplot.boxplot <- function(x, y, j, value, main, dat, ...) {
  bg <- subset(dat, group == 'bg')
  gs <- subset(dat, group == 'geneset')
  plt <- plot_ly(bg, y=val, type='box') %>%
    add_trace(y=val, type='box', boxpoints='all', jitter=0.3, pointpos=0,
              data=gs)

}

iplot.boxplot.rbokeh <- function(x, y, j, value, main, dat, ...) {
  gs <- subset(dat, group == 'geneset') %>%
   transform(jgrp=catjitter(group, 0.5), stringsAsFactors=FALSE)
  bg <- subset(dat, group != 'geneset')
  n.gs <- sum(dat$group == 'geneset')
  p <- figure(xlab=sprintf("Gene Set Group (%d genes)", n.gs)) %>%
    ly_boxplot(x="group", y="val", data=bg, fill_color='white') %>%
    ly_boxplot(x="group", y="val", data=gs,  fill_color='white', with_outliers=FALSE) %>%
    ly_points(x="jgrp", y="val", data=gs, color="significant",
              hover=list(symbol, logFC, padj))
  p
}


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

## density ---------------------------------------------------------------------
## This is super slow compared to the ggplot version?
iplot.density <- function(x, y, j, value, main, dat, ...) {
  stopifnot(is(x, 'MultiGSEAResult'))

  dens <- dat[, {
    xd <- density(val)
    list(x=xd$x, y=xd$y)
  }, by='group']

  dens.cols <- c('bg'='black', 'geneset'='red')
  point.cols <- c('grey', 'red')
  gs.dat <- filter(dat, group == 'geneset')
  gs.dat$y <- jitter(rep(0.1, nrow(gs.dat)), amount=0.05)

  pp <- suppressWarnings({
    ## plotly is complaining about RColorBrewer not being able to produce
    ## a palette for two colors (bg and geneset) even though I am sending in
    ## a predefined color palette.
    plot_ly() %>%
      add_trace(x=x, y=y, color=group, colors=dens.cols, mode='lines', data=dens,
                hoverinfo='none', line=list(width=2.5)) %>%
      add_trace(x=val, y=y, mode='markers',
                # symbol=significant,
                color=significant,
                colors=point.cols,
                hoverinfo='text',
                text=sprintf("symbol: %s<br>logFC: %.2f<br>FDR: %.2f",
                             gs.dat$symbol, gs.dat$logFC, gs.dat$padj),

                data=gs.dat) %>%
      layout(title=main,
             xaxis=list(title=sprintf("%s<br>(%d genes)", value, nrow(gs.dat))),
             yaxis=list(title="Density"))
  })
  pp
}

iplot.density.gg <- function(x, y, j, value, main, dat, ...) {
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

iplot.density.rbokeh <- function(x, y, j, value, main, dat, ...) {
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

  # gs <- filter(dat, group == 'geneset') %>%
  #   mutate(jgrp=catjitter(group, 0.5))
  bg <- filter(dat, group == 'bg')

  p <- figure() %>%
    ly_density(x="val", data=bg, color="black", width=2) %>%
    ly_density(x="val", data=gs.dat, color="red", width=2) %>%
    ly_points(x="val", y="y", data=gs.dat, color="significant",
               hover=list(symbol, logFC, padj))
}


iplot.cdf.gg <- function(x, y, j, value, main, dat, ...) {
  stopifnot(is(x, 'MultiGSEAResult'))
  gs.dat <- subset(dat, group == 'geneset')
  gs.dat$significant <- ifelse(gs.dat$significant, 'sig', 'notsig')
  gs.dat$significant <- ifelse(gs.dat$padj < 0.10 & gs.dat$significant == 'notsig',
                               'psig', gs.dat$significant)
  gs.dat$y <- c('notsig'=0.1, 'psig'=0.2, 'sig'=0.3)[gs.dat$significant]
  gs.dat$y <- 0.15
  cols <- c('bg'='black', 'geneset'='red', 'notsig'='grey', 'psig'='lightblue', sig='darkblue')

  gg <- ggplot(dat, aes(val, color=group)) +
    stat_ecdf(lwd=1, geom='step') +
    xlab(sprintf("%s (%d genes)", value, sum(dat$group == 'geneset'))) +
    ggtitle(main) +
    scale_color_manual(values=cols)

  ggplotly(gg)

}


iplot.volcano <- function(x, y, j, value, main, dat,  ...) {
  stopifnot(is(x, 'MultiGSEAResult'))
}




