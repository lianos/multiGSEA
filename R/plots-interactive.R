##' Create interactive GSEA plots
##'
##' export
iplot <- function(x, y, j, value=c('logFC', 't'), type=c('density', 'volcano'),
                  main=NULL, ...) {
  stopifnot(is(x, 'MultiGSEAResult'))
  type <- match.arg(type)
  value <- match.arg(value)

  if (is.null(main)) {
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

  if (type == 'volcano') {
    out <- iplot.volcano(x, y, j, value, main, dat=dat, ...)
  } else if (type == 'density') {
    out <- iplot.density.gg(x, y, j, value, main, dat=dat, ...)
  }

  out
}

## This is super slow compared to the ggplot version?
iplot.density <- function(x, y, j, value, main, dat, ...) {
  if (FALSE) {
    x <- xmg; y <- 'h'; j <- 'HALLMARK_E2F_TARGETS'; value <- 'logFC'
  }
  stopifnot(is(x, 'MultiGSEAResult'))

  dens <- dat[, {
    xd <- density(val)
    list(x=xd$x, y=xd$y)
  }, by='group']

  dens.cols <- c('bg'='black', 'geneset'='red')
  point.cols <- c('grey', 'red')
  gs.dat <- filter(dat, group == 'geneset')
  gs.dat$y <- jitter(rep(0.1, nrow(gs.dat)), amount=0.05)
  gs.dat$significant <- ifelse(gs.dat$significant, 'sig', 'notsig')


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
  if (FALSE) {
    x <- xmg; y <- 'h'; j <- 'HALLMARK_E2F_TARGETS'; value <- 'logFC'
  }
  stopifnot(is(x, 'MultiGSEAResult'))

  gs.dat <- subset(dat, group == 'geneset')
  gs.dat$significant <- ifelse(gs.dat$significant, 'sig', 'notsig')
  gs.dat$significant <- ifelse(gs.dat$padj < 0.10 & gs.dat$significant == 'notsig',
                               'psig', gs.dat$significant)
  gs.dat$y <- c('notsig'=0.1, 'psig'=0.2, 'sig'=0.3)[gs.dat$significant]
  gs.dat$y <- 0.15
  cols <- c('bg'='black', 'geneset'='red', 'notsig'='grey', 'psig'='lightblue', sig='darkblue')

  gg <- ggplot(dat, aes(val, color=group)) +
    stat_density(lwd=1, geom='line') +
    geom_point(aes(val, y, color=significant,# shape=significant,
                   text=sprintf("symbol: %s<br>logFC: %.2f<br>FDR: %.2f", symbol, logFC, padj)),
               data=gs.dat,
               position=position_jitter(width=0, height=0.1)) +
    xlab(sprintf("%s (%d genes)", value, sum(dat$group == 'geneset'))) +
    ggtitle(main) +
    scale_color_manual(values=cols)

  ggplotly(gg)
}

iplot.volcano <- function(x, y, j, value, main, dat,  ...) {
  if (FALSE) {
    x <- xmg; y <- 'h'; j <- 'HALLMARK_E2F_TARGETS'; value <- 'logFC'
  }
  stopifnot(is(x, 'MultiGSEAResult'))
}




