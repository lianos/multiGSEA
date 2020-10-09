#' Visualize gene level behavior of genes within a geneset across a contrast.
#'
#' @description
#' It is informative to look at the individual log fold changes of the genes
#' within a gene set to explore the degree to which they (1) are coherent with
#' respect to each other; and (2) see how the compare to the background
#' distribution of log fold changes of the entire transcriptome.
#'
#' You can visualize this behavior via a `type = "density"` plot, or a
#' `type = "boxplot". It is also common to plot either the individual
#' log fold changes `value = "logFC"` or t-statistics `value = "t"`.
#'
#' @rdname iplot
#' @export
#' @importFrom ggplot2 theme_bw
#'
#' @param x A [MultiGSEAResult()] object
#' @param y the name of the gene set collection
#' @param j the name of the gene set name
#' @param value A string indicating the column name for the value of the
#'   gene-level metadata to plot. Default is `"logFC"`. Anoter often used choice
#'   might also be `"t"`, to plot t-statistics (if they're in the result). But
#'   this can be any numeric column found in the data.frame returned by
#'   `geneSet(x, y, j)`. If this is a named string (vector), then the value in
#'   `names(value)` will be used on the axis when plotted.
#' @param type plot the distributions as a `"density"` plot or `"boxplot"`.
#' @param tools the tools to display in the rbokeh plot
#' @param main A title to display. If not specified, the gene set name
#'   will be used, otherwise you can pass in a custom title, or `NULL`
#'   will disable the title altogether.
#' @param with.legend Draws a legend to map point color to meaning. There are
#'   three levels a point (gene level statistic) can be color as, "notsig",
#'   "psig", and "sig". "notsig" implies that the FDR >= 10%, "psig" means that
#'   FDR <= 10%, but the logFC is "unremarkable" (< 1), and "sig" means
#'   that both the FDR <= 10% and the logFC >= 1
#' @return the ploty plot ojbect
#' @examples
#' mgr <- exampleMultiGSEAResult()
#' iplot(mgr, "c2", "BURTON_ADIPOGENESIS_PEAK_AT_2HR", c("t-statistic" = "t"),
#'       type = "density")
#' iplot(mgr, "c2", "BURTON_ADIPOGENESIS_PEAK_AT_2HR", c("t-statistic" = "t"),
#'       type = "gsea")
iplot <- function(x, y, j, value = "logFC",
                  type=c("density", "gsea", "boxplot"),
                  tools=c('wheel_zoom', 'box_select', 'reset', 'save'),
                  main=NULL, with.legend=TRUE,
                  shiny_source='mggenes', width=NULL, height=NULL,
                  ggtheme=theme_bw(), trim=0.005, ...) {
  if (FALSE) {
    x <- xmg; y <- 'H'; j <- 'HALLMARK_E2F_TARGETS'; value <- 'logFC';
    main <- NULL; type <- 'boxplot'; with.legend <- TRUE
  }

  stopifnot(is(x, 'MultiGSEAResult'))
  type <- match.arg(type)
  lfc <- copy(logFC(x, as.dt = TRUE)[, group := "bg"])
  match.arg(value, colnames(lfc))
  stopifnot(is.numeric(lfc[[value]]))

  if (missing(main)) {
    main <- sprintf("%s (%s)", j, y)
  }

  if (type == "gsea") {
    gset <- geneSet(x, collection = y, name = j)
    return(iplot.gsea.plot(lfc, gset, rank_by = value, title = main, ...))
  }

  # silence R CMD check NOTEs
  val <- NULL
  dat <- local({
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
    out <- iplot.density.plotly(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, tools=tools,
                                shiny_source=shiny_source,
                                ggtheme=ggtheme, trim=trim, ...)
  } else if (type == 'boxplot') {
    out <- iplot.boxplot.plotly(x, y, j, value, main, dat=dat,
                                with.legend=with.legend, tools=tools,
                                shiny_source=shiny_source,
                                width=width, height=height, ggtheme=ggtheme,
                                trim=trim, ...)
  } else if (type == 'volcano') {
  }

  out
}

## plotly ======================================================================

#' Reimplementation of fgsea::plotEnrichment so we can add more interactive
#' bits, and also to highlight genes on the leading edge, and what not.
#'
#' Lots of code here is copied from fgsea
#'
#' @noRd
iplot.gsea.plot <- function(lfc, geneset, rank_by, title, gseaParam = 1,
                            ticksSize = 0.2, ..., .plot_default = FALSE,
                            .plot_static = FALSE) {
  if (!requireNamespace("fgsea")) stop("'fgsea' package required")

  # Setup params so we can just copy and paste fgsea::plotEnrichment code
  pathway <- geneset[["feature_id"]]
  stats <- setNames(lfc[[rank_by]], lfc[["feature_id"]])

  if (.plot_default) {
    return(fgsea::plotEnrichment(pathway, stats, gseaParam, ticksSize))
  }

  # fgsea::plotEnrichment ......................................................
  rnk <- rank(-stats)
  ord <- order(rnk)

  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))

  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)

  gseaRes <- fgsea::calcGseaStat(statsAdj, selectedStats = pathway,
                                  returnAllExtremes = TRUE)

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops

  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x=c(0, xs, n + 1), y=c(0, ys, 0))

  diff <- (max(tops) - min(bottoms)) / 8

  # Getting rid of NOTEs
  x <- y <- NULL

  # Creates a data.frame to use for the line segments drawn for each  # :custom
  # feature. This allows us to add useful hover information.          # :custom
  features <- data.frame(                                             # :custom
    x = pathway,                                                      # :custom
    y = -diff / 2,                                                    # :custom
    xend = pathway,                                                   # :custom
    yend = diff / 2,                                                  # :custom
    feature_id = names(statsAdj)[pathway])                            # :custom

  add.labels <- c("feature_id", "statistic", "statistic_adj")         # :custom
  xref <- match(features[["feature_id"]], geneset[["feature_id"]])    # :custom

  features[["statistic"]] <- geneset[[rank_by]][xref]                 # :custom
  features[["statistic_adj"]] <- unname(statsAdj)[pathway]            # :custom

  label <- intersect(c("symbol", "name"), colnames(geneset))          # :custom
  if (length(label)) {                                                # :custom
    features[["name"]] <- geneset[[label[1L]]][xref]                  # :custom
    add.labels <- c("name", add.labels)                               # :custom
  }

  stat.cols <- c("x", "y", "xend", "yend")                            # :custom
  for (cname in setdiff(colnames(features), stat.cols))   {           # :custom
    if (is.numeric(features[[cname]])) {                              # :custom
      features[[cname]] <- sprintf("%0.3f", features[[cname]])        # :custom
    } else {
      features[[cname]] <- as.character(features[[cname]])            # :custom
    }
  }

  features[["label"]] <- sapply(1:nrow(features), function(i) {       # :custom
    f <- features[i, add.labels]                                      # :custom
    paste(names(f), ":", unname(f[1,]), collapse = "<br>")            # :custom
  })

  if (is.null(names(rank_by))) {                                      # :custom
    xlabel <- sprintf("rank\n(by: %s)", rank_by)                      # :custom
  } else {
    xlabel <- sprintf("rank\n(by: %s)", names(rank_by))               # :custom
  }

  g <- ggplot2::ggplot(toPlot, ggplot2::aes(x=x, y=y)) +
    ggplot2::geom_point(color = "green", size = 0.1) +
    ggplot2::geom_hline(yintercept = max(tops), colour = "red",
                        linetype = "dashed") +
    ggplot2::geom_hline(yintercept = min(bottoms), colour = "red",
                        linetype = "dashed") +
    ggplot2::geom_hline(yintercept = 0, colour = "black") +
    ggplot2::geom_line(color = "green") + theme_bw() +
    # ggplot2::geom_segment(
    #   data = data.frame(x = pathway, label = ),
    #   mapping = ggplot2::aes(x = x, y = -diff / 2, xend = x, yend = diff / 2),
    #   size = ticksSize) +
    suppressWarnings({
      ggplot2::geom_segment(                                          # :custom
        data = features,                                              # :custom
        mapping = ggplot2::aes(                                       # :custom
          x = x, y = y, xend = xend, yend = yend, text = label),      # :custom
        size = ticksSize)                                             # :custom
    }) +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(
      x = xlabel,
      y = "enrichment score",
      title = title)

  if (!.plot_static) {
    g <- plotly::ggplotly(g, tooltip = "label") %>%
      layout(dragmode="select") %>%
      config(displaylogo=FALSE)
  }

  g
}

#' @noRd
#' @importFrom plotly add_markers add_lines config layout plot_ly
iplot.density.plotly <- function(x, y, j, value, main, dat, with.legend=TRUE,
                                 with.points=TRUE, shiny_source='mggenes',
                                 legend.pos=c('inside', 'outside'),
                                 height=NULL, width=NULL, trim=0.02,
                                 square=TRUE, ...) {
  stopifnot(is(x, 'MultiGSEAResult'))
  legend.pos <- match.arg(legend.pos)
  gs.dat <- subset(dat, group == 'geneset') %>% setDF
  cols <- c('bg'='black', 'geneset'='red',
            'notsig'='grey', 'psig'='lightblue', 'sig'='darkblue')

  if (value == 't') {
    gs.dat$y <- 0.0015 + runif(nrow(gs.dat), 0, 0.004)
    jitter <- 0.005
  } else {
    ## gs.dat$y <- c('notsig'=0.1, 'psig'=0.2, 'sig'=0.3)[gs.dat$significant]
    gs.dat$y <- 0.005 + runif(nrow(gs.dat), 0, 0.040)
    jitter <- 0.05
  }

  if (is.null(names(value))) {
    label <- if (value == "t") "t-statistic" else value
  } else {
    label <- names(value)
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
    layout(xaxis=list(title=label, range=xrange),
           yaxis=list(title="Density"),
           showlegend = with.legend, title=main, dragmode="select")
  if ('symbol' %in% names(gs.dat) && with.points) {
    p <- add_markers(p, x=~val, y=~y, key=~feature_id, data=gs.dat, name="Genes",
                     hoverinfo='text',
                     text=~paste0('Symbol: ', symbol, '<br>',
                                  'logFC: ', sprintf('%.3f', logFC), '<br>',
                                  'FDR: ', sprintf('%.3f', padj)))
  } else if (with.points) {
    p <- add_markers(p, x=~val, y=~y, key=~feature_id, data=gs.dat, name="Genes",
                     text=~paste0('feature_id: ', feature_id, '<br>',
                                  'logFC: ', logFC, '<br>',
                                  'FDR: ', padj))
  }
  if (legend.pos == 'inside') {
    p <- layout(p, legend=list(x=0.75, y=1))
  }

  config(p, displaylogo=FALSE)
}

#' @noRd
#' @importFrom plotly config ggplotly layout plotly_build
#' @importFrom ggplot2 aes geom_boxplot geom_jitter ggplot
iplot.boxplot.plotly <- function(x, y, j, value, main, dat, with.legend=TRUE,
                                 with.points=TRUE, shiny_source='mggenes',
                                 height=NULL, width=NULL, ggtheme=theme_bw(),
                                 trim=0.02, ...) {
  is.gs <- dat[['group']] == 'geneset'
  gs <- subset(dat, is.gs) %>% setDF
  bg <- subset(dat, !is.gs) %>% setDF
  n.gs <- sum(is.gs)

  if (is.null(names(value))) {
    label <- if (value == "t") "t-statistic" else value
  } else {
    label <- names(value)
  }

  # silence R CMD check NOTEs
  val <- symbol <- NULL
  gg <- ggplot(bg, aes(group, val)) +
    geom_boxplot(data=bg) +
    geom_boxplot(outlier.shape=NA, data=gs)
  if ('symbol' %in% names(bg) && with.points) {
    gg <- gg +
      suppressWarnings({
        geom_jitter(aes(key=feature_id,
                        text=paste0('Symbol: ', symbol, '<br>',
                                    'logFC: ', sprintf('%.3f', logFC), '<br>',
                                    'FDR: ', sprintf('%.3f', padj))),
                    data=gs, width=0.2)
      })

  } else if (with.points) {
    gg <- gg +
      suppressWarnings({
        geom_jitter(aes(key=feature_id,
                        text=paste(paste0('feature_id: ', feature_id, '<br>',
                                          'logFC: ', logFC, '<br>',
                                          'FDR: ', padj))),
                    data=gs, width=0.2)
      })
  }

  if (is(ggtheme, 'theme')) {
    gg <- gg + ggtheme
  }

  ## ggplotly keeps suggesting to use the github/ggplot2
  p <- suppressMessages({
    ggplotly(gg, width = width, height = height, tooltip = "text")
  })
  p <- layout(p, yaxis = list(title=label),
              dragmode = "select", showlegend = with.legend)
  p <- plotly_build(p)

  p$x$source <- shiny_source
  ## Hacks to hide hover events on boxplots and remove outliers
  ## https://plot.ly/ggplot2/box-plots/#outliers
  p$x$data[[1]]$hoverinfo <- 'none'
  p$x$data[[1]]$marker <- list(opacity=0)
  p$x$data[[2]]$hoverinfo <- 'none'
  p$x$data[[2]]$marker <- list(opacity=0)
  config(p, displaylogo = FALSE)
}
