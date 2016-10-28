##' correlation colors for corplot
col.pairs <- circlize::colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'), 0.5)

##' Plots the correlation of genes in signature.
##'
##' @export
##' @rdname corplot
##' @param E sample-by-gene matrix from \code{\link{fetchExpression}} to plot.
##'   Note that this is transpose of "the usual" gene x sample matrix.
##' @param title The title of the plot
##' @param cluster \code{logical} indicating whether or not to shuffle genes
##'   around into some clustering.
##' @param col.point the color of the points in the scatterplots
##' @param diag.distro show the distribution of values on the diagnols?
##' @return nothing, just creates the plot
corplot <- function(E, title, cluster=FALSE, col.point='#00000066',
                    diag.distro=TRUE) {
  if (missing(title)) {
    title <- 'Pairs Plot'
  }
  
  if (is.null(E)) {
    return(plot(3, 3, pch=16, xlim=c(1, 5), ylim=c(1, 5)))
  }
  
  if (cluster) {
    cors <- cor(E, method='spearman', use='na.or.complete')
    dists <- as.dist((1 - cors) / 2)
    hc <- hclust(dists)
    dendro <- as.dendrogram(hc)
    idxs <- rev(order.dendrogram(dendro))
    E <- E[, idxs]
  }
  
  suppressWarnings({
    pairs(E, col.point=col.point,
          lower.panel=panel.cor,
          upper.panel=panel.spoints,
          diag.panel=if (diag.distro) panel.hist else NULL,
          gap=0.2, pch=16,
          main=title)
  })
  
  invisible(NULL)
}

## Helper functions ------------------------------------------------------------

##' Plots the correlation stats with colored background for the lower triangle
##' of pairs plot.
##' @rdname corplot
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par(c("usr"))
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method='spearman', use='na.or.complete')
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if (missing(cex.cor)) {
    cex.cor <- 0.8/strwidth(txt)
  }
  
  if (!is.na(r)) {
    bg.col <- col.pairs(r)
  } else {
    bg.col <- "#3f3f3f33"
  }
  
  rect(0, 0, 1, 1, col=bg.col)
  text(0.5, 0.5, txt, cex = cex.cor)
}

##' Plots the expression distribution as histogrm (mai diagonal of pair plot)
##' @rdname corplot
panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr=c(usr[1:2], 0, 1.5))
  h <- hist(x, breaks=20, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  suppressWarnings(rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...))
}

##' Scatter plot for upper triangle of pairs plot
##' @rdname corplot
panel.spoints <- function(x, y, col=par("col"), bg=NA, pch=par("pch"), cex=1,
                          col.smooth="red", span=2/3, iter=3,
                          col.point="#00000022", ...) {
  ## points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  points(x, y, pch=16, col=col.point, bg=bg, cex=cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok))
    suppressWarnings({
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
            col = col.smooth, ...)
    })
}
