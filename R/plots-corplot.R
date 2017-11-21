## correlation colors for corplot, red is high, blue is low.
col.pairs <- circlize::colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'), 0.5)

##' Plots the correlation among the columns of a numeric matrix.
##'
##' We assume that this is a sample x gene expression matrix, but it can
##' (of course) be any numeric matrix of your choosing. The column names appear
##' in the main diagonal of the plot. Note that you might prefer the corrplot
##' package for similar functionality, and this functionality is intentionally
##' named different from that..
##'
##' TODO: Add with.signature parameter to allow a box to plot the signature
##'   score of all genes in E.
##'
##' @export
##' @importFrom graphics smoothScatter text
##' @seealso \url{http://cran.r-project.org/package=corrplot}
##'
##' @param E the matrix used to plot a pairs correlation plot. The vectors used
##'   to assess all pairwise correlation should be \emph{in the columns} of the
##'   matrixl.
##' @param title The title of the plot
##' @param cluster \code{logical} indicating whether or not to shuffle genes
##'   around into some clustering.
##' @param col.point the color of the points in the scatterplots
##' @param diag.distro show the distribution of values on the diagnols?
##' @return nothing, just creates the plot
##'
##' @examples
##' x <- matrix(rnorm(1000), ncol=5)
##' corplot(x)
corplot <- function(E, title, cluster=FALSE, col.point='#00000066',
                    diag.distro=TRUE, smooth.scatter=nrow(E) > 400, ...) {
  E <- as_matrix(E)
  if (missing(title)) {
    title <- 'Pairs Plot'
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
          main=title, smooth.scatter=smooth.scatter, ...)
  })

  invisible(NULL)
}

## Helper functions ------------------------------------------------------------

## Helper function that fills in the lower left part of the output plot. The
## correlation between the two vectors is calculated and written in the
## corresponding panel. The background of the panel is colored in accordance
## to the strength and direction of the correlation.
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

## Helper function to draw the histogram on the diagonal of the corplot when
## `diag.distro == TRUE'
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

## Helper function to that generates the scatter plot of points in the upper
## right diagonal of the corplot.
panel.spoints <- function(x, y, col=par("col"), bg=NA, pch=par("pch"), cex=1,
                          col.smooth="red", span=2/3, iter=3,
                          col.point="#00000022", smooth.scatter=FALSE, ...) {
  ## points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  if (smooth.scatter) {
    smoothScatter(x, y, nrpoints = 0, add=TRUE)
  } else {
    points(x, y, pch=16, col=col.point, bg=bg, cex=cex)
  }

  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) {
    suppressWarnings({
      lines(stats::lowess(x[ok], y[ok], f = span, iter = iter),
            col = col.smooth, ...)
    })
  }
}
