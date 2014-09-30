design.params.name <- function(design, contrast) {
  if (length(contrast) == 1) {
    design.cols <- contrast
    coef.vals <- ''
  } else {
    design.cols <- which(contrast != 0)
    coef.vals <- contrast[design.cols]
    if (all(as.integer(coef.vals) == coef.vals)) {
      coef.vals <- as.character(coef.vals)
    } else {
      coef.vals <- sprintf("%.2f", coef.vals)
    }
  }
  coef.names <- colnames(design)[design.cols]
  paste(coef.vals, coef.names, sep='', collapse='')
}

cache.data.fn <- function(method, design, contrast, extra.args=list(),
                          outdir=NULL) {
  fn <- paste(method,
              sprintf('design(%s)', design.params.name(design, contrast)),
              sep='-')
  if (length(extra.args)) {
    extra <- paste(names(extra.args), as.character(extra.args),
                   sep='-', collapse=',')
    fn <- paste(fn, extra, sep=',')
  }
  if (!is.null(outdir)) {
    file.path(outdir, 'cache', fn)
  }
  fn
}
