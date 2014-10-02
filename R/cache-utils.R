design.params.name <- function(design, contrast, with.space=FALSE) {
  if (length(contrast) == 1) {
    design.cols <- contrast
    coef.vals <- ''
  } else {
    design.cols <- which(contrast != 0)
    coef.vals <- contrast[design.cols]
    ## Ensure that both positive and negative coefficients are printed with
    ## their sign
    coef.vals <- paste(
      ifelse(coef.vals > 0, '+', '-'),
      abs(coef.vals), sep='')
  }
  coef.names <- colnames(design)[design.cols]
  ## Ensure that positive coef values are printed with `+`
  collapse <- if (with.space) ' ' else ''
  paste(coef.vals, coef.names, sep='', collapse=collapse)
}

cache.data.fn <- function(method, design, contrast, extra.args=list(),
                          outdir=NULL, ext=NULL) {
  fn <- paste(method,
              sprintf('design(%s)', design.params.name(design, contrast)),
              sep='-')
  if (length(extra.args)) {
    extra <- paste(names(extra.args), as.character(extra.args),
                   sep='-', collapse=',')
    fn <- paste(fn, extra, sep=',')
  }
  if (!is.null(outdir)) {
    fn <- file.path(outdir, 'cache', fn)
  }
  if (!is.null(ext)) {
    ext <- sub('\\.+', '', ext)
    fn <- paste0(fn, '.', ext)
  }
  fn
}
