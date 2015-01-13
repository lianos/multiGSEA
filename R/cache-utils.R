design.params.name <- function(design, contrast, with.space=FALSE) {
  if (length(contrast) == 1) {
    if (is.character(contrast)) {
      design.cols <- which(colnames(design) == contrast)
    } else {
      design.cols <- contrast
    }
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

cache.data.fn <- function(method, design, contrast, extra.args=character(),
                          caller.env=sys.frame(sys.nframe() - 1L),
                          outdir=NULL, ext=NULL) {
  force(caller.env)
  if (missing(design)) {
    design <- matrix(nrow=0,ncol=1, dimnames=list(NULL, 'nodesign'))
  }
  if (missing(contrast)) {
    contrast <- ncol(design)
  }
  fn <- paste(method,
              sprintf('design(%s)', design.params.name(design, contrast)),
              sep='-')
  if (length(extra.args)) {
    extra <- sapply(extra.args, function(arg) {
      val <- get(arg,  envir=caller.env, inherits=FALSE)
      sprintf("%s-%s", arg, as.character(val))
    }, simplify=FALSE)
    extra <- paste(extra, collapse=',')
    fn <- paste(fn, extra, sep=',')
  }
  if (is.character(outdir) && isTRUE(file.exists(outdir))) {
    fn <- file.path(outdir, 'cache', fn)
  }
  if (!is.null(ext)) {
    ext <- sub('\\.+', '', ext)
    fn <- paste0(fn, '.', ext)
  }
  fn
}
