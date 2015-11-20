.onLoad <- function(libname, pkgname) {
  ## Setup default option values
  opts <- options()

  pkg.opts <- list(
    multiGSEA.df.return='data.frame')
  toset <- !(names(pkg.opts) %in% names(opts))
  if (any(toset)) {
    options(pkg.opts[toset])
  }

  df.return <- getOption('multiGSEA.df.return')
  if (!df.return %in% c('data.table', 'data.frame')) {
    warning("Invalid value for options(multiGSEA.df.return) (",
            df.return, "), setting to 'data.table'", immediate.=TRUE)
  }

  invisible()
}
