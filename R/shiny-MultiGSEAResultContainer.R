##' A wrapper for a MultiGSEAResult that makes explicity some shiny bits
##'
##' @export
##' @param x A \code{\link[multiGSEA]{MultiGSEAResult}} object, or a path to
##'   an rds-serliazed one
##' @return a \code{MultiGSEAResultContainer} object
MultiGSEAResultContainer <- function(x) {
  if (is.character(x)) {
    ## Assume this is a file
    if (!file.exists(x)) {
      stop("file does not exist: ", x)
    }
    mg <- readRDS(x)
  } else if (is(x, 'MultiGSEAResult')) {
    mg <- x
  } else {
    mg <- NULL
  }

  if (!is(mg, 'MultiGSEAResult')) {
    stop("Don't know how to create container from: ", class(x)[1L])
  }

  methods <- local({
    tmp <- resultNames(mg)
    if (length(tmp) == 0L) {
      stop("No GSEA methods found in MultiGSEAResult")
    }
    ## I am biased and prefer to show these methods first, if available
    pref <- c('camera', 'goseq', 'goseq.up', 'goseq.down')
    pref <- pref[pref %in% tmp]
    c(pref, setdiff(tmp, pref))
  })

  gs.choices <- gs.select.choices(mg)

  out <- list(mg=mg, methods=methods, choices=gs.choices)
  class(out) <- c('MultiGSEAResultContainer', class(out))
  out
}
