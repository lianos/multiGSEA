##' @include validateInputs.R
NULL

validate.inputs.goana <- .validate.inputs.full.design

do.goana <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
                     use.cache=TRUE, ...) {
  stopifnot(is.conformed(gsd, x))
  stop("limma::goana not yet implemented")
}


