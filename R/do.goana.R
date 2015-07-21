##' @include validateInputs.R
NULL

validate.inputs.goana <- .validate.inputs.full.design
validate.x.goana <- validate.X

do.goana <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
                     use.cache=TRUE, ...) {
  stop("not yet implemented")
  stopifnot(is.conformed(gsd, x))
}


