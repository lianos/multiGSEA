##' @include validateInputs.R
NULL

validate.inputs.camera <- .validate.inputs.full.design
validate.x.camera <- validate.X

##' Worker function to run camera from within a multiGSEA pipeline
##'
##' @description
##' camera as originally implemented tends to be very conservative, this is
##' due to the intra-geneset correlation that it originally estimates from
##' the data. It has been argued that this is hard to estimate and becomes
##' very conservative for experiments with small N.
##'
##' As a result, an \code{inter.gene.cor} parameter was introduced to camera
##' in limma v3.24.14, and as of Bioc3.3, its value is set to 0.01 by default
##' instead of letting camera approximate it from the data.
##'
##' For the original discussion around this topic, please refer to the following
##' thread in the bioconductor support forum:
##'
##' \url{https://support.bioconductor.org/p/70005/#70195}
##'
##' \strong{This function is not meant to be called directly, it should only be
##' called internally within \code{multiGSEA}}
do.camera <- function(gsd, x, design, contrast=ncol(design),
                      gs.idxs=as.list(gsd, active.only=TRUE, value='x.idx'),
                      ...) {
  stopifnot(is.conformed(gsd, x))

  args <- list(...)
  if (!is.null(args$inter.gene.cor)) {
    ## Preset inter.gene.cor values were implemented in limma v3.24.14
    if (packageVersion('limma') < '3.24.14') {
      warning("inter.gene.cor values for camera require limma >= 3.24.14, ",
              "the parameter is ignored.", immediate.=TRUE)
    }
  }

  call.args <- as.list(formals(limma::camera.default))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  call.args[['y']] <- x
  call.args[['index']] <- gs.idxs
  call.args[['design']] <- design
  call.args[['contrast']] <- contrast
  call.args[['sort']] <- FALSE
  call.args[['...']] <- NULL

  res <- do.call(camera, call.args)
  setattr(res, 'rawresult', TRUE)
}

mgres.camera <- function(res, gsd, ...) {
  if (!isTRUE(attr(res, 'rawresult'))) return(res)
  out <- cbind(
    geneSets(gsd, as.dt=TRUE)[, list(collection, name)],
    as.data.table(res))
  NGenes <- NULL # silence R CMD check NOTEs
  out[, NGenes := NULL]

  # camera result doesn't have an FDR column if we only tested on geneset
  # https://github.com/lianos/multiGSEA/issues/7
  if (is.null(out[["FDR"]])) out[, FDR := p.adjust(PValue)]

  setnames(out, c('PValue', 'FDR'), c('pval', 'padj'))
}
