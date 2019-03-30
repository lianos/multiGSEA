validate.inputs.cameraPR <- .validate.inputs.preranked
validate.x.cameraPR <- validate.X

#' @noRd
#' @importFrom limma cameraPR
do.cameraPR <- function(gsd, x, design, contrast=ncol(design),
                        score.by=c('t', 'logFC', 'pval'), logFC=NULL,
                        gs.idxs=as.list(gsd, active.only=TRUE, value='x.idx'),
                        ...) {
  score.by <- match.arg(score.by)
  stopifnot(is.conformed(gsd, x))

  stats <- extract_preranked_stats(x, design, contrast, score.by=score.by,
                                   logFC=logFC, ...)

  call.args <- as.list(formals(limma::cameraPR.default))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  call.args[['statistic']] <- stats
  call.args[['index']] <- gs.idxs
  call.args[['sort']] <- FALSE
  call.args[['...']] <- NULL
  res <- do.call(limma::cameraPR.default, call.args)
  setattr(res, 'rawresult', TRUE)
}

mgres.cameraPR <- mgres.camera
