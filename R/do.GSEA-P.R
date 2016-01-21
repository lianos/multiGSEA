##' @include validateInputs.R
NULL

validate.inputs.gseap <- .validate.inputs.logFC.only
validate.x.gseap <- validate.X

##' Runs a method analagous to GSEA-P
##'
##' See the following link for reference
##' http://guangchuangyu.github.io/2015/11/comparison-of-clusterprofiler-and-gsea-p/
##'
##' @importFrom clusterProfiler GSEA
do.gseap <- function(gsd, x, design, contrast,
                     logFC=NULL, score.by=c('t', 'logFC'),
                     robust.fit=FALSE, robust.eBayes=FALSE, ...) {
  score.by <- match.arg(score.by)
  stopifnot(is.conformed(gsd, x))
  if (!missing(design) && missing(contrast)) {
    contrast <- ncol(design)
  }

  if (ncol(x) > 1) {
    if (is.null(logFC)) {
      logFC <- calculateIndividualLogFC(x, design, contrast, robust.fit,
                                        robust.eBayes, ..., .external=FALSE)
    } else {
      is.logFC.like(logFC, x, as.error=TRUE)
    }
    stats <- setNames(logFC[[score.by]], logFC$featureId)
  } else {
    ## This already a column matrix of precomputed things (logFC, perhaps)
    ## to rank
    stats <- setNames(x[, 1L], rownames(x))
  }

  args <- list(...)
  call.args <- as.list(formals(clusterProfiler::GSEA))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  stats <- stats[rownames(x)]

  ## Double checking that feature <-> expression IDs match up here. This is
  ## supposed to be redundant, but measure twice & cut once.
  if (any(is.na(stats))
      || length(stats) != nrow(x)
      || any(rownames(x) != names(stats))) {
    stop("individual stats do not match up with rownames of x")
  }

  gs <- geneSets(gsd, .external=FALSE)
  gsd.df <- as.data.frame(gsd)
  gsd.df$ont <- paste(gsd.df$collection, gsd.df$name, sep=';;')
  gsd.df$gene <- gsd.df$featureId

  call.args[['geneList']] <- stats
  call.args[['TERM2GENE']] <- gsd.df[, c('ont', 'gene')]
  call.args[['pvalueCutoff']] <- 1
  call.args[['minGSSize']] <- min(gs$n)

  res <- do.call(clusterProfiler::GSEA, call.args)
  out <- gs[, list(collection, name)]
  ## TODO: Finish do.gseap function

  ## out[, pval := pvals]
  ## out[, padj := p.adjust(pval, 'BH')]
}
