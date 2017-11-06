## Note: I won't implement limma::goana since they do not accept their own
##       universe. The methodology in goseq is equivalent except for minor
##       detail when low number of DE genes (cf ?goana)
##' @include validateInputs.R
NULL

validate.x.goseq <- validate.X
validate.inputs.goseq <- function(x, design, contrast, feature.bias, ...) {
  default <- .validate.inputs.full.design(x, design, contrast)
  if (length(default)) {
    return(default)
  }
  ## Ensure that caller provides a named feature.bias vector
  errs <- list()
  if (missing(feature.bias)) {
    errs <- paste('feature.bias vector is required, use "hyperGeometricTest"',
                  'if you do not have one')
    return(errs)
  }
  if (!is.numeric(feature.bias)) {
    errs <- 'feature.bias must be a numeric vector'
  }
  if (!all(rownames(x) %in% names(feature.bias))) {
    errs <- c(errs, 'some rownames(x) not in names(feature.bias)')
  }
  return(errs)
}

##' Performs goseq analysis significance of gene set membership.
##'
##' Genes are selected for testing against each geneset by virture of them
##' passing a maximum FDR and minimum log fold change as perscribed by the
##' \code{min.logFC} and \code{max.padj} parameters, respectfully.
##'
##' Note that we are intentionally adding a hyperG.selected column by reference
##' so that this information is kicked back to the caller multiGSEA function
##' and included in downstream reporting.
##'
##' @param gsd The \code{\link{GeneSetDb}} for analysis
##' @param x The expression object
##' @param design Experimental design
##' @param contrast The contrast to test
##' @param feature.bias a named vector as long as \code{nrow(x)} that has the
##'   "bias" information for the features/genes tested (ie. vector of gene
##'   lengths). \code{names(feature.bias)} should equal \code{rownames(x)}.
##'   The caller MUST provide this. The goseq package provides a
##'   \code{\link[goseq]{getlength}} function which facilitates getting default
##'   values for these if you do not have the correct values used in your
##'   analysis. If there is no way for you to get this information, then use
##'   \code{method='hyperGeometricTest'}.
##' @param method The method to use to calculate the unbiased category
##'   enrichment scores
##' @param repcnt Number of random samples to be calculated when random sampling
##'   is used. Ignored unless \code{method="Sampling"}.
##' @param use_genes_without_cat A boolean to indicate whether genes without a
##'   categorie should still be used. For example, a large number of gene may
##'   have no GO term annotated. If this option is set to FALSE, those genes
##'   will be ignored in the calculation of p-values (default behaviour). If
##'   this option is set to TRUE, then these genes will count towards the total
##'   number of genes outside the category being tested.
##' @param direction Same as direction in \code{GOstats}
##' @param plot.fit To plot (or not) the bias in selected genes vs.
##'   \code{feature.bias}.
##' @param logFC The logFC data.table from \code{calculateIndividualLogFC}
##' @param ... arguments to pass down into \code{calculateIndividualLogFC}
##' @return A data.table of goseq results. The "pval" column here refers to
##'   pval.over, for simplicity in other places. If \code{split.updown=TRUE},
##'   a list of data.table's are returned named 'goseq', 'goseq.up', and
##'   'goseq.down' which are the results of running goseq three independent
##'   times.
do.goseq <- function(gsd, x, design, contrast=ncol(design),
                     feature.bias,
                     method="Wallenius",
                     repcnt=2000, use_genes_without_cat=TRUE,
                     split.updown=TRUE,
                     direction=c('over', 'under'),
                     plot.fit=FALSE, use.treat=FALSE,
                     feature.min.logFC=if (use.treat) log2(1.25) else 1,
                     feature.max.padj=0.10, logFC=NULL, ...) {
  stopifnot(is.conformed(gsd, x))
  direction <- match.arg(direction)

  if (is.null(logFC)) {
    treat.lfc <- if (use.treat) feature.min.logFC else NULL
    logFC <- calculateIndividualLogFC(x, design, contrast, treat.lfc=treat.lfc,
                                      ..., as.dt=TRUE)
    if (use.treat) {
      logFC[, significant := padj <= feature.max.padj]
    }
  }
  is.logFC.like(logFC, x, as.error=TRUE)
  logFC <- setDT(copy(logFC))
  if (is.null(logFC$significant)) {
    logFC[, significant := {
      padj <= feature.max.padj & abs(logFC) >= feature.min.logFC
    }]
  }

  do <- c('all', if (split.updown) c('up', 'down') else NULL)

  out <- sapply(do, function(dge.dir) {
    drawn <- switch(dge.dir,
                    all=logFC[significant == TRUE]$featureId,
                    up=logFC[significant == TRUE & logFC > 0]$featureId,
                    down=logFC[significant == TRUE & logFC < 0]$featureId)
    res <- suppressWarnings({
      multiGSEA::goseq(gsd, drawn, rownames(x), feature.bias, method,
                       repcnt, use_genes_without_cat, plot.fit=plot.fit,
                       do.conform=FALSE, as.dt=FALSE, .pipelined=TRUE)
    })
  }, simplify=FALSE)
  if (length(out) == 1L) {
    out <- out[[1L]]
  } else {
    setattr(out, 'mgunlist', TRUE)
  }
  out
}

##' Perform goseq Enrichment tests across a GeneSetDb.
##'
##' Note that we do not import things from goseq directly, and only load
##' it if this function is fired. I can't figure out a way to selectively
##' import functions from the goseq package without it having to load its
##' dependencies, which take a long time -- and I don't want loading multiGSEA
##' to take a long time. So, the goseq package has moved to Suggests and then
##' is loaded within this function when necessary.
##'
##' @export
##'
##' @param gsd The \code{GeneSetDb} object to run tests against
##' @param selected The ids of the selected features
##' @param universe The ids of the universe
##' @param feature.bias a named vector as long as \code{nrow(x)} that has the
##'   "bias" information for the features/genes tested (ie. vector of gene
##'   lengths). \code{names(feature.bias)} should equal \code{rownames(x)}.
##'   If this is not provided, all feature lengths are set to 1 (no bias).
##'   The goseq package provides a \code{\link[goseq]{getlength}} function which
##'   facilitates getting default values for these if you do not have the
##'   correct values used in your analysis.
##' @param method The method to use to calculate the unbiased category
##'   enrichment scores
##' @param repcnt Number of random samples to be calculated when random sampling
##'   is used. Ignored unless \code{method="Sampling"}.
##' @param use_genes_without_cat A boolean to indicate whether genes without a
##'   categorie should still be used. For example, a large number of gene may
##'   have no GO term annotated. If this option is set to FALSE, those genes
##'   will be ignored in the calculation of p-values (default behaviour). If
##'   this option is set to TRUE, then these genes will count towards the total
##'   number of genes outside the category being tested.
##' @param do.conform By default \code{TRUE}: does some gymnastics to conform
##'   the \code{gsd} to the \code{universe} vector. This should neber be set
##'   to \code{FALSE}, but this parameter is here so that when this function
##'   is called from the \code{\link{multiGSEA}} codepath, we do not have to
##'   reconform the \code{GeneSetDb} object, because it has already been done.
##' @param active.only If \code{TRUE}, only "active" genesets are used
##' @param value The featureId types to extract from \code{gsd}
##' @param .pipelined If this is being external to a multiGSEA pipeline, then
##'   some additional cleanup of columns name output will be done. Otherwise
##'   the column renaming and post processing is left to the do.goseq caller
##'   (Default: \code{FALSE}).
##' @return A \code{data.table} of results, similar to goseq output. The output
##'   from \code{\link[goseq]{nullp}} is added to the outgoing data.table as
##'   an attribue named \code{"pwf"}.
goseq <- function(gsd, selected, universe, feature.bias,
                  method=c("Wallenius", "Sampling", "Hypergeometric"),
                  repcnt=2000, use_genes_without_cat=TRUE,
                  plot.fit=TRUE, do.conform=TRUE, as.dt=FALSE,
                  .pipelined=FALSE) {
  gseq <- tryCatch(loadNamespace("goseq"), error=function(e) NULL)
  if (is.null(goseq)) {
    stop("You must install the goseq package to enable this functionality")
  }
  stopifnot(is(gsd, 'GeneSetDb'))
  stopifnot(is.character(selected))
  stopifnot(is.character(universe))
  selected <- unique(selected)
  universe <- unique(universe)
  method <- match.arg(method)
  stopifnot(all(selected %in% universe))
  stopifnot(is.numeric(feature.bias) && all(universe %in% names(feature.bias)))
  feature.bias <- feature.bias[universe]

  ## This needs to be conformed to work
  if (!is.conformed(gsd) || !do.conform) {
    gsd <- conform(gsd, universe)
  }
  stopifnot(is.conformed(gsd, universe))

  # silence R CMD check NOTEs
  category <- padj_over <- pval_over <- padj_under <- pval_under <-  NULL

  gs <- geneSets(gsd, active.only=TRUE, as.dt=TRUE)
  gs[, category := encode_gskey(gs)]
  g2c <- as.data.frame(gsd, active.only=TRUE, value='x.id')
  g2c <- transform(g2c, category=encode_gskey(g2c))
  g2c <- g2c[, c('category', 'featureId')]

  selected <- intersect(selected, universe)
  de.genes <- setNames(integer(length(universe)), universe)
  de.genes[selected] <- 1L

  pwf <- gseq$nullp(de.genes, bias.data=feature.bias, plot.fit=plot.fit)
  res <- gseq$goseq(pwf, gene2cat=g2c, method=method, repcnt=repcnt,
                    use_genes_without_cat=use_genes_without_cat)

  ## resort output to match gs ordering
  res <- res[match(gs$category, res$category),,drop=FALSE]
  res <- if (as.dt) setDT(res) else setDF(res)
  setattr(res, 'pwf', pwf)
  setattr(res, 'rawresult', TRUE)
  res
}

mgres.goseq <- function(res, gsd, ...) {
  if (!isTRUE(attr(res, 'rawresult'))) return(res)
  stopifnot(is.data.frame(res), is(gsd, "GeneSetDb"))
  gs <- geneSets(gsd, active.only=TRUE, as.dt=TRUE)
  category <- NULL ## silence "no visible binding" note
  gs[, category := encode_gskey(gs)]
  xref <- match(gs$category, res$category)
  if (any(is.na(xref))) {
    stop("The conformed GeneSetDb used in analysis is not this one")
  }
  res <- res[xref,,drop=FALSE]
  out <- gs[, list(collection, name, n, n.drawn=res$numDEInCat)]
  rcols <- setdiff(names(res), c('category', 'numInCat', 'numDEInCat'))
  for (rcol in rcols) {
    out[, (rcol) := res[[rcol]]]
  }
  setnames(out, c('over_represented_pvalue', 'under_represented_pvalue'),
           c('pval', 'pval.under'))
  padj.under <- pval.under <- NULL # silence R CMD check NOTEs
  out[, padj := p.adjust(pval, 'BH')]
  out[, padj.under := p.adjust(pval.under, 'BH')]
}
