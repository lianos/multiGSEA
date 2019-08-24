validate.x.enrichtest <- validate.X
vaidate.inputs.enrichtest <- function(x, design, contrast, feature.bias,
                                      xmeta. = NULL, ...) {
  if (!is.data.frame(xmeta.)) {
    default <- .validate.inputs.full.design(x, design, contrast)
    if (length(default)) {
      return(default)
    }
  }

  ## Ensure that caller provides a named feature.bias vector
  errs <- list()
  if (missing(feature.bias) || is.null(feature.bias)) {
    # This is actually OK, a normal enrichment will be run w/ no bias
  } else {
    if (is.character(feature.bias)) {
      if (!is.data.frame(xmeta.) || !is.numeric(xmeta.[[feature.bias]])) {
        errs <- c(
          errs,
          paste("when feature.bias is a string, xmeta. needs to be a data.frame ",
                "with a numeric column named `feature.bias`"))
        return(errs)
      }
      feature.bias <- setNames(xmeta.[[feature.bias]], xmeta.[["featureId"]])
    }

    if (!is.numeric(feature.bias)) {
      errs <- 'feature.bias must be a numeric vector'
    }
    if (!all(rownames(x) %in% names(feature.bias))) {
      errs <- c(errs, 'some rownames(x) not in names(feature.bias)')
    }
  }

  return(errs)
}


#' This is a generic wrapper around limma::kegga to perform "biased enrichment"
#'
#' This, in principle, works similarly to goseq but uses [limma::kega()] as its
#' engine. If you don't want to calculate any type of biased enrichment, then
#' explicitly set feature.bias and prior.prob to `NULL`.
#'
#' Running goseq would sometimes throw errors in the `makespline` call from
#' [goseq::nullp()], so I jumped over to this given this insight:
#' https://support.bioconductor.org/p/65789/#65914
#'
#' @references
#' Young, M. D., Wakefield, M. J., Smyth, G. K., Oshlack, A. (2010).
#' Gene ontology analysis for RNA-seq: accounting for selection bias.
#' *Genome Biology* 11, R14. http://genomebiology.com/2010/11/2/R14
#'
#' @param feature.bias we will try to extract the average expression of the
#'   gene as teh default bias, but you can send in gene length, or
#'   what-have-you
do.enrichtest <- function(gsd, x, design, contrast = ncol(design),
                          feature.bias = "AveExpr", prior.prob = NULL,
                          restrict.universe = FALSE,
                          split.updown = TRUE, use.treat = FALSE,
                          feature.min.logFC = if (use.treat) log2(1.25) else 1,
                          feature.max.padj = 0.10, logFC = NULL, ...) {
  # 1. Specify up and down genes as a list of identifiers:
  #      list(Up = sigup, Down = sigdown)
  # 2a. Process feature.bias parameter so that it is a numeric bias vector
  # 2b. Pass feature.bias into .calc_prior_prob
  # 2c. Pass 2b into kegga.default(covariate = NULL, prior.prob = 2b)
  #
  # 3. if there are no degenes, do not run anything and set outgoing results
  #    with pvals hammered to 1.

  stopifnot(is.conformed(gsd, x))
  # stop("testing graceful method failure in multiGSEA call")
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
  if (is.null(logFC[["significant"]])) {
    logFC[, significant := {
      padj <= feature.max.padj & abs(logFC) >= feature.min.logFC
    }]
  }

  do <- c('all', if (split.updown) c('up', 'down') else NULL)
  if (any(c("up", "down") %in% do && is.numeric(logFC[["logFC"]]))) {
    if (!is.character(logFC[["direction"]])) {
      logFC[["direction"]] <- ifelse(logFC[["logFC"]] > 0, "up", "down")
    }
  }

  if (!(is.character(logFC[["direction"]]) && split.updown)) {
    split.direction <- NULL
  } else {
    split.direction <- "direction"
  }

  res <- enrichtest(gsd, logFC, selected = "significant",
                    groups = split.direction,
                    feature.bias = feature.bias,
                    restrict.universe = restrict.universe)

  # TODO: plit into all/up/down results like do.goseq
  #
  if (!is.numeric(prior.prob) && !is.null(feature.bias)) {
    if (is.character(feature.bias)) {
      fb <- logFC[[feature.bias]]
    }
  }
}

#' Runs biased enrichment test over a data.frame of features
#'
#' This function wraps [limma::kegga()] to perform biased enrichment tests over
#' gene set collection stored in a GeneSetDb (`gsd`) object. Its easiest to
#' use this function when the biases and selection criteria are stored as
#' columns of the input data.frame `dat`.
#'
#' @export
#' @importFrom limma kegga
#'
#' @param gsd The GeneSetDb
#' @param dat A data.frame with feature-level statistics. Minimally, this should
#'   have a `"featureId"` (character) column, but read on ...
#' @param selected Either the name of a logical column in `dat` used to subset
#'   out the features to run the enrichement over, or a character vector of
#'   `"featureId"`s that are selected from `dat[["featureId"]]`.
#' @param groups Encodes groups of features that we can use to test selected
#'   features individual, as well as "all" together. This can be specified by:
#'   (1) specifying a name of a column in `dat` to split the enriched features
#'   into subgroups. (2) A named list of features to intersect with `selected`.
#'   By default this is `NULL`, so we only run enrichment over
#'   all elements in `selected`. See examples for details.
#' @param feature.bias If `NULL` (default), no bias is used in enrichment
#'   analysis. Otherwise, can be the name of a column in `dat` to extract
#'   a numeric bias vector (gene length, GC content, average expression, etc.)
#'   or a named (using featureIds) numeric vector of the same.
#' @param universe Defaults to all elements in `dat[["featureId"]]`.
#' @param restrict.universe See same parameter in [limma::kegga()]
#' @param plot.bias See `plot` parameter in [limma::kega()]
#' @return A data.frame of pathway enrichment. The last N colums are enrichment
#'   statistics per pathway, grouped by the `groups` parameter. `P.all` are the
#'   stats for all selected features, and the remaingin `P.*` columns are for
#'   the features specifed by `groups`.
#' @examples
#' dgestats <- exampleDgeResult("human", "ensembl")
#' gdb <- getMSigGeneSetDb("h", "human", "ensembl")
#'
#' # Run enrichmnent without accounting for any bias
#' nobias <- enrichtest(gdb, dgestats, selected = "selected",
#'                      groups = "direction",
#'                      feature.bias = NULL)
#'
#' # Run enrichment and account for gene length
#' lbias <- enrichtest(gdb, dgestats, selected = "selected",
#'                     feature.bias = "effective_length")
#'
#' # plot length bias with DGE status
#' plot_enrichtest_bias(dgestats, "selected", "effective_length")
#'
#' # induce length bias and see what is the what ...............................
#' biased <- dgestats[order(dgestats$pval),]
#' biased$effective_length <- sort(biased$effective_length, decreasing = TRUE)
#' plot_enrichtest_bias(biased, "selected", "effective_length")
#' etest <- enrichtest(gdb, biased, selected = "selected",
#'                     groups = "direction",
#'                     feature.bias = "effective_length")
enrichtest <- function(gsd, dat, selected = "significant",
                       groups = NULL,
                       feature.bias = NULL, universe = NULL,
                       restrict.universe = FALSE,
                       plot.bias = FALSE, ...) {
  dat <- validate.xmeta(dat) # enforse featureId column
  if (is.null(universe)) universe <- dat[["featureId"]]
  gsd <- conform(gsd, universe, ...)

  if (test_string(selected) && test_logical(dat[[selected]])) {
    selected. <- dat[["featureId"]][dat[[selected]]]
  } else if (test_character(selected, min.len = 1L)) {
    selected. <- intersect(selected, universe)
    if (!setequal(selected., selected)) {
      warning("Only ", length(selected.), " / ", length(selected),
              "features found in 'dat'")
    }
  }

  if (!is.character(selected.)) {
    stop("Illegal argument type of `selected`: ", class(selected)[1L])
  }

  de <- list(all = selected.)

  if (!is.null(groups)) {
    if (test_string(groups) && test_character(dat[[groups]])) {
      groups <- split(dat[["featureId"]], dat[[groups]])
    }
  }

  if (is.list(groups)) {
    if (is.null(names(groups))) {
      names(groups) <- paste0("group_", seq_along(groups))
    }
    for (group in names(groups)) {
      features <- intersect(selected., groups[[group]])
      if (length(features)) {
        if (group == "all") group <- "all2"
        de[[group]] <- features
      }
    }
  }

  if (!is.null(feature.bias)) {
    if (test_string(feature.bias) && test_numeric(dat[[feature.bias]])) {
      feature.bias <- setNames(dat[[feature.bias]], dat[["featureId"]])
    }
    if (!is.numeric(feature.bias) && all(universe %in% names(feature.bias))) {
      warning("feature.bias vector does not cover universe: running unbiased ",
              "enrichment")
      feature.bias <- NULL
    }
  }

  if (is.numeric(feature.bias)) {
    feature.bias <- feature.bias[universe]
    if (!requireNamespace("BiasedUrn", quietly = TRUE)) {
      warning("BiasedUrn package not installed, running in unbiased mode ...",
              immediate. = TRUE)
      feature.bias <- NULL
    }
  }

  # Transform GeneSetDb into required kegga bits ...............................

  # gene.pathway is a 2d data.frame like so:
  #    GeneID      PathwayID
  #     10327  path:hsa00010
  #       124  path:hsa00010
  #       125  path:hsa00010
  #       126  path:hsa00010
  gene.pathway <- local({
    gp <- as.data.frame(gsd, active.only = TRUE)
    data.frame(GeneID = gp[["featureId"]], PathwayID = encode_gskey(gp),
               stringsAsFactors = FALSE)
  })

  # pathway.names is a 2col data.frame, like so:
  #        PathwayID                                  Description
  #    path:hsa00010             Glycolysis / Gluconeogenesis ...
  #    path:hsa00020                Citrate cycle (TCA cycle) ...
  #    path:hsa00030                Pentose phosphate pathway ...
  pathway.names <- local({
    gs <- geneSets(gsd, active.only = TRUE)
    gs$key <- encode_gskey(gs)
    data.frame(PathwayID = gs$key, Description = gs$key,
               stringsAsFactors = FALSE)
  })

  kres <- NULL
  if (length(de[["all"]])) {
    kres <- limma::kegga(de, universe = universe, covariate = feature.bias,
                         restrict.universe = restrict.universe,
                         gene.pathway = gene.pathway,
                         pathway.names = pathway.names,
                         plot = plot.bias && !is.null(feature.bias))
  }

  setattr(kres, "rawresult", TRUE)
  kres
}

#' Plots bias of coviarate to DE / selected status
#'
#' This meat and potatoes of this function code was extracted from limma::kegga,
#' originally written by Gordon Smyth and Yifang Hu.
#'
#' @export
#' @importFrom stats approx
#' @importFrom limma barcoddeplot
plot_enrichtest_bias <- function(x, selected, feature.bias,
                                 title = "DE status vs bias") {
  assert_multi_class(x, c("data.frame", "tibble"))
  if (test_string(selected)) {
    selected <- x[[selected]]
  }
  assert_logical(selected)
  if (test_string(feature.bias)) {
    feature.bias <- x[[feature.bias]]
  }
  assert_numeric(feature.bias, len = length(selected))

  span <- approx(x = c(20,200), y = c(1, 0.5),
                 xout = sum(selected),
                 rule = 2, ties = list("ordered", mean))$y
  limma::barcodeplot(feature.bias,
                     index = selected,
                     worm = TRUE, span.worm = span,
                     main = "DE status vs covariate (manual)")

}
