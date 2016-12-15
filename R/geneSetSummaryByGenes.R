## The geneSetSummaryByGenes function (and its utility helper funcionts) are
## too hairy to separate, so putting them here. I especially don't intend the
## rename.feature.columns function, in particular, to be reused.
##
## Although the code looks hairy, note that there are unit tests in
## test-geneSetSummaryByGenes.R

##' @rdname geneSetSummaryByGenes
setMethod("geneSetSummaryByGenes", c(x="GeneSetDb"),
function(x, features, with.features=TRUE, feature.rename=NULL, ...,
         .external=TRUE) {
  stopifnot(is.character(features))
  unk.f <- setdiff(features, featureIds(x))
  if (length(unk.f)) {
    warning(length(unk.f), "/", length(features), " do not exist in GeneSetDb")
    features <- setdiff(features, unk.f)
  }
  if (length(features) == 0) {
    stop("No selected features in GeneSetDb")
  }
  x.sub <- subsetByFeatures(x, features)
  x.db <- x.sub@db[featureId %in% features]
  x.gs <- geneSets(x.sub, .external=FALSE)

  gs.cols <- c('collection', 'name', 'active')
  if (is.conformed(x)) {
    out <- x.gs[, c(gs.cols, 'n'), with=FALSE]
    setnames(out, 'n', 'N')
  } else {
    out <- x.gs[, c(gs.cols, 'N'), with=FALSE]
  }

  ## Each geneset row will be annotated with the number of features it
  ## it has, even if caller doesn't ask for `with.features`. To do so we first
  ## create geneset <-> featureId contingency table.
  ## we turn x.dt$featureId into a factor with levels == features because there
  ## maybe be some features that don't appear anywhere (due to conformation(?))

  x.dt <- dcast.data.table(x.db[, present := TRUE],
                           collection + name ~ featureId,
                           value.var='present', fill=FALSE)

  ## I was once given a strange GeneSetDb object with genesets in gdb@db that
  ## were not in gdb@table. Such an object would fail validObject(gdb), so this
  ## should never be. Let's dance around it, but also warn when it happens.
  nfids <- intersect(features, tail(names(x.dt), -2))
  if (!setequal(nfids, features)) {
    warning("You are in a scenario you thought should be impossible",
            immediate.=TRUE)
    features <- nfids
  }
  x.dt$n <- rowSums(as.matrix(x.dt[, -(1:2), with=FALSE]))
  setcolorder(x.dt, c(setdiff(names(x.dt), features), features))
  if (with.features) {
    x.dt <- rename.feature.columns(x.dt, x.sub, feature.rename)
  } else {
    ## If the caller didn't want the features, we just return the n
    x.dt <- x.dt[, list(collection, name, n)]
  }
  out <- out[x.dt, nomatch=0]
  ret.df(out, .external=.external)
})

##' @rdname geneSetSummaryByGenes
setMethod("geneSetSummaryByGenes", c(x="MultiGSEAResult"),
function(x, features, with.features=TRUE, feature.rename=NULL,
         method=NULL, max.p=0.3, p.col=c('padj', 'padj.by.collection', 'pval'),
         ..., .external=TRUE) {
  stopifnot(is.character(features))
  unk.f <- setdiff(features, featureIds(x))
  if (length(unk.f)) {
    warning(length(unk.f), "/", length(features), " do not exist in GeneSetDb")
    features <- setdiff(features, unk.f)
  }
  if (length(features) == 0) {
    stop("No selected features in GeneSetDb")
  }
  if (is.character(method)) {
    p.col <- match.arg(p.col)
    method <- match.arg(method, resultNames(x))
    stopifnot(max.p >= 0 & max.p <= 1)
    r <- result(x, method, .external=FALSE)
  }

  gdb <- copy(geneSetDb(x))
  gs <- geneSets(x)

  res <- geneSetSummaryByGenes(gdb, features, with.features,
                               feature.rename=FALSE, .external=FALSE, ...)

  if (with.features) {
    ## Account for impossible scenario that I "dance around" in the function
    ## upstairs
    features <- intersect(features, names(res))

    fstart <- ncol(res) - length(features) + 1
    meta <- res[, 1:(fstart - 1), with=FALSE]
    fcols <- res[, fstart:ncol(res), with=FALSE]
  } else {
    meta <- res
  }

  ## add logFC and pvalues
  xref <- match(paste(res$collection, res$name), paste(gs$collection, gs$name))
  meta$logFC <- gs$mean.logFC.trim[xref]
  if (is.character(method)) {
    rxref <- match(paste(res$collection, res$name), paste(r$collection, r$name))
    meta$padj <- r$padj[rxref]
  }

  ## When this is called on a result object, instead of a TRUE/FALSE
  ## geneset <-> feature contingency table, we replace the TRUE's with the
  ## logFC of that feature
  if (with.features) {
    ## spank on the logFCs of the genes
    lfc <- setkeyv(copy(logFC(x, .external=FALSE)), 'featureId')
    for (fid in features) {
      replace.me <- fcols[[fid]]
      fcols[, (fid) := ifelse(replace.me, lfc[fid]$logFC, 0)]
    }

    if (is.character(feature.rename)) {
      feature.rename <- feature.rename[1]
      ## Where do we find the featureId <-> renamed xref? If we find this column
      ## in the logFC(x) table, then that trumps all.
      if (is.character(lfc[[feature.rename]])) {
        xref <- match(gdb@db$featureId, lfc$featureId)
        gdb@db[[feature.rename]] <- lfc[[feature.rename]][xref]
      }
    }

    fcols <- rename.feature.columns(fcols, gdb, feature.rename)
    res <- cbind(meta, fcols)
  } else {
    res <- meta
  }

  setkeyv(res, c('collection', 'name'))

  if (is.character(method)) {
    keep <- r[[p.col]] <= max.p
    keep <- r[keep][, list(collection, name)]
    res <- res[keep, nomatch=0]
  }

  ret.df(res, .external=.external)
})

# utility function. accepts a pivoted geneset <-> feature table, and renames
# the columns associated with features in gdb. if feature.names is NULL, then
# the column names are prefixed with 'featureId_' to avoid column names that
# are numbers.
#
# columns in `out` are identified as features by virtue of them matching values
# returned from `featureids(gdb)`
#
# when an internal function calls this with feature.rename=FALSE, then no
# rename will happen at all.
rename.feature.columns <- function(x, gdb, feature.rename=NULL) {
  stopifnot(is.data.table(x))
  stopifnot(is(gdb, "GeneSetDb"))
  if (is.logical(feature.rename[1L]) && !feature.rename[1L]) {
    return(x)
  }
  rename.cols <- colnames(x)[colnames(x) %in% featureIds(gdb)]
  if (length(rename.cols) == 0) {
    return(out)
  }

  if (!is.null(feature.rename)) {
    feature.rename <- feature.rename[1L]
    if (is.null(gdb@db[[feature.rename]])) {
      warning("Can't rename features by: ", feature.rename,
              " (not found as a column in gdb@db")
      feature.rename <- NULL
    }
  }

  default.names <- setNames(paste0('featureId_', rename.cols), rename.cols)

  if (is.null(feature.rename)) {
    new.names <- default.names
  } else {
    db <- gdb@db[featureId %in% rename.cols]
    db <- unique(db, by=c('featureId'))
    new.names <- setNames(db[[feature.rename]], db[['featureId']])
    isna <- is.na(new.names)
    new.names[isna] <- default.names[isna]
  }
  setnames(x, names(new.names), new.names)
  setcolorder(x, c(setdiff(names(x), new.names), sort(new.names)))
}
