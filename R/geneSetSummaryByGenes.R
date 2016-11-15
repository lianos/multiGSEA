## The geneSetSummaryByGenes function (and its utility helper funcionts) are
## too hairy to separate, so putting them here. I especially don't intend the
## rename.feature.columns function, in particular, to be reused.
##
## Although the code looks hairy, note that there are unit tests in
## test-geneSetSummaryByGenes.R
setMethod("geneSetSummaryByGenes", c(x="GeneSetDb"),
function(x, features, with.features=TRUE, feature.rename=NULL, ...,
         .external=TRUE) {
  stopifnot(is.character(features))
  unk.f <- setdiff(features, featureIds(x))
  if (length(unk.f)) {
    stop("These features do not exist in x: ", paste(unk.f, collapse=','))
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

  if (with.features) {
    ## Create geneset <-> featureId contingency table and spank it on to the
    ## end of the outoing data.table
    x.dt <- dcast.data.table(x.db[, present := TRUE],
                             collection + name ~ featureId,
                             value.var='present', fill=FALSE)
    x.dt <- rename.feature.columns(x.dt, x.sub, feature.rename)
    out <- out[x.dt]
  }

  ret.df(out, .external=.external)
})

setMethod("geneSetSummaryByGenes", c(x="MultiGSEAResult"),
function(x, features, with.features=TRUE, feature.rename=NULL,
         name=NULL, max.p=0.3, p.col=c('padj', 'padj.by.collection', 'pval'),
         ..., .external=TRUE) {
  if (is.character(name)) {
    p.col <- match.arg(p.col)
    name <- match.arg(name, resultNames(x))
    stopifnot(max.p >= 0 & max.p <= 1)
  }

  gdb <- geneSetDb(x)
  res <- geneSetSummaryByGenes(gdb, features, with.features,
                               feature.rename=FALSE, .external=FALSE, ...)
  ## When this is called on a result object, instead of a TRUE/FALSE
  ## geneset <-> feature contingency table, we replace the TRUE's with the
  ## logFC of that feature
  if (with.features) {
    ## spank on the logFCs of the genes
    lfc <- setkeyv(copy(logFC(x, .external=FALSE)), 'featureId')
    for (fid in features) {
      replace.me <- res[[fid]]
      res[, (fid) := ifelse(replace.me, lfc[fid]$logFC, 0)]
    }

    if (is.character(feature.rename)) {
      feature.rename <- feature.rename[1]
      ## Where do we find the featureId <-> renamed xref? If it's not a column
      ## in the gdb object, can we make it one from the logFC table?
      if (is.null(gdb@db[[feature.rename]]) && is.character(lfc[[feature.rename]])) {
        xref <- lfc[, c('featureId', feature.rename), with=FALSE]
        gdb@db[, (feature.rename) := {
          xref[[feature.rename]][match(featureId, xref$featureId)]
        }]
      }
    }

    res <- rename.feature.columns(res, gdb, feature.rename)
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
      warning("Can't rename features by: ", feature.rename)
      feature.rename <- NULL
    }
  }

  default.names <- setNames(paste0('featureId_', rename.cols), rename.cols)
  # browser()
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
