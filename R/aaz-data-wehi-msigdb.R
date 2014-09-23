##' WEHI MSigDB data for human and mouse
##'
##' This data is downloaded from the
##' \href{http://bioinf.wehi.edu.au/software/MSigDB/}{WEHI MSigDB for R page},
##' Which includes \code{*.RData} objects for human and mouse of the
##' \href{http://www.broadinstitute.org/gsea/msigdb/}{MSigDB database}.
##'
##' These datasets are versioned, and the last two datasets are kept available
##' inside this package (starting from v4.0, which was last modified on
##' 21 October 2013). For all analyses using these pre-generated genesets,
##' the default option is to use the latest version of this data.
##'
##' The \code{*.RData} objects are simply named lists of Entrez IDs. The names
##' of the elements of the list correspond to their gene set ID, ie. "chr5q23"
##' or "KEGG_GLYCOLYSIS_GLUCONEOGENESIS". Each data object holds only one
##' such list.
##'
##' @docType data
##' @name MSigDB
NULL

.wehi.msigdb.species <- c('human', 'mouse')
.wehi.msigdb.current <- 'v4.0'
.wehi.msigdb.datadir <- system.file('extdata', 'MSigDB', package='multiGSEA')

.available.msigdb.versions <- function() {
  dir(.wehi.msigdb.datadir)
}

.parse.msigdb.version <- function(v) {
  if (is.numeric(v)) {
    v <- sprintf('v%.1f', v)
  }
  stopifnot(is.character(v))
  valid <- v %in% .available.msigdb.versions()
  if (!valid) {
    msg <- "Illegal MSigDB version number (%s). Valid versions are:\n    %s"
    msg <- sprintf(msg, version,
                   paste(.available.msigdb.versions(), collapse="\n    "))
    stop(msg)
  }
  v
}

.parse.msigdb.ids <- function(id, species=.wehi.msigdb.species) {
  species <- match.arg(species)
  valid.ids <- paste0('c', if (species == 'human') 1:7 else 2:7)
  if (is.numeric(id)) {
    id <- paste0('c', id)
  }
  stopifnot(is.character(id))
  bad.ids <- setdiff(id, valid.ids)
  if (length(bad.ids)) {
    msg <- sprintf('Unknown MSigDB %s gene set IDs (%s)', species,
                   paste(bad.ids, collapse=','))
    stop(msg)
  }
  id
}

##' Fetches the requested MSigDB gene sets.
##'
##' @export
##'
##' @param id A character vector (or integer vector) indicating the desired
##' genesets
##' @param species Character vector indicating which species desired,
##' \code{"human"} or \code{"mouse"}
##' @param version The version of the MSigDB desired. Defaults to latest.
##'
##' @return A list of genesets, named by \code{id}.
getMSigDBset <- function(id, species=.wehi.msigdb.species,
                         version=.wehi.msigdb.current) {
  species <- match.arg(species)
  version <- .parse.msigdb.version(version)
  id <- .parse.msigdb.ids(id, species)
  data.dir <- file.path(.wehi.msigdb.datadir, version)
  rdata <- sapply(id, function(.id) {
    fn.regex <- sprintf('%s_%s_?.*\\.rdata', species, .id)
    fn <- dir(data.dir, fn.regex, ignore.case=TRUE, full.names=TRUE)
    if (length(fn) != 1) {
      msg <- "Filename regex (%s) didn't match just one file in data dir (%s)"
      stop(sprintf(msg, fn.regex, data.dir))
    }
    get(load(fn))
  }, simplify=FALSE)
  rdata
}
