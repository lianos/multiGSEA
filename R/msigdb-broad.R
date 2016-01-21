.msigdb.species <- c('human', 'mouse')
.msigdb.collections <- list(
  'v4.0'=paste0('c', 1:7),
  'v5.0'=c(paste0('c', 1:7), 'h'),
  'v5.1'=c(paste0('c', 1:7), 'h'))
.msigdb.version.current <- tail(names(.msigdb.collections), 1L)

##' Fetches a \code{GeneSetDb} from geneset collections defined in MSigDB.
##'
##' @description
##' This provides versioned genesets from gene set collections defined in
##' \href{http://software.broadinstitute.org/gsea/msigdb}{MSigDB}. The following
##' versions are included in this package:
##'
##' \itemize{
##'   \item v4.0
##'   \item v5.0
##'   \item v5.1
##' }
##'
##' Starting from version 5.1, the \code{GeneSetDb} includes the symbols for
##' the features in the genesets (in \code{gdb@db}) as well as additional meta
##' information about the genesets (in \code{geneSets(gdb)}) for "subcategory"
##' (which may be most useful for GO gene sets, eg. CC vs BP vs MF) and
##' "organism" which defines the organism in which the geneset was derived
##' from ("Homo sapiens" or "Mus musculus").
##'
##' @export
##' @rdname MSigDB
##' @param collection character vector specifying the collections you want
##'   (c1, c2, ..., c7, h)
##' @param species human or mouse?
##' @param version the version of the MSigDB database to use.
##' @return a \code{GeneSetDb} object
getMSigGeneSetDb <- function(collection, species='human',
                             version=.msigdb.version.current) {
  species <- resolve.species(species)
  version <- match.arg(version, names(.msigdb.collections))
  avail.cols <- .msigdb.collections[[version]]
  if (species == 'Mus_musculus') {
    avail.cols <- setdiff(avail.cols, 'c1')
  }
  if (is.null(collection) || missing(collection)) {
    collection <- avail.cols
  } else {
    bad.col <- setdiff(collection, avail.cols)
    if (length(bad.col)) {
      stop("Illegal MSigDB collection IDs in ", species, ":\n  ",
           paste(bad.col, collapse=','))
    }
  }

  if (version == 'v4.0') {
    warning("Retrieving MSigDB definitions from WEHI for v4.0", immediate.=TRUE)
    spec <- c('Mus_musculus'='mouse', 'Homo_sapiens'='human')[species]
    return(getMSigDBset.wehi(collection, spec, version=version))
  }

  fn <- sprintf('MSigDB.%s.GeneSetDb.rds', species)
  fn <- system.file('extdata', 'MSigDB', version, fn, package='multiGSEA')
  out <- readRDS(fn)

  if (!setequal(collection, avail.cols)) {
    gs <- geneSets(out, .external=FALSE)
    keep <- gs$collection %in% collection
    out <- subset.GeneSetDb(out, gs$collection %in% collection)
  }

  out
}

##' @rdname MSigDB
##' @export getMSigDBset
getMSigDBset <- getMSigGeneSetDb

resolve.species <- function(x) {
  stopifnot(is.character(x) && length(character) == 1L)
  xx <- tolower(x)
  opts <- c(human='Homo_sapiens',
            homo_sapiens='Homo_sapiens',
            hsa='Homo_sapiens',
            ## mouse
            mouse='Mus_musculus',
            mus_musculus='Mus_musculus',
            mmu='Mus_musculus')
  out <- opts[xx]
  if (is.na(out)) {
    stop("Illegal species value: ", x)
  }
  out
}
