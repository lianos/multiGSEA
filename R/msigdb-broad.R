.msigdb.species <- c('human', 'mouse')
.msigdb.collections <- list(
  'v4.0'=paste0('c', 1:7),
  'v5.0'=c(paste0('c', 1:7), 'h'))
.msigdb.version.current <- tail(names(.msigdb.collections), 1L)

##' Fetches an MSigDB GeneSetDb
##'
##' The gene sets for version 5.0 were parsed straight from MSigDB's XML dumps.
##' Version 4.0 was based off of the data provided by WEHI:
##'   http://bioinf.wehi.edu.au/software/MSigDB/
##'
##' @export
##'
##' @param collection character vector specifying the collections you want
##'   (c1, c2, ..., c7, h)
##' @param species human or mouse?
##' @param version the version of the MSigDB database to use.
##' @return a \code{GeneSetDb} object
getMSigDBset <- function(collection, species='human',
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

  if (species == 'Mus_musculus') {
    cmd <- paste0('getMSigDBset(c(',
                  paste(sprintf('"%s"', collection), collapse=","),
                  '), species="mouse", version="v4.0")')
    stop("MSigDB v5.0 GeneSetDb has not yet been built for v5.0.\n",
         "  The v4.0 signatures can still be retrieved like so:\n\n    ",
         cmd, "\n")
  }

  fn <- sprintf('MSigDB.%s.GeneSetDb.rds', species)
  fn <- system.file('extdata', 'MSigDB', version, fn, package='multiGSEA')
  out <- readRDS(fn)

  if (!setequal(collection, avail.cols)) {
    keep <- geneSets(out)$collection %in% collection
    out <- subset.GeneSetDb(out, geneSets(out)$collection %in% collection)
  }

  out
}

resolve.species <- function(x) {
  stopifnot(is.character(x) && length(character) == 1L)
  xx <- tolower(x)
  opts <- c(human='Homo_sapiens', homo_sapiens='Homo_sapiens',
            mouse='Mus_musculus', mus_musculus='Mus_musculus')
  out <- opts[xx]
  if (is.na(out)) {
    stop("Illegal species value: ", x)
  }
  out
}
