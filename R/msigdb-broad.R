.msigdb.species <- c('human', 'mouse')
.msigdb.collections <- list(
  'v5.1'=c(paste0('c', 1:7), 'h'),
  'v5.2'=c(paste0('c', 1:7), 'h'),
  'v6.1'=c(paste0('c', 1:7), 'h'))
.msigdb.version.current <- tail(names(.msigdb.collections), 1L)

##' Fetches a \code{GeneSetDb} from geneset collections defined in MSigDB.
##'
##' @description
##' This provides versioned genesets from gene set collections defined in
##' \href{http://software.broadinstitute.org/gsea/msigdb}{MSigDB}. We strive
##' to keep the latest version of MSigDB included in this package. The current
##' version provided is v5.2.
##'
##' Although the primariy identifiers used in MSigDB are human entrez IDs, we
##' also include gene symbols, and mouse versions of each gene set (except c1)
##' that were created from simple orthology mapping via biomaRt.
##'
##' We also include metadata at the geneset level by adding extra columns
##' to the \code{geneSets(gdb)} data.frame. "subcategory" may be most useful
##' for GO gene sets, eg. CC vs BP vs MF, and "organism" defines the organism in
##' which the geneset was derived from ("Homo sapiens" or "Mus musculus").
##'
##' Note that this function by default will exclude the KEGG genesets from the
##' MSigDB c2 collection. Set \code{with.kegg=TRUE} to include them. It is the
##' end user's responsibility to ensure that they are appropriately adhering
##' to the licensing restricions of the genesets provided herein.
##'
##' @export
##' @rdname MSigDB
##' @param collection character vector specifying the collections you want
##'   (c1, c2, ..., c7, h)
##' @param species human or mouse?
##' @param with.kegg The Broad distributes the latest versions of the KEGG
##'   genesets as part of the c2 collection. These genesets come with a
##'   restricted license, so by default we do not return them as part of the
##'   GeneSetDb. To include the KEGG gene sets when asking for the c2
##'   collection, set this flag to \code{TRUE}.
##' @param species.specific Many of the genesets defined in MSigDB are annotated
##'   with the organism from which the experiment was conducted and which the
##'   geneset was extracted from. If this is set to \code{TRUE}, then only
##'   gene sets that match \code{species} will be included. Default is
##'   \code{FALSE}.
##' @param version the version of the MSigDB database to use.
##' @return a \code{GeneSetDb} object
getMSigGeneSetDb <- function(collection, species='human', with.kegg=FALSE,
                             species.specific=FALSE,
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

  fn <- sprintf('MSigDB.%s.GeneSetDb.rds', species)
  fn <- system.file('extdata', 'MSigDB', version, fn, package='multiGSEA')
  out <- readRDS(fn)

  gs <- geneSets(out, as.dt=TRUE)

  if (!setequal(collection, avail.cols)) {
    keep <- gs$collection %in% collection
    out <- subset.GeneSetDb(out, gs$collection %in% collection)
    gs <- geneSets(out, as.dt=TRUE)
  }

  if ('c2' %in% gs$collection && !with.kegg) {
    is.c2 <- gs$collection == 'c2'
    keep <- !is.c2 | (is.c2 & !grepl('^KEGG_', gs$name))
    out <- subset.GeneSetDb(out, keep)
  }

  if (species.specific) {
    spec <- sub('_', ' ', species)
    org <- geneSets(out)$organism
    keep <-  !is.na(org) & org == spec
    out <- subset.GeneSetDb(out, keep)
  }

  out
}

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
