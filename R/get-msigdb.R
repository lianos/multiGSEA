.msigdb.species <- c('human', 'mouse')
.msigdb.collections <- list(
  # 'v5.1'=c(paste0('c', 1:7), 'h'),
  'v5.2'=c(paste0('c', 1:7), 'h'),
  'v6.1'=c(paste0('c', 1:7), 'h'))
.msigdb.version.current <- tail(names(.msigdb.collections), 1L)

#' Fetches a `GeneSetDb` from geneset collections defined in MSigDB.
#'
#' @description
#' This provides versioned genesets from gene set collections defined in
#' [MSigDB](http://software.broadinstitute.org/gsea/msigdb). We strive to keep
#' the latest version of MSigDB available via this function. The current version
#' provided is v6.1.
#'
#' Although the primariy identifiers used in MSigDB are human entrez IDs, we
#' also include gene symbols, and mouse versions of each gene set (except c1)
#' that were created from simple orthology mapping via biomaRt.
#'
#' We also include metadata at the geneset level by adding extra columns
#' to the `geneSets(gdb)` data.frame. "subcategory" may be most useful
#' for GO gene sets, eg. CC vs BP vs MF, and "organism" defines the organism in
#' which the geneset was derived from ("Homo sapiens" or "Mus musculus").
#'
#' Note that this function by default will exclude the KEGG genesets from the
#' MSigDB c2 collection. Set `with.kegg = TRUE` to include them. It is the
#' end user's responsibility to ensure that they are appropriately adhering
#' to the licensing restricions of the genesets provided herein.
#'
#' Note that this function is a wrapper that requires the GeneSetDb.MSigDB*
#' package (data files from AnnotationHub) to be installed/accessible.
#'
#' @export
#' @rdname MSigDB
#'
#' @param collection character vector specifying the collections you want
#'   (c1, c2, ..., c7, h)
#' @param species human or mouse?
#' @param with.kegg The Broad distributes the latest versions of the KEGG
#'   genesets as part of the c2 collection. These genesets come with a
#'   restricted license, so by default we do not return them as part of the
#'   GeneSetDb. To include the KEGG gene sets when asking for the c2
#'   collection, set this flag to `TRUE`.
#' @param species.specific Many of the genesets defined in MSigDB are annotated
#'   with the organism from which the experiment was conducted and which the
#'   geneset was extracted from. If this is set to `TRUE`, then only
#'   gene sets that match `species` will be included. Default is
#'   `FALSE`.
#' @param version the version of the MSigDB database to use.
#' @return a `GeneSetDb` object
#' \dontrun{
#'   gdb.h.entrez <- getMSigGeneSetDb(c("h", "c2"), "human", "entrez")
#'   gdb.h.ens <- getMSigGeneSetDb(c("h", "c2"), "human", "ensembl")
#'   gdb.m.entrez <- getMSigGeneSetDb(c("h", "c2"), "mouse", "entrez")
#' }
getMSigGeneSetDb <- function(collection, species='human',
                             id.type = c("entrez", "ensembl", "symbol"),
                             with.kegg=FALSE, species.specific=FALSE,
                             version=.msigdb.version.current) {
  species <- resolve.species(species)
  version <- match.arg(version, names(.msigdb.collections))
  id.type <- match.arg(id.type)
  avail.cols <- .msigdb.collections[[version]]

  pkg.species <- if (species == 'Mus_musculus') 'Mmusculus' else 'Hsapiens'
  pkg.version <- sub('\\.', '', version)

  ## Note: This will be updated when AnnotationHub files come oneline to
  ## support that properly, but right now I require these packages to be
  ## installed as the "more traditional" data packages.

  ## For the first releas, v5.2 is stored in the v6.1 package, so the general
  ## code below is special cased. This will change in the future after relase.
  # pkg <- sprintf("GeneSetDb.MSigDB.%s.%s", pkg.species, pkg.version)
  # pkg.dir <- find.package(pkg, quiet=TRUE)
  # if (length(pkg.dir) == 0L) {
  #   stop("GeneSetDb.* package not installed: ", pkg)
  # }
  # fn <- paste0(pkg, '.rds')
  # fn <- file.path(pkg.dir, 'extdata', fn)
  # out <- readRDS(fn)
  pkg <- sprintf("GeneSetDb.MSigDB.%s.v61", pkg.species)
  ## Check that package is installed so we can get the extdata
  pkg.dir <- find.package(pkg, quiet=TRUE)
  if (length(pkg.dir) == 0L) {
    stop("GeneSetDb.* package not installed: ", pkg)
  }

  use.symbols <- id.type == "symbol"
  if (use.symbols) {
    id.type <- "entrez" # this is the primary identifiers MSigDB provides
  }

  fn <- sub('\\.v61', paste0("-", id.type, ".", pkg.version, ".rds"), pkg)
  fn <- file.path(pkg.dir, 'extdata', fn)

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

  # Note: this functionality should really use the featureIdMap cross
  # referencing that you haven't made  yet.
  if (use.symbols) {
    out@db[[id.type]] <- out@db[["featureId"]]
    out@db[["featureId"]] <- out@db[["symbol"]]
    # Sometimes identifiers can't be mapped to symbols, ie. the symbol is NA
    # or different identifier map to the same symbol, so we handle that here.
    out@db <- out@db[!is.na(featureId)]
    out@db <- unique(out@db, by = c("collection", "name", "featureId"))
    # clean up ou@table be ensuring that N is correct
    new.N <- out@db[, list(N = .N), by = c("collection", "name")]
    out@table <- out@table[new.N]
    out@table[, N := i.N]
    out@table[, i.N := NULL]
    out@collectionMetadata <- out@collectionMetadata[collection %in% out@table$collection]
  }

  # NOTE: remove count collectionMetadata
  out@collectionMetadata <- out@collectionMetadata[name != "count"]

  out
}

#' @noRd
resolve.species <- function(x) {
  stopifnot(is.character(x) && length(character) == 1L)
  xx <- tolower(x)
  opts <- c(
    fly="Drosophila_melanogaster",
    human='Homo_sapiens',
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
n
