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
#' @section Citing the Molecular Signatures Database:
#' To cite your use of the Molecular Signatures Database (MSigDB), please
#' reference Subramanian, Tamayo, et al. (2005, PNAS 102, 15545-15550) and one
#' or more of the following as appropriate:
#'
#' * Liberzon, et al. (2011, Bionformatics);
#' * Liberzon, et al. (2015, Cell Systems); and
#' * The source for the gene set as listed on the gene set page.
#'
#' @export
#' @importFrom GeneSetDb.MSigDB msigdb_retrieve
#'
#' @param collection character vector specifying the collections you want
#'   (c1, c2, ..., c7, h)
#' @param species human or mouse?
#' @param with.kegg The Broad distributes the latest versions of the KEGG
#'   genesets as part of the c2 collection. These genesets come with a
#'   restricted license, so by default we do not return them as part of the
#'   GeneSetDb. To include the KEGG gene sets when asking for the c2
#'   collection, set this flag to `TRUE`.
#' @param allow_multimap,min_ortho_sources configure how to handle orthology
#'   mapping (allow multimappers, and what type of level of db suport required).
#'   See help in [GeneSetDb.MSigDB::msigdb_retrieve()]
#' @param version the version of the MSigDB database to use.
#' @return a `GeneSetDb` object
#' @examples
#' \dontrun{
#'   gdb.h.entrez <- getMSigGeneSetDb(c("h", "c2"), "human", "entrez")
#'   gdb.h.ens <- getMSigGeneSetDb(c("h", "c2"), "human", "ensembl")
#'   gdb.m.entrez <- getMSigGeneSetDb(c("h", "c2"), "mouse", "entrez")
#' }
getMSigGeneSetDb <- function(collection = c("H", paste0("C", 1:7)),
                             species = "human",
                             id.type = c("ensembl", "entrez", "symbol"),
                             with.kegg = FALSE,
                             allow_multimap = TRUE, min_ortho_sources = 2,
                             version = NULL, ...) {
  id.type <- match.arg(id.type)
  msig.db <- msigdb_retrieve(species, collection, id.type, ...)
  if (!with.kegg) {
    msig.db <- subset(msig.db, !subcategory %in% "CP:KEGG")
  }
  gdb <- GeneSetDb(msig.db)

  # Beef up collectionMetadata
  if (id.type == "ensembl") {
    idtype <- GSEABase::ENSEMBLIdentifier()
  } else if (id.type == "entrez") {
    idtype <- GSEABase::EntrezIdentifier()
  } else {
    idtype <- GSEABase::SymbolIdentifier()
  }

  url.fn <- function(collection, name) {
    url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
    sprintf(url, name)
  }
  for (col in unique(geneSets(gdb)$collection)) {
    geneSetCollectionURLfunction(gdb, col) <- url.fn
    featureIdType(gdb, col) <- idtype
    gdb <- addCollectionMetadata(gdb, col, 'source',
                                 attr(msig.db, "msigdb_version"))
  }

  org(gdb) <- attr(msig.db, "species_info")[["species_name_"]]
  gdb@collectionMetadata <- gdb@collectionMetadata[name != "count"]
  gdb
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

