
#' Fetches a `GeneSetDb` from geneset collections defined in MSigDB.
#'
#' This provides versioned genesets from gene set collections defined in
#' [MSigDB](http://software.broadinstitute.org/gsea/msigdb). Collections can
#' be retrieved by their collection name, ie `c("H", "C2", "C7")`.
#'
#' Some subsets of curated genesets from within C2 can be retrieved by name,
#' like `"reactome"`, `"kegg"`, `"biocarta"`, and `"pid"`. You can, for
#' instance, call this function with`collection = c("reactome", "H")`, and
#' the reactome subset of C2 will be returned, along with all of the hallmark
#' genesets. When invoked like this, these "blessed" subsets of collections
#' will be promoted out of the C2 collection and into its own. This happens
#' when `promote_subcategory_to_collection = FALSE` (the default).
#'
#' The GO collection (C5) will also be promoted out of C5 and into their own
#' `"GO_MP"`, `"GO_BP"`, and `"GO_MF"` collections.
#'
#' @section KEGG Gene Sets:
#' Due to the licensing restrictions over the KEGG collections, they are not
#' returned from this function unless they are explicitly asked for. You can
#' ask for them through this function by either (i) querying for the `"c2"`
#' collection while setting `with.kegg = TRUE`; or (ii) explicitly calling with
#' `collection = "kegg"`.

#' @section MSigDB Versions:
#' We recently switched to using the msigdbr package as the source of truth for
#' these, so v7 is the earliest version of the MSigDB collections we make
#' available. Version 6 are available in the following (deprecated) packages:
#'
#' * https://github.com/lianos/GeneSetDb.MSigDB.Mmusculus.v61
#' * https://github.com/lianos/GeneSetDb.MSigDB.Hsapiens.v61
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
#' @importFrom msigdb.data msigdb_retrieve
#'
#' @param collection character vector specifying the collections you want
#'   (c1, c2, ..., c7, h). By default we load just the hallmark collecitons.
#'   Setting this to `NULL` loads all collections. Alternative you can also
#'   include named subsets of collections, like `"reactome"`. Refer to the
#'   Details section for more information.
#' @param species human or mouse?
#' @param with.kegg The Broad distributes the latest versions of the KEGG
#'   genesets as part of the c2 collection. These genesets come with a
#'   restricted license, so by default we do not return them as part of the
#'   GeneSetDb. To include the KEGG gene sets when asking for the c2
#'   collection, set this flag to `TRUE`.
#' @param allow_multimap,min_ortho_sources configure how to handle orthology
#'   mapping (allow multimappers, and what type of level of db suport required).
#'   See help in [msigdb.data::msigdb_retrieve()]
#' @param version the version of the MSigDB database to use.
#' @return a `GeneSetDb` object
#' @examples
#' \dontrun{
#'   gdb <- getMSigGeneSetDb(c("h", "reactome"), "human", "entrez")
#'   gdb.h.entrez <- getMSigGeneSetDb(c("h", "c2"), "human", "entrez")
#'   gdb.h.ens <- getMSigGeneSetDb(c("h", "c2"), "human", "ensembl")
#'   gdb.m.entrez <- getMSigGeneSetDb(c("h", "c2"), "mouse", "entrez")
#' }
getMSigGeneSetDb <- function(collection = "H",
                             species = "human",
                             id.type = c("ensembl", "entrez", "symbol"),
                             with.kegg = FALSE,
                             allow_multimap = TRUE, min_ortho_sources = 2,
                             promote_subcategory_to_collection = TRUE,
                             version = NULL, ...) {
  id.type <- match.arg(id.type)
  msig.db <- msigdb_retrieve(
    collection, species, id.type, allow_multimap = allow_multimap,
    min_ortho_sources = min_ortho_sources,
    promote_subcategory_to_collection = promote_subcategory_to_collection, ...)

  if (!with.kegg && !"kegg" %in% tolower(collection)) {
    # If we didn't ask for kegg explicitly, we will remove these collections
    # due to their licensing policy. Better to be conservative here than not.
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

  url.fn <- function(collection, name, ...) {
    url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
    sprintf(url, name)
  }

  # this is for when we promote certain subcategories to their own collections,
  # ie. promote_subcategory_to_collection = TRUE
  promoted.url.fn <- function(collection, name, ...) {
    name.prefix <- sub("_.*$", "", collection)
    new.name <- paste(name.prefix, name, sep = "_")
    url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
    sprintf(url, new.name)
  }

  for (col in unique(geneSets(gdb)$collection)) {
    is.promoted <- !toupper(col) %in% c("H", paste0("C", 1:7))
    fn <- if (is.promoted) promoted.url.fn else url.fn
    geneSetCollectionURLfunction(gdb, col) <- fn
    featureIdType(gdb, col) <- idtype
    gdb <- addCollectionMetadata(gdb, col, 'source',
                                 attr(msig.db, "msigdb_version"))
  }

  org(gdb) <- attr(msig.db, "species_info")[["species_name_"]]
  gdb@collectionMetadata <- gdb@collectionMetadata[name != "count"]
  gdb
}

