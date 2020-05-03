#' Converts the feature identifiers in a GeneSetDb to a set of new ones.
#'
#' The various GeneSetDb data providers limit the identifier types that they
#' provide. Use this function to map the given identifiers to whichever type
#' you like.
#'
#' You need to provie a data.frame that has a column for the current identifiers
#' and another column for the target identifiers.
#'
#' When there are multiple target id's for the source id, they will all be
#' returned. When there is no target id for the source id, the soure feature
#' will be axed.
#'
#' By default, the first column of the `xref` data.frame is assumed to be
#' the source identifiers, and the second is the target_id you want to remap
#' the identifiers to.
#'
#' @export
#' @param x The GeneSetDb with identifiers to convert
#' @param xref If a species (`"mouse"` or `"human"`) we will load the internal
#'   ensembl <-> entrez data.frame, or you can provide your own.
#' @return a remapped GeneSetDb object
#' @examples
#' gdb.entrez <- exampleGeneSetDb()
#' gdb.ens <- remap_identifiers(gdb.entrez, "human",
#'                              original_id = "entrezgene_id",
#'                              target_id = "ensembl_gene_id")
remap_identifiers <- function(x, xref, original_id = colnames(xref)[1L],
                              target_id = colnames(xref)[2L], ...) {
  assert_class(x, "GeneSetDb")
  if (test_string(xref)) xref <- load_id_xref(xref)
  assert_multi_class(xref, c("data.frame", "data.table", "tbl"))
  assert_string(original_id)
  assert_string(target_id)
  assert_subset(c(original_id, target_id), colnames(xref))

  db <- merge(x@db, xref, by.x = "feature_id", by.y = original_id,
              suffixes = c(".original", ""))
  setnames(db, "feature_id", "original_id")
  setnames(db, target_id, "feature_id")

  db <- db[!is.na(feature_id) & nchar(feature_id) > 0]
  db <- unique(db, by = c("collection", "name", "feature_id"))
  gs.dt <- merge(db, geneSets(x, as.dt = TRUE), by = c("collection", "name"))
  gs.dt[, N := NULL]
  gs.dt[, n := NULL]
  gs.dt[, active := NULL]

  out <- GeneSetDb(gs.dt)
  out@collectionMetadata <- x@collectionMetadata
  out
}

#' Load the mouse or human ensembl <-> entrez id maps.
#'
#' For convenience, entrez <-> ensembl maps come bundled with the is package.
#' While there are other sactioned ways to achieve this within the bioconductor
#' ecosystem, we bundle it here to make our lives easier.
#'
#' These tables are gnerated by the `inst/scripts/genereate-id-maps.R`
#'
#' @export
#' @param species The name of the species to load the entrez <-> ensembl xref
#'   table for. Currently we only provide tables for human and mouse ... and we
#'   shouldn't even provide those in here, but ...
#' @return a table of ensembl and entrez ids for each species
load_id_xref <- function(species, ..., as.dt = FALSE) {
  assert_string(species)
  sinfo <- species_info(species)
  fn <- sprintf("%s-entrez-ensembl.csv.gz", sinfo[["alias"]])
  fn <- system.file("extdata", "identifiers", fn, package = "multiGSEA")
  if (!test_file_exists(fn)) {
    stop("Can not find `extdata/idenifiers` xref table for: ", species)
  }
  out <- data.table::fread(fn)
  out[, entrezgene_id := as.character(entrezgene_id)]
  if (!as.dt) out <- setDF(copy(out))
  out
}
