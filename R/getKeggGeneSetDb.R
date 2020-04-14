#' Retrieves a GeneSetDb from KEGG's REST API
#'
#' Uses [limma::getGeneKEGGLinks()] and [limma::getKEGGPathwayNames()]
#' internally.
#'
#' Currently we just support the pathway database, and only entrez ids.
#'
#' Note that **it is your responsibility** to ensure that you can use the KEGG
#' database according to their licensing requirements.
#'
#' @export
#' @param species `"human"`, `"mouse"` or any of the bioconductor or kegg-style
#'   abbreviations.
#' @param id.type entrez, enembl?
#' @return A GeneSetDb of the kegg stuffs
getKeggGeneSetDb <- function(species = "human", id.type = "entrez",
                             database = "pathway", ...) {
  sinfo <- species_info(species)
  database <- match.arg(database)
  fn <- getFunction(sprintf(".get_kegg_%s_db", database))
  gdb <- .get_kegg_pathway_db(sinfo$kegg)
  gdb
}

.get_kegg_pathway_db <- function(species.code, id.type = NULL, ...) {
  if (FALSE) {
    species.code <- "hsa" # human
    species.code <- "mmu" # mouse
    species.code <- "mcf" # cyno
  }
  pnames <- limma::getKEGGPathwayNames(species.code, remove.qualifier = TRUE)
  colnames(pnames) <- c("pathway_id", "name")
  membership <- limma::getGeneKEGGLinks(species.code)
  colnames(membership) <- c("feature_id", "pathway_id")

  df. <- merge(membership, pnames, by = "pathway_id")
  df.[["pathway_id"]] <- sub("path:", "", df.[["pathway_id"]])
  df.[["collection"]] <- "KEGG"
  out <- GeneSetDb(df.)

  fn <- function(collection, name, gdb, ...) {
    if (missing(gdb) || !is(gdb, "GeneSetDb")) {
      return("https://www.kegg.jp/kegg/pathway.html")
    }
    gs <- geneSets(gdb)
    pathway_id <- gs[gs[["name"]] == name,][["pathway_id"]]
    sprintf("https://www.genome.jp/dbget-bin/www_bget?%s", pathway_id)
  }
  geneSetCollectionURLfunction(out, "KEGG") <- fn
  featureIdType(out, "KEGG") <- GSEABase::EntrezIdentifier()
  out
}
