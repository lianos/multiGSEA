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
#' @param id.type Gene identifiers are returned by the REST service as
#'   entrez identifiers. Set this to `"ensembl"` to translate them internally
#'   using [remap_identifiers()]. If `species`is not `"human"` or `"mouse"`,
#'   you need to provide an idxref table that works with [remap_identifiers()].
#' @return A GeneSetDb of the kegg stuffs
#' @examples
#' \dontrun{
#' mouse.entrez <- getKeggGeneSetDb("mouse")
#' mouse.ens <- getKeggGeneSetDb("mouse", id.type = "ensembl")
#' }
getKeggGeneSetDb <- function(species = "human",
                             id.type = c("ensembl", "entrez"),
                             database = "pathway",
                             idxref = NULL, ...) {
  sinfo <- species_info(species)
  database <- match.arg(database)
  id.type <- match.arg(id.type)
  fn <- getFunction(sprintf(".get_kegg_%s_db", database))
  gdb <- .get_kegg_pathway_db(sinfo, id.type, idxref)
  gdb
}

.get_kegg_pathway_db <- function(sinfo, id.type = c("ensembl", "entrez"),
                                 idxref = NULL, ...) {
  if (FALSE) {
    species.code <- "hsa" # human
    species.code <- "mmu" # mouse
    species.code <- "mcf" # cyno
  }
  species.code <- sinfo[["kegg"]]
  id.type <- match.arg(id.type)
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

  if (id.type == "ensembl") {
    # I am doing the args -> do.call mojo so that we can use the original_id
    # and target_id parameters if supplied by the user when they supply their
    # own idxref table, otherwise we use the defaults
    args <- list(...)
    if (is.null(idxref)) {
      args[["xref"]] <- sinfo[["alias"]]
      args[["original_id"]] <- "entrezgene_id"
      args[["target_id"]] <- "ensembl_gene_id"
    } else {
      args[["xref"]] <- idxref
    }
    args[["x"]] <- out
    out <- do.call(remap_identifiers, args)
  } else {
    featureIdType(out, "KEGG") <- GSEABase::EntrezIdentifier()
  }

  out
}
