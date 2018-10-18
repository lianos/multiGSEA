#' Retrieves a housekeeping GeneSetDb for a given species
#'
#' The genesets provided here were extracted from the genes used in the
#' Qiagen RT2 Profiler PCR Arrays.
#'
#' Sources for the mouse and human genesets can be found here:
#'
#' * mouse: https://www.qiagen.com/us/shop/pcr/primer-sets/rt2-profiler-pcr-arrays/?catno=PAMM-000Z#geneglobe
#' * human: https://www.qiagen.com/us/shop/pcr/primer-sets/rt2-profiler-pcr-arrays/?catno=PAHS-000Z#geneglobe
#'
#' @export
#' @importFrom GSEABase EntrezIdentifier ENSEMBLIdentifier
#' @md
#'
#' @param species `"human"` or `"mouse"`
#' @param id.type defines the type of identifers used in the `featureId` column
#'   for the GeneSetDb, `"entrez"` or `"ensembl"`.
#' @return a GeneSetDb of housekeeping genes
getHousekeepingGeneSetDb <- function(species=c("human", "mouse"),
                                     id.type = c("entrez", "ensembl"), ...) {
  species <- match.arg(species)
  id.type <- match.arg(id.type)

  fn <- system.file("extdata", "signatures", "housekeeping",
                    "qiagen-rt2-pcr-array.csv", package = "multiGSEA",
                    mustWork = TRUE)

  gsets <- read.csv(fn, colClasses = "character")
  gsets[["featureId"]] <- gsets[[id.type]]
  gsets[["collection"]] <- "qiagen"
  gsets[["name"]] <- "housekeeping"

  gsets <- subset(gsets, organism == species)
  gdb <- GeneSetDb(gsets[, c("collection", "name", "featureId", "symbol")])

  geneSetCollectionURLfunction(gdb, 'qiagen') <- qiagen.url.fn(species)

  idtype <- switch(id.type,
                   entrez = EntrezIdentifier(),
                   ensembl = ENSEMBLIdentifier())
  featureIdType(gdb, 'qiagen') <- idtype

  organism <- switch(species, human = "Homo_sapiens", mouse = "Mus_musculus")
  org(gdb, "qiagen") <- organism

  addCollectionMetadata(gdb, "qiagen", "source",
                        "Qiagen RT2 Profiler (Housekeeping) PCR Arrays")
}

#' @noRd
qiagen.url.fn <- function(species = c("human", "mouse")) {
  species <- match.arg(species)

  if (species == "human") {
    catno <- c(housekeeping = "PAHS-000Z")
  } else {
    catno <- c(housekeeping = "PAMM-000Z")
  }

  url.fn <- function(collection, name) {
    if (name %in% names(catno)) {
      url <- paste0("https://www.qiagen.com/us/shop/pcr/primer-sets",
                    "/rt2-profiler-pcr-arrays/?catno=", catno[name],
                    "#geneglobe")
    } else {
      url <- paste0("https://www.qiagen.com/us/shop/pcr/",
                    "real-time-pcr-enzymes-and-kits/two-step-qrt-pcr/",
                    "rt2-profiler-pcr-arrays/#orderinginformation")
    }
    url
  }

  url.fn
}

