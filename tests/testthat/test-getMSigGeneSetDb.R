context("MSigDB GeneSetDb")

test_that("MSigDB retrieval respects collection subsets", {
  gdb.all <- getMSigGeneSetDb()
  expect_setequal(geneSets(gdb.all)$collection, c("H", paste0("C", 1:7)))
  gdb.sub <- getMSigGeneSetDb(c("H", "C6"))
  expect_setequal(geneSets(gdb.sub)$collection, c("H", "C6"))
})

test_that("with.kegg honors inclusion/exclusion of KEGG gene sets", {
  with.kegg <- getMSigGeneSetDb("c2", with.kegg = TRUE)
  no.kegg <- getMSigGeneSetDb("c2", with.kegg = FALSE)

  gs.kegg <- subset(geneSets(with.kegg), subcategory == "CP:KEGG")
  expect_true(nrow(gs.kegg) > 0L)

  gs.nokegg <- subset(geneSets(no.kegg), subcategory == "CP:KEGG")
  expect_true(nrow(gs.nokegg) == 0L)
})

test_that("url function stored correctly", {
  go.df <- GeneSetDb.MSigDB::msigdb_retrieve("human", "C5", "entrez")
  gdb.o <- GeneSetDb(go.df)

  # Let's change the collection to GO_MP, GO_BP, GO_CC and fix custom url
  # functions
  go.2 <- go.df %>%
    mutate(collection = paste0("GO_", subcategory),
           name = sub("^GO_", "", name)) %>%
    select(-subcategory)
  gdb.2 <- GeneSetDb(go.2)

  go.url.fn <- function(collection, name) {
    name. <- paste("GO", name, sep = "_")
    sprintf("http://www.broadinstitute.org/gsea/msigdb/cards/%s.html", name.)
  }

  for (col in unique(go.2$collection)) {
    geneSetCollectionURLfunction(gdb.2, col) <- go.url.fn
    featureIdType(gdb.2, col) <- GSEABase::EntrezIdentifier()
    gdb.2 <- addCollectionMetadata(gdb.2, col, "source", "v7")
  }

  expect_equal(
    basename(geneSetURL(gdb.2, "GO_MF", "ENZYME_ACTIVATOR_ACTIVITY")),
    "GO_ENZYME_ACTIVATOR_ACTIVITY.html", info = "MF")

  expect_equal(
    basename(geneSetURL(gdb.2, "GO_BP", "LIVER_REGENERATION")),
    "GO_LIVER_REGENERATION.html", info = "BP")

  expect_equal(
    basename(geneSetURL(gdb.2, "GO_CC", "GOLGI_APPARATUS")),
    "GO_GOLGI_APPARATUS.html", info = "CC")
})
