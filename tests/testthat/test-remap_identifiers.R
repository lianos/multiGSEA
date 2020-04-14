context("remap_identifiers")

test_that("entrez id's remapped to ensemble", {
  xref.fn <- system.file("extdata", "identifiers",
                        "human-entrez-ensembl.csv.gz",
                        package = "multiGSEA")
  xref <- data.table::fread(xref.fn, colClasses = "character")
  gdb.entrez <- exampleGeneSetDb()
  gdb.ens <- remap_identifiers(gdb.entrez, xref,
                              original_id = "entrezgene_id",
                              target_id = "ensembl_gene_id")

  gs.entrez <- geneSet(gdb.entrez, name = "REACTOME_RAF_MAP_KINASE_CASCADE")
  gs.ens <- geneSet(gdb.ens, name = "REACTOME_RAF_MAP_KINASE_CASCADE")
  expect_true(all(substr(gs.ens$feature_id, 1, 4) == "ENSG"))
  expect_set_equal(gs.ens$original_id, gs.entrez$feature_id)
})
