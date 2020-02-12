context("PANTHER.db")

# test_that("PANTHER.db GOSLIM categories work", {
#   ## PANTHER.db v1.0.2 had messed up internal database mappings, so just want
#   ## to ensure that fetching ontologies works.
#   mdb <- getGOslimGeneSetDb('mouse')
#   hdb <- getGOslimGeneSetDb('human')
#   expect_false(any(mdb@db$feature_id %in% hdb@db$feature_id))
#   expect_true(setequal(geneSets(mdb)$GOID, geneSets(hdb)$GOID))
# })
