context("Retrieiving MSigDB GeneSetDb objects")

test_that("Specifying id.type = 'symbol' keeps GeneSetDb honest", {
  species <- c("human", "mouse")

  for (s in species) {
    gdbo <- getMSigGeneSetDb(c("h", "c7"), species = s, id.type = "entrez")

    gdbs <- getMSigGeneSetDb(c("h", "c7"), species = s, id.type = "symbol")
    expect_true(validObject(gdbs))

    gso <- geneSets(gdbo, as.dt = TRUE)
    gss <- geneSets(gdbs, as.dt = TRUE)

    # everything is the same in these databases except the geneset sizes
    expect_equal(
      gss[, setdiff(colnames(gss), "N"), with = FALSE],
      gso[, setdiff(colnames(gso), "N"), with = FALSE],
      info = s)

    # Compare gene counts per gene set. The gene sets in the symbol db
    # can only be the same size as the original, or smaller.
    gst <- gso[gss[, list(collection, name, N)]]
    expect_true(all(gst$N >= gst$i.N), info = s)
    expect_true(any(gst$N >  gst$i.N), info = s)

    # featureId's in the new @db are symbols from original one
    omap <- unique(gdbo@db[, list(featureId, symbol)])
    setkeyv(omap, "featureId")
    smap <- unique(gdbs@db[, list(featureId = entrez, symbol = featureId)])
    setkeyv(smap, "featureId")
    xmap <- smap[omap]
    # $symbol has values from the "symbol mapped" gdb
    # $i.symbol has the symbol value for the feature in the "enetrez mapped" gdb

    # where is.na($symbol): there was no symbol provided.
    no.symbol <- xmap[is.na(symbol)]

    # features with no symbol should not appear in smap$featureId
    expect_true(sum(no.symbol$featureId %in% smap$featureId) == 0)

    # symbols in smap should be same as the original GeneSetDb
    has.symbol <- xmap[!is.na(symbol)]
    # Note that in mouse, some ids are mapped to, ie.
    # entrez_id 67118 is mapped to both 3110001I22Rik and Bfar.
    # Let's remove these for this test.
    multi.id <- has.symbol$featureId[duplicated(has.symbol$featureId)]
    has.symbol <- xmap[!is.na(symbol) & !featureId %in% multi.id]
    multi.map <- subset(has.symbol, featureId %in% multi.id)

    expect_equal(has.symbol$symbol, has.symbol$i.symbol)
  }
})
