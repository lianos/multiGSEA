context("GeneSetDb")

## TODO: These tests need to be added to test-GeneSetDb
##  * test the URL functions
##  * test various collectionMetadata manipulation, ie:
##      - collectionMetadata(x)
##      - collectionMetadata(x, collection)
##      - collectionMetadata(x, collection, name)
##  * test collectionMetadata<- ensures single collection,name pairs

gdb.h <- getMSigGeneSetDb(c('h'))
gdb.c6 <- getMSigGeneSetDb(c('c6'))

test_that("GeneSetDb constructor preserves featureIDs per geneset", {
  ## This test exercise both the single list and list-of-lists input for geneset
  ## membership info.
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  expect_is(gsd, 'GeneSetDb')

  ## ids in gsd@db should match input lists
  ## This test is not exercising the API that fetches featureIds, but rather
  ## retrieves them using the back door -- this is intentional.
  for (xgrp in names(gsl)) {
    for (xid in names(gsl[[xgrp]])) {
      gsd.ids <- gsd@db[J(xgrp, xid)]$featureId
      expected.ids <- gsl[[xgrp]][[xid]]
      info.lol <- sprintf('collection %s id %s', xgrp, xid)
      expect_true(setequal(expected.ids, gsd.ids), info=info.lol)
      expect_false(any(duplicated(gsd.ids)), info=info.lol)
    }
  }
})

test_that("GeneSetDb constructor works with an input data.frame", {
  gdb0 <- GeneSetDb(exampleGeneSets())
  df <- as.data.frame(gdb0)

  ## Adding a fake symbol here, just to see what is the what
  meta <- data.frame(featureId=unique(df$featureId), stringsAsFactors=FALSE)
  faux <- replicate(nrow(meta),
                    paste(sample(letters, 5, replace=TRUE), collapse=""))
  meta$symbol <- faux
  df.in <- merge(df, meta, by='featureId')

  ## A warning is fired if merging extra columns (symbol, here) hoses something
  ## in the GeneSetDb, so let's make sure there is no such warning here.
  gdb <- GeneSetDb(df.in[sample(nrow(df.in)),]) ## randomize rows for fun
  expect_equal(gdb, gdb0, features.only=TRUE)

  ## Check that the symbol column from df.in was added to gdb@db
  expect_is(gdb@db$symbol, 'character')
})


test_that("GeneSetDb contructor converts GeneSetCollection properly", {
  gsc <- as(gdb.h, 'GeneSetCollection')
  gdbn <- GeneSetDb(gsc, collectionName='h')
  expect_equal(gdb.h, gdbn, features.only=TRUE)
})

test_that("GeneSetDb contructor converts list of GeneSetCollection properly", {
  gdbo <- append(gdb.h, gdb.c6)

  gscl <- list(h=as(gdb.h, 'GeneSetCollection'),
               c6=as(gdb.c6, 'GeneSetCollection'))
  gdbn <- GeneSetDb(gscl)

  ## Ensure that collection names are preserved, since gscl is a named list
  ## of collections
  expect_equal(gdbn, gdbo, features.only=TRUE)
})

test_that("GeneSetDb constructor honors custom collectionName args", {
  gdbo <- append(gdb.h, gdb.c6)

  lol <- as.list(gdbo, nested=TRUE)
  new.cnames <- setNames(c('x1', 'x2'), names(lol))

  ## Change collectionName from h,c6 to c2,c1
  gdbn <- GeneSetDb(lol, collectionName=new.cnames)
  gso <- geneSets(gdbo)
  for (oname in names(new.cnames)) {
    nname <- new.cnames[oname]
    gs.names <- subset(geneSets(gdbo), collection == oname)$name
    for (gs.name in gs.names) {
      oids <- featureIds(gdbo, oname, gs.name)
      nids <- featureIds(gdbn, nname, gs.name)
      expect_true(setequal(nids, oids),
                  info=sprintf("featureId parity for (%s:%s, %s)",
                               oname, nname, gs.name))
    }
  }
})

test_that("as(gdb, 'GeneSetCollection') preserves featureIds per GeneSet", {
  gdb <- getMSigGeneSetDb(c('h', 'c6'))
  gsc <- as(gdb, 'GeneSetCollection')
  for (gs in gsc) {
    gs.info <- strsplit(GSEABase::setName(gs), ';')[[1]]
    coll <- gs.info[1]
    name <- gs.info[2]
    expect_true(setequal(GSEABase::geneIds(gs), featureIds(gdb, coll, name)),
                info=sprintf("featureId match for geneset (%s,%s)", coll, name))
  }
})

test_that("featureIds(GeneSetDb, i, j) accessor works", {
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)

  for (group in names(gsl)) {
    for (id in names(gsl[[group]])) {
      expected.ids <- gsl[[group]][[id]]
      gsd.ids <- featureIds(gsd, group, id)
      msg <- sprintf("unexpected ids returned featureIds(gsd, %s, %s)", group, id)
      expect_true(setequal(expected.ids, gsd.ids), info=msg)
    }
  }
})

## This test and the conform,GeneSetDb test below are testing similar things
test_that("featureIds(GeneSetDb, i, j) removes 'unconformable' featureIds", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsc <- conform(gsd, vm)

  ## This assumes all genesets passed in are "active"
  for (group in names(gsl)) {
    for (id in names(gsl[[group]])) {
      all.ids <- gsl[[group]][[id]]
      expected.ids <- intersect(all.ids, rownames(vm))
      gsc.ids.all <- featureIds(gsc, group, id, fetch.all=TRUE)
      gsc.ids <- featureIds(gsc, group, id)
      msg <- "unexpected ids returned featureIds(gsd, %s, %s), fetch.all=%s"
      expect_true(setequal(expected.ids, gsc.ids),
                  info=sprintf(group, id, FALSE))
      expect_true(setequal(all.ids, gsc.ids.all),
                  info=sprintf(group, id, TRUE))
    }
  }
})

test_that("featureIds(GeneSetDb, i, MISSING) gets all features in a collection", {
  vm <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsc <- conform(gsd, vm)

  cols <- unique(geneSets(gsd)$collection)
  for (col in cols) {
    fids <- featureIds(gsd, col)
    expected <- unique(subset(gsd@db, collection == col)$featureId)
    expect_true(setequal(expected, fids), info=paste("collection:", col))
  }

  ## Uncformable features dropped
  for (col in cols) {
    fids <- featureIds(gsc, col)
    expected <- intersect(subset(gsd@db, collection == col)$featureId,
                          rownames(vm))
    expect_true(setequal(expected, fids),
                info=paste("active collection:", col))
  }
})

test_that("conform,GeneSetDb follows row permutation in expression object", {
  es <- exampleExpressionSet(do.voom=FALSE)
  es.mixed <- es[sample(nrow(es))]
  es.sub <- es[sample(nrow(es), 100)]

  gsd <- GeneSetDb(exampleGeneSets())

  gsd.es <- conform(gsd, es)
  gsd.mixed <- conform(gsd, es.mixed)

  ## gsd.sub is super small, so we expect a warning due to not being able to
  ## match many featureIds to the expression object
  expect_warning({
    gsd.sub <- conform(gsd, es.sub)
  }, "^fraction .* low:", ignore.case=TRUE)

  ## gsd.es and gsd.mixed should have the same features in them but different
  ## x.id
  expect_equal(geneSets(gsd.es), geneSets(gsd.mixed))
  gt <- geneSets(gsd.es)
  gt.sub <- geneSets(gsd.sub)

  for (i in 1:nrow(gt)) {
    grp <- gt$collection[i]
    xid <- gt$name[i]
    n <- gt$n[i]
    label <- sprintf("(%s, %s)", gt$collection[i], gt$name[i])
    ## The feature IDs of fids.es and fids.mix must be the same
    fids.es <- featureIds(gsd.es, grp, xid)
    fids.mix <- featureIds(gsd.mixed, grp, xid)
    expect_equal(n, length(fids.es))
    expect_true(setequal(fids.es, fids.mix))

    ## Ensure that the $x.idx's for each match the rownames of the expression
    ## object
    es.rn <- rownames(es)[featureIds(gsd.es, grp, xid, 'x.idx')]
    es.mixed.rn <- rownames(es.mixed)[featureIds(gsd.mixed, grp, xid, 'x.idx')]
    expect_true(setequal(es.rn, es.mixed.rn), info=label)
    expect_true(setequal(es.rn, fids.es), info=label)

    ## Was this geneset deactivated in the subset gt?
    if (is.active(gsd.sub, grp, xid)) {
      fids.sub <- featureIds(gsd.sub, grp, xid)
      n.sub <- gsd.sub@table[J(grp, xid)]$n
      expect_equal(n.sub, length(fids.sub))
      expect_true(length(fids.sub) <= length(fids.es))
      ## sub features are subset of original features
      expect_true(length(setdiff(fids.sub, fids.es)) == 0)

      ## check that the rownames of the expression object match the features
      ## returned "by index"
      es.sub.rn <- rownames(es.sub)[featureIds(gsd.sub, grp, xid, 'x.idx')]
      expect_true(setequal(fids.sub, es.sub.rn),
                  info=sprintf("gsd.sub: %s,%s", grp, xid))
    }
  }
})

test_that("append,GeneSetDb works", {
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  extra <- list(more=list(first=head(letters, 10), second=tail(letters, 10)))
  gst2 <- append(gsd, extra)
  expect_is(gst2, 'GeneSetDb')
  expect_true(validObject(gst2))

  all.gsl <- c(gsl, extra)

  ## Ensure that all new and old features are in the new GeneSetDb
  for (group in names(all.gsl)) {
    for (id in names(all.gsl[[group]])) {
      info <- sprintf('appended collection: %s, name %s', group, id)
      expected <- all.gsl[[group]][[id]]
      fids <- featureIds(gst2, group, id)
      expect_true(setequal(expected, fids), info=info)
    }
  }
})

test_that("append,GeneSetDb honors geneset metadata in columns of geneSets()", {
  m <- getMSigGeneSetDb('h')
  r <- getReactomeGeneSetDb()

  a <- append(r, m)

  ## Check that all columns are there
  expect_true(setequal(c(names(geneSets(m)), names(geneSets(r))),
                       names(geneSets(a))))

  ## Add a species column to m@table to check if it carries through after append
  m@table$species <- 'human'
  a2 <- append(m, r)
  expect_true(setequal(c(names(geneSets(m)), names(geneSets(r))),
                       names(geneSets(a2))))

})

test_that("as.list.GeneSetDb returns gene sets in same order as GeneSetDb", {
  es <- exampleExpressionSet()
  gsd <- conform(exampleGeneSetDb(), es)

  gs.idxs <- as.list(gsd)
  info <- strsplit(names(gs.idxs), ';;')
  res <- data.table(collection=sapply(info,'[[',1L), name=sapply(info,'[[',2L))
  expected <- geneSets(gsd, active.only=TRUE, .external=FALSE)
  expect_equal(res, expected[, list(collection, name)], check.attributes=FALSE)
})

test_that("as.list.GeneSetDb returns proper indexes into conformed object", {
  es <- exampleExpressionSet()
  gsi <- exampleGeneSets(es)

  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), es)
  indexes <- as.list(gsd, 'x.idx', nested=TRUE)

  for (xgrp in names(gsl)) {
    for (xid in names(gsl[[xgrp]])) {
      expected <- match(gsl[[xgrp]][[xid]], rownames(es))
      expected <- expected[!is.na(expected)]
      gsd.idxs <- indexes[[xgrp]][[xid]]
      expect_true(setequal(expected, gsd.idxs),
                  info=sprintf("%s,%s", xgrp, xid))
    }
  }
})

test_that("conformed & unconformed GeneSetDb,incidenceMatrix is kosher", {
  es <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- GeneSetDb(gsl)
  gsdc <- conform(gsd, es)

  im <- incidenceMatrix(gsd)
  imc <- incidenceMatrix(gsdc)
  g.cols <- sub(';;.*', '', rownames(im))
  g.names <- sub('.*?;;', '', rownames(im))
  for (i in 1:nrow(im)) {
    col <- g.cols[i]
    name <- g.names[i]
    fidx <- im[i,] == 1
    fidxc <- imc[i,] == 1

    fids <- gsl[[col]][[name]]
    fidsc <- intersect(fids, rownames(es))

    im.fids <- colnames(im)[fidx]
    imc.fids <- colnames(imc)[fidxc]
    ## Check ids from incidence matrix
    expect_true(setequal(im.fids, fids), info=paste(i, "unconformed GeneSetDb"))
    expect_true(setequal(imc.fids, fidsc), info=paste(i ,"conformed GeneSetDb"))
  }
})

test_that("annotateGeneSetMembership works", {
  vm <- exampleExpressionSet(do.voom=TRUE)
  gdb <- GeneSetDb(exampleGeneSets())
  mg <- multiGSEA(gdb, vm, vm$design, ncol(vm$design), NULL)
  lfc <- logFC(mg)

  ## Test that annotation is consistent with pre-conformed gdb vs uncormed
  lfc.anno.u <- annotateGeneSetMembership(lfc, gdb)    ## unconformed
  lfc.anno.p <- annotateGeneSetMembership(lfc, mg@gsd) ## preconformed
  expect_equal(lfc.anno.u, lfc.anno.p)

  ## ensure that annotateGeneSetMembership guessed the right column
  lfc.anno.x <- annotateGeneSetMembership(lfc, gdb, x.ids=lfc$featureId)
  expect_equal(lfc.anno.x, lfc.anno.u)

  ## ensure that x.ids specified by column name in lfc works
  lfc.anno.c <- annotateGeneSetMembership(lfc, gdb, x.ids='featureId')
  expect_equal(lfc.anno.x, lfc.anno.c)
})

test_that("subsetByFeatures returns correct genesets for features", {
  set.seed(0xBEEEF)
  gdb <- exampleGeneSetDb()
  features <- sample(featureIds(gdb), 10)
  gdb.sub <- subsetByFeatures(gdb, features)

  db.all <- gdb@db
  db.sub <- gdb.sub@db
  db.rest <- anti_join(db.all, db.sub, by=c('collection', 'name'))

  ## db.sub + db.rest should == db.all
  expect_equal(nrow(db.sub) + nrow(db.rest), nrow(db.all))

  ## 1. Ensure that each geneset in subsetted gdb (gdb.sub) has >= 1
  ##    of requested features in its featureId column.
  has.1 <- db.sub[, {
    list(N=.N, n=sum(featureId %in% features))
  }, by=c('collection', 'name')]
  expect_true(all(has.1$n >= 1))

  ## 2. Ensure that any geneset not in the subsetted GeneSetDb doesn't have
  ##    any of the requested features
  has.0 <- db.rest[, {
    list(N=.N, n=sum(featureId %in% features))
  }, by=c('collection', 'featureId')]
  expect_true(all(has.0$n) == 0)
})

test_that('subset.GeneSetDb ("[".GeneSetDb) creates valid result', {
  # es <- exampleExpressionSet()
  set.seed(1234)
  gdb <- exampleGeneSetDb()
  keep <- sample(c(TRUE, FALSE), length(gdb), replace=TRUE)
  sdb <- gdb[keep]
  expect_equal(length(sdb), sum(keep))
  expect_equal(geneSets(sdb)$name, geneSets(gdb)$name[keep])

  ## check counts
  sdb.coll.counts <- collectionMetadata(sdb, .external=FALSE)[name == 'count']
  sdb.coll.counts[, value := unlist(value)]
  exp.coll.counts <- geneSets(sdb, .external=FALSE)[, {
    .(name='count', value=.N)
  }, by='collection']
  setkeyv(exp.coll.counts, c('collection', 'name'))
  expect_equal(sdb.coll.counts, exp.coll.counts)
  expect_true(validObject(sdb))
})

test_that("GeneSetDb indexing `[` works", {
  ## TODO: Test indexing GeneSetDbs
})
