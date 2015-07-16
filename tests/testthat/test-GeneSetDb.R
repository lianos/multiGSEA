context("GeneSetDb")

## TODO: These tests need to be added to test-GeneSetDb
##  * test the URL functions
##  * test various collectionMetadata manipulation, ie:
##      - collectionMetadata(x)
##      - collectionMetadata(x, collection)
##      - collectionMetadata(x, collection, name)
##  * test collectionMetadata<- ensures single collection,name pairs

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

test_that("GeneSetDb contructor converts GeneSetCollection properly", {
  gdbo <- getMSigDBGeneSetDb('h')
  gsc <- as(gdbo, 'GeneSetCollection')
  gdbn <- GeneSetDb(gsc, collectionName='h')
  expect_equal(gdbo, gdbn, features.only=TRUE)
})

test_that("GeneSetDb contructor converts list of GeneSetCollection properly", {
  gdb.h <- getMSigDBGeneSetDb(c('h'))
  gdb.c6 <- getMSigDBGeneSetDb(c('c6'))
  gdbo <- append(gdb.h, gdb.c6)

  gscl <- list(h=as(gdb.h, 'GeneSetCollection'),
               c6=as(gdb.c6, 'GeneSetCollection'))
  gdbn <- GeneSetDb(gscl)

  ## Ensure that collection names are preserved, since gscl is a named list
  ## of collections
  expect_true(setequal(geneSets(gdbn)$collection, geneSets(gdbo)$collection),
              info='Collection names preserved from named list of collections')
  expect_equal(gdbn, gdbo, features.only=TRUE,
               info="feature parity between GeneSetDb")
})

test_that("GeneSetDb constructor honors custom collectionName args", {
  gdb.h <- getMSigDBGeneSetDb(c('h'))
  gdb.c6 <- getMSigDBGeneSetDb(c('c6'))
  gdbo <- append(gdb.h, gdb.c6)

  lol <- as.list(gdbo, nested=TRUE)
  new.cnames <- setNames(c('c1', 'c2'), names(lol))

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

test_that("as(gdb, 'GeneSetDb') preserves featureIds per GeneSet", {
  gdb <- getMSigDBGeneSetDb(c('h', 'c6'))
  gsc <- as(gdb, 'GeneSetCollection')
  for (gs in gsc) {
    gs.info <- strsplit(setName(gs), ';')[[1]]
    coll <- gs.info[1]
    name <- gs.info[2]
    expect_true(setequal(geneIds(gs), featureIds(gdb, coll, name)),
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
  expect_warning(gsd.sub <- conform(gsd, es.sub),
                 "^fraction .* low:", ignore.case=TRUE)

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
  m <- getMSigDBset('h')
  r <- getReactomeGeneSetDb()
  a <- append(r, m)

  ## Check that all columns are there
  expect_true(setequal(c(names(geneSets(m)), names(geneSets(r))),
                       names(geneSets(a))))

  m@table$species <- 'human'
  a2 <- append(m, r)
  expect_true(setequal(c(names(geneSets(m)), names(geneSets(r))),
                       names(geneSets(a2))))

})

test_that("GeneSetDb returns proper limma lists via as.expression.indexes", {
  es <- exampleExpressionSet()
  gsi <- exampleGeneSets(es)

  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), es)
  indexes <- multiGSEA:::as.expression.indexes(gsd, 'x.idx')

  for (xgrp in names(gsl)) {
    for (xid in names(gsl[[xgrp]])) {
      expected <- match(gsl[[xgrp]][[xid]], rownames(es))
      expected <- expected[!is.na(expected)]
      gsd.idxs <- featureIds(gsd, xgrp, xid, value='x.idx')
      expect_true(setequal(expected, gsd.idxs),
                  info=sprintf("%s,%s", xgrp, xid))
    }
  }
})

test_that("GeneSetDb,incidenceMatrix is correct", {
  es <- exampleExpressionSet()
  gsl <- exampleGeneSets()
  gsd <- conform(GeneSetDb(gsl), es)

  im <- incidenceMatrix(gsd, es)
  g.cols <- sub(';.*', '', rownames(im))
  g.names <- sub('.*?;', '', rownames(im))
  for (i in 1:nrow(im)) {
    col <- g.cols[i]
    name <- g.names[i]
    fids <- featureIds(gsd, col, name)
    expected <- intersect(gsl[[col]][[name]], rownames(es))
    expect_true(setequal(expected, fids))
  }
})

test_that("subset.GeneSetDb works", {
  ## TODO: Test subset.GeneSetDb
})

test_that("GeneSetDb indexing `[` works", {
  ## TODO: Test indexing GeneSetDbs
})
