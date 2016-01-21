## This script is meant to be run in an R session with a working directory
## that is set to this directory.
##
## This script parses the msigdb_v5.0.xml file that has all of the geneset
## information into:
##   1. A GeneSetCollection (which is saved); and
##   2. A GeneSetDb version of (1): also saved.
##
## This file is downloaded from:
##   http://www.broadinstitute.org/gsea/msigdb/download_file.jsp?\
##     filePath=/resources/msigdb/5.0/msigdb_v5.0.xml
## Set the `fn` argument to the path that the xml file was saved to
##
## TODO: We need to do ortholog mapping to mouse IDs
library(GSEABase)
library(multiGSEA)
library(data.table)

fn <- '~/Downloads/msigdb_v5.0.xml'
fn <- '~/Downloads/msigdb_v5.1.xml'
gsc <- getBroadSets(fn, membersId="MEMBERS_EZID")
saveRDS(gsc ,'GeneSetCollection-v51-MEMBERS_EZID.rds')

g.go <- gsc[[1]]
g.h <- gsc[['HALLMARK_ANGIOGENESIS']]
sapply(slotNames(g.h), function(x) slot(g.h, x), simplify=FALSE)

info <- lapply(gsc, function(gs) {
  data.table(collection=bcCategory(collectionType(gs)),
             name=setName(gs),
             organism=organism(gs),
             subcategory=gs@collectionType@subCategory,
             featureId=geneIds(gs))
})
info <- rbindlist(info)
for (col in names(info)) {
  info[, (col) := as.character(info[[col]])]
}

library(org.Hs.eg.db)
info[, symbol := mapIds(org.Hs.eg.db, featureId, 'SYMBOL', 'ENTREZID')]

gdb <- GeneSetDb(info[, list(collection, name, featureId, symbol)])
gdb@table <- local({
  gs <- unique(info[, list(collection, name, subcategory, organism)])
  stopifnot(sum(duplicated(gs$name)) == 0)
  gst <- merge(gdb@table, gs, by=c('collection', 'name'))
  setkeyv(gst, key(gdb@table))
  stopifnot(all.equal(gst[, names(gdb@table), with=FALSE], gdb@table))
  gst
})

## Beef up collectionMetadata --------------------------------------------------
## URL function
url.fn <- function(collection, name) {
  url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  sprintf(url, name)
}
for (col in unique(geneSets(gdb)$collection)) {
  geneSetCollectionURLfunction(gdb, col) <- url.fn
  featureIdType(gdb, col) <- EntrezIdentifier()
  gdb <- addCollectionMetadata(gdb, col, 'source', 'MSigDB_v5.1')
}

org(gdb) <- 'Homo_sapiens'

gdb.fn <- sprintf('MSigDB.Homo_sapiens.GeneSetDb.rds', species)
saveRDS(gdb, gdb.fn)

## Convert GeneSetDb to mouse --------------------------------------------------
## I use IGIS here, so let's do it on rescomp
## TODO: This isn't done yet
library(rmd.plugins)
library(GSEABase)
library(multiGSEA)
library(igis)

stopifnot(rmd.plugins:::on.rescomp())
gdbh <- readRDS('MSigDB.Homo_sapiens.GeneSetDb.rds')
hs.ids <- unique(gdbh@db$featureId)
ortho <- local({
  ## the ortholog mapping returns a 1:many result from Hs to Mm, but this is OK
  ## for now
  o <- orthologs(geneIds=paste0('GeneID:', hs.ids), toSpecies='mouse')
  dt <- data.table(featureId=sub('GeneID:', '', o$Entrez),
                   featureId.Mm=sub('GeneID:', '', o$mouse_ensembl_gene))
  unique(dt)[!is.na(featureId.Mm) & nchar(featureId.Mm) > 0]
})

## I want to keep around the sourge human entrez ID that mapped to mouse
mm.db <- merge(gdbh@db, ortho, by='featureId')
setnames(mm.db, c('featureId', 'featureId.Mm'), c('featureId.Hs', 'featureId'))
setcolorder(mm.db, c(names(gdbh@db), 'featureId.Hs'))

## we don't want to make a geneset for c1, since position info doesn't translate
## to mouse
lol <- sapply(setdiff(mm.db$collection, 'c1'), function(xcol) {
  with(subset(mm.db, collection == xcol), {
    split(featureId, name)
  })
}, simplify=FALSE)

gdb <- GeneSetDb(lol)

url.fn <- function(collection, name) {
  url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  sprintf(url, name)
}
for (col in names(lol)) {
  multiGSEA:::geneSetCollectionURLfunction(gdb, col) <- url.fn
}

## Add id_type for the identifiers and species
gdb@collectionMetadata <- local({
  cm <- gdb@collectionMetadata
  more <- lapply(names(lol), function(col) {
    data.table(
      collection=col,
      name=c('organism', 'id_type', 'source'),
      value=list('Mus_musculus', EntrezIdentifier(), 'MSigDB_v5.0'))
  })
  more <- rbindlist(more)
  out <- rbind(cm, more)
  setkeyv(out, key(cm))
})

## Add the human entrezId each feature was mapped from in gdb@db
xmap <- data.table(featureId=ortho$featureId.Mm, featureId.Hs=ortho$featureId)
xmap <- unique(xmap, by='featureId')
udb <- merge(gdb@db, xmap, by='featureId')
setcolorder(udb, c(names(gdb@db), 'featureId.Hs'))
udb <- udb[order(collection, name, featureId)]
setkeyv(udb, key(gdb@db))

## Check that these are same as original db
odb <- gdb@db[order(collection, name, featureId)]
stopifnot(all.equal(udb$collection, odb$collection))
stopifnot(all.equal(udb$name, odb$name))
stopifnot(all.equal(udb$featureId, odb$featureId))

gdb@db <- udb
stopifnot(validObject(gdb))

gdb.fn <- sprintf('MSigDB.Mus_musculus.GeneSetDb.rds', species)
saveRDS(gdb, gdb.fn)
