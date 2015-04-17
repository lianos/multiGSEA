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

fn <- '~/Downloads/msigdb_v5.0.xml'
gsc <- getBroadSets(fn, membersId="MEMBERS_EZID")
saveRDS(gsc ,'GeneSetCollection-v5-MEMBES_EZID.rds')

g.go <- gsc[[1]]
g.h <- gsc[['HALLMARK_ANGIOGENESIS']]
sapply(slotNames(g.h), function(x) slot(g.h, x), simplify=FALSE)

infos <- lapply(gsc, function(gs) {
  data.table(collection=as.character(bcCategory(collectionType(gs))),
             name=as.character(setName(gs)),
             featureId=as.character(geneIds(gs)))
})
info <- rbindlist(infos)
table(unique(info, by='name')$collection)

lol <- sapply(unique(info$collection), function(col) {
  with(subset(info, collection == col), {
    split(featureId, name)
  })
}, simplify=FALSE)

gdb <- GeneSetDb(lol)

## Beef up collectionMetadata --------------------------------------------------
## URL function
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
      value=list('Homo_sapiens', EntrezIdentifier(), 'MSigDB_v5.0'))
  })
  more <- rbindlist(more)
  out <- rbind(cm, more)
  setkeyv(out, key(cm))
})

gdb.fn <- sprintf('MSigDB.Homo_sapiens.GeneSetDb.rds', species)
saveRDS(gdb, gdb.fn)

## Convert GeneSetDb to mouse --------------------------------------------------
