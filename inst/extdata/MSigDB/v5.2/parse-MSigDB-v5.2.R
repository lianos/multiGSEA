## This script is meant to be run in an R session with a working directory
## that is set to this directory.
##
## This script parses the msigdb_v5.0.xml file that has all of the geneset
## information into:
##   1. A GeneSetCollection (which is saved); and
##   2. A GeneSetDb version of (1): also saved.
##
## This file is downloaded from:
##   http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/5.2/msigdb_v5.2.xml
## Set the `fn` argument to the path that the xml file was saved to
##
## TODO: We need to do ortholog mapping to mouse IDs
library(GSEABase)
library(multiGSEA)
library(data.table)
library(XML)

fn <- 'msigdb_v5.2.xml'

## taken from GSEABase::getBroadSets.
## I am mofifying it here, because some 'factories' are not recognized
## (for category == 'archived' for instance) and I don't want the whole thing
## to fail
myXMLNodeToGeneSetFactory <- function(file, membersId = "MEMBERS_SYMBOLIZED") {
  .mkSplit <- function(x) {
    if (is.null(x))
      character(0)
    else unlist(strsplit(x, GSEABase:::.BROAD_SEPARATOR))
  }
  url <- NULL
  if (length(file) == 1) {
    isUri <- grep("^(http|ftp|file)://", file)
    if (length(isUri) == 1 && isUri == 1)
      url <- file
    else if (file.exists(file)) {
      full <- path.expand(file)
      if (length(grep("^/", full)) == 1)
        url <- paste("file:/", full, sep = "")
      else if (.Platform$OS.type == "windows")
        url <- paste("file:///", full, sep = "")
      else url = full
    }
  }
  symbolId <- SymbolIdentifier()
  function(node) {
    attrs <- as.list(xmlAttrs(node))
    args <- list(symbolId, setName = attrs[["STANDARD_NAME"]],
                 setIdentifier = attrs[["SYSTEMATIC_NAME"]], geneIds = unique(.mkSplit(attrs[[membersId]])),
                 organism = attrs[["ORGANISM"]], urls = c(getBroadSet = url,
                                                          .mkSplit(attrs[["EXTERNAL_DETAILS_URL"]])), collectionType = {
                                                            categories <- .mkSplit(attrs[["CATEGORY_CODE"]])
                                                            subcategories <- .mkSplit(attrs[["SUB_CATEGORY_CODE"]])
                                                            category <- subcategory <- as.character(NA)
                                                            if (length(categories) >= 1) category <- tolower(categories[[1]])
                                                            if (length(subcategories) >= 1) subcategory <- subcategories[[1]] else if (length(categories) >=
                                                                                                                                       2) subcategory <- categories[[2]]
                                                            if (length(categories) > 2 || (length(categories) >
                                                                                           1 && length(subcategories) != 0)) {
                                                              fmt <- "Broad 'CATEGORY_CODE' too long: '%s'"
                                                              txt <- paste(categories, collapse = "' '")
                                                              warning(sprintf(fmt, txt))
                                                            }
                                                            MyBroadCollection(category = mkScalar(category),
                                                                            subCategory = mkScalar(subcategory))
                                                          }, contributor = attrs[["CONTRIBUTOR"]], pubMedIds = attrs[["PMID"]],
                 shortDescription = attrs[["DESCRIPTION_BRIEF"]],
                 longDescription = attrs[["DESCRIPTION_FULL"]], TAGS = NULL,
                 MESH = NULL, CHIP = NULL, MEMBERS = NULL)
    args <- args[!sapply(args, is.null)]
    do.call(GeneSet, args)
  }
}

MyBroadCollection <- function(category = "c1", subCategory = NA, ...)  {
  # if (length(category) != 1 || !(category %in% c("c1", "c2",
  #                                                "c3", "c4", "c5", "c6", "c7", "h")))
  #   stop(sprintf("invalid BroadCollection category: '%s'",
  #                paste(category, collapse = "', '")))
  new("BroadCollection", category = mkScalar(category), subCategory = mkScalar(as.character(subCategory)))
}

.fromXML <- function(file, node, handler, ...) {
  res <- xmlTreeParse(file, useInternalNodes = TRUE, ...)
  geneSets <- getNodeSet(res, node, fun = handler)
  free(res)
  geneSets
}

getMsigSets <- function (uri, ..., membersId = c("MEMBERS_SYMBOLIZED", "MEMBERS_EZID"))
{
  membersId <- match.arg(membersId)
  factories <- sapply(uri, myXMLNodeToGeneSetFactory, membersId = membersId)

  # tryCatch({
  #   geneSets <- unlist(mapply(GSEABase:::.fromXML, uri, "//GENESET",
  #                             factories, SIMPLIFY = FALSE, USE.NAMES = FALSE))
  # }, error = function(err) {
  #   msg <- paste("'getMsigSets' failed to create gene sets:\n  ",
  #                conditionMessage(err))
  #   msg
  # })
  # # GeneSetCollection(geneSets)
  # geneSets
  unlist(mapply(GSEABase:::.fromXML, uri, "//GENESET",
                factories, SIMPLIFY = FALSE, USE.NAMES = FALSE))

}


gs.list <- getMsigSets(fn, membersId="MEMBERS_EZID")
saveRDS(gs.list ,'GeneSet-list-v52-MEMBERS_EZID.rds')

gsc <- GeneSetCollection(gs.list)

if (FALSE) {
  g.go <- gsc[[1]]
  g.h <- gsc[['HALLMARK_ANGIOGENESIS']]
  sapply(slotNames(g.h), function(x) slot(g.h, x), simplify=FALSE)
}

info <- lapply(gsc, function(gs) {
  data.table(collection=bcCategory(collectionType(gs)),
             name=setName(gs),
             organism=organism(gs),
             subcategory=gs@collectionType@subCategory,
             featureId=geneIds(gs))
})
info.all <- rbindlist(info)
for (col in names(info.all)) {
  info.all[, (col) := as.character(info.all[[col]])]
}

info <- subset(info.all, collection %in% c('h', paste0('c', 1:7)))
setdiff(info.all$collection, info$collection) ## "archived"

library(org.Hs.eg.db)
info[, symbol := mapIds(org.Hs.eg.db, featureId, 'SYMBOL', 'ENTREZID')]

gdb <- GeneSetDb(info[, list(collection, name, featureId, symbol)])

## Update the gdb@table
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
  gdb <- addCollectionMetadata(gdb, col, 'source', 'MSigDB_v5.2')
}

org(gdb) <- 'Homo_sapiens'

gdb.fn <- sprintf('MSigDB.Homo_sapiens.GeneSetDb.rds', species)
saveRDS(gdb, gdb.fn)

## Convert GeneSetDb to mouse --------------------------------------------------
## I use IGIS here, so let's do it on rescomp
## TODO: This isn't done yet
# library(rmd.plugins)
# library(GSEABase)
# library(multiGSEA)
library(igis)

# stopifnot(rmd.plugins:::on.rescomp())
#gdbh <- readRDS('MSigDB.Homo_sapiens.GeneSetDb.rds')

gdbh <- gdb
rm(gdb)
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
setnames(mm.db, c('featureId', 'symbol', 'featureId.Mm'), c('featureId.Hs', 'symbol.Hs', 'featureId'))

library(org.Mm.eg.db)
mm.db$symbol <- mapIds(org.Mm.eg.db, mm.db$featureId, 'SYMBOL', 'ENTREZID')
setcolorder(mm.db, c(names(gdbh@db), setdiff(names(mm.db), names(gdbh@db))))

## we don't want to make a geneset for c1, since position info doesn't translate
## to mouse
mm.db.all <- mm.db
mm.db <- subset(mm.db.all, collection != 'c1')


gdb <- GeneSetDb(mm.db)

## Spot check that the featureId symbol featureId.Hs symbol.Hs columns in
## gdb@db match up, then remove the *.Hs columns (justs add more HD space)
gdb@db[sample(nrow(gdb@db), 20)]
gdb@db <- gdb@db[, 1:4, with=FALSE]

url.fn <- function(collection, name) {
  url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  sprintf(url, name)
}
for (col in unique(geneSets(gdb)$collection)) {
  geneSetCollectionURLfunction(gdb, col) <- url.fn
  featureIdType(gdb, col) <- EntrezIdentifier()
  gdb <- addCollectionMetadata(gdb, col, 'source', 'MSigDB_v5.2')
}

org(gdb) <- 'Mus_musculus'

## Update the gdb@table
gdb@table <- local({
  otable <- gdbh@table[, c('collection', 'name', setdiff(names(gdbh@table), names(gdb@table))), with=FALSE]
  gst <- merge(gdb@table, otable, by=c('collection', 'name'))
  stopifnot(sum(duplicated(gst$name)) == 0)
  setkeyv(gst, key(gdb@table))
  stopifnot(all.equal(gst[, names(gdb@table), with=FALSE], gdb@table))
  gst
})

gdb.fn <- sprintf('MSigDB.Mus_musculus.GeneSetDb.rds', species)
saveRDS(gdb, gdb.fn)
