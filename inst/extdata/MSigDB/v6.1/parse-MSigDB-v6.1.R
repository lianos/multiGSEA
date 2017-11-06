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

fn <- '~/Downloads/msigdb_v6.1.xml'

## taken from GSEABase::getBroadSets.
## I am modifying it here, because some 'factories' are not recognized
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
saveRDS(gs.list ,'~/Downloads/GeneSet-list-v62-MEMBERS_EZID.rds')

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
  gdb <- addCollectionMetadata(gdb, col, 'source', 'MSigDB_v6.1')
}

org(gdb) <- 'Homo_sapiens'

# gdb.fn <- sprintf('MSigDB.Homo_sapiens.GeneSetDb.rds', species)
gdb.fn <- 'MSigDB.Homo_sapiens.GeneSetDb.rds'
saveRDS(gdb, gdb.fn)

## Convert GeneSetDb to mouse --------------------------------------------------
## We use the biomaRt::getLDS (get linked datasets) function to map human
## entrez ids to mouse. Look at ?getLDS function in biomaRt
library(biomaRt)
library(dplyr)

hdf <- gdb %>%
  as.data.frame %>%
  group_by(collection, name) %>%
  mutate(n=n()) %>%
  ungroup
hids <- hdf %>%
  select(entrezgene=featureId, symbol) %>%
  distinct(entrezgene, .keep_all=TRUE) %>%
  mutate(featureId=entrezgene)

human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
xref <- getLDS(attributes = c("hgnc_symbol", "entrezgene", "ensembl_gene_id"),
               filters = "entrezgene", values = hids$entrezgene, mart = human,
               attributesL = c("ensembl_gene_id", "entrezgene", "mgi_symbol"),
               martL = mouse)
xfer <- xref %>%
  select(featureId=NCBI.gene.ID, mm_featureId=NCBI.gene.ID.1, mm_symbol=MGI.symbol) %>%
  distinct(featureId, mm_featureId, .keep_all = TRUE) %>%
  mutate(featureId=as.character(featureId), mm_featureId=as.character(mm_featureId))

mdf <- hdf %>%
  filter(collection != 'c1') %>%
  inner_join(xfer, by="featureId") %>%
  group_by(collection, name) %>%
  mutate(nm=n()) %>%
  ungroup

mdf.go <- mdf %>%
  select(collection, name, featureId=mm_featureId, symbol=mm_symbol) %>%
  filter(!is.na(featureId), !is.na(symbol),
         nchar(featureId) > 0, nchar(symbol) > 0)

mgdb <- GeneSetDb(mdf.go)
## add metadata from gdb@table
take.cols <- c('collection', 'name', setdiff(colnames(gdb@table), colnames(mgdb@table)))
meta <- gdb@table[mgdb@table, take.cols, with=FALSE]
mnew <- mgdb@table[meta]
stopifnot(
  all.equal(mgdb@table[, list(collection, name)], mnew[, list(collection, name)]),
  all.equal(mgdb@table$N, mnew$N))
mgdb@table <- mnew

url.fn <- function(collection, name) {
  url <- "http://www.broadinstitute.org/gsea/msigdb/cards/%s.html"
  sprintf(url, name)
}
for (col in unique(geneSets(mgdb)$collection)) {
  geneSetCollectionURLfunction(mgdb, col) <- url.fn
  featureIdType(mgdb, col) <- EntrezIdentifier()
  mgdb <- addCollectionMetadata(mgdb, col, 'source', 'MSigDB_v6.1')
}

org(mgdb) <- 'Mus_musculus'
# gdb.fn <- sprintf('MSigDB.Mus_musculus.GeneSetDb.rds', species)
gdb.fn <- 'MSigDB.Mus_musculus.GeneSetDb.rds'
saveRDS(mgdb, gdb.fn)
