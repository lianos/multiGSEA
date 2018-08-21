##' Get pathways/GOslim information from PANTHER.db Biocondcutor package.
##'
##' @description
##' This is a convience function that orchestrates the PANTHER.db package to
##' return GeneSetDb objects for either pathway or GOslim information for
##' human or mouse.
##'
##' Note that for some reason the \code{PANTHER.db} package needs to be
##' installed in a user-writable package location for this to work properly.
##' If you see an error like \code{Error in resqlite_send_query ... attempt to
##' write a readonly database}, this is the problem. Please install another
##' version of the \code{PANTHER.db} package in a user-writable directory using
##' \code{\link[BiocInstaller]{biocLite}}.
##'
##' @section GOSLIM:
##' \href{http://geneontology.org/page/go-slim-and-subset-guide}{GO Slims} are
##' "cut down" versions of the GO ontology that contain a subset of the terms in
##' the whole GO.
##'
##' PANTHER provides their own set of
##' \href{GO slims}{http://www.pantherdb.org/panther/ontologies.jsp}, although
##' it's not clear how often these get updated.
##'
##' @rdname getPantherGeneSetDb
##' @export
##' @importFrom GSEABase EntrezIdentifier
##' @param type "pathway" or, "goslim"
##' @param species "human" or "mouse"
##'
##' @return A wired up GeneSetDb
getPantherGeneSetDb <- function(type=c('pathway', 'goslim'),
                                species=c('human', 'mouse')) {
  species <- match.arg(species)
  type <- match.arg(type)

  if (!requireNamespace('PANTHER.db', quietly=TRUE)) {
    stop("The PANTHER.db bioconductor package is required")
  }
  if (species == 'human') {
    org.pkg <- 'org.Hs.eg.db'
    xorg <- 'Homo_sapiens'
  } else {
    org.pkg <- 'org.Mm.eg.db'
    xorg <- 'Mus_musculus'
  }
  on.exit({
    unloadNamespace('PANTHER.db')
    unloadNamespace(org.pkg)
  })
  if (!requireNamespace(org.pkg, quietly=TRUE)) {
    stop(org.pkg, " bioconductor package required for this species query")
  }

  p.db <- PANTHER.db::PANTHER.db
  PANTHER.db::pthOrganisms(p.db) <- toupper(species)
  org.db <- getFromNamespace(org.pkg, org.pkg)

  out <- switch(type,
                pathway=getPantherPathways(p.db, org.db),
                goslim=getPantherGOSLIM(p.db, org.db))
  mapIds <- getFromNamespace('mapIds', 'AnnotationDbi')
  out@db$symbol <- mapIds(org.db, out@db$featureId, 'SYMBOL', 'ENTREZID')
  org(out) <- xorg
  out
}

##' @rdname getPantherGeneSetDb
##' @export
getGOslimGeneSetDb <- function(species=c('human', 'mouse')) {
  getPantherGeneSetDb('goslim', species)
}

getPantherPathways <- function(p.db, org.db) {
  aselect <- getFromNamespace('select', 'AnnotationDbi')
  p.all <- aselect(p.db, AnnotationDbi::keys(p.db, keytype="PATHWAY_ID"),
                   columns=c("PATHWAY_ID", "PATHWAY_TERM", "UNIPROT"),
                   'PATHWAY_ID')
  ## Map uniprot to entrez
  umap <- aselect(org.db, p.all$UNIPROT, c('UNIPROT', 'ENTREZID'), 'UNIPROT')
  m <- merge(p.all, umap, by='UNIPROT')
  m <- m[!is.na(m[['ENTREZID']]),,drop=FALSE]
  lol <- list(`panther pathway`=split(m$ENTREZID, m$PATHWAY_TERM))

  idxref <- unique(as.data.table(p.all)[, c('PATHWAY_TERM', 'PATHWAY_ID'), with=FALSE])
  setkeyv(idxref, 'PATHWAY_TERM')

  url.fn <- function(coll, name) {
    ## Captures the lookup table for future use from parent.frame(!)
    ## Not sure if this will cause any serious memory leak, but ...
    pid <- idxref[list(name)]$PATHWAY_ID
    if (is.na(pid)) {
      return("http://pantherdb.org/panther/prowler.jsp?reset=1&selectedView=5")
    }
    paste0('http://pantherdb.org/pathway/pathwayDiagram.jsp?catAccession=', pid)
  }

  gdb <- GeneSetDb(lol)
  geneSetCollectionURLfunction(gdb, names(lol)) <- url.fn
  featureIdType(gdb, names(lol)) <- EntrezIdentifier()
  gdb
}

getPantherGOSLIM <- function(p.db, org.db) {
  if (!require("GO.db")) {
    stop("GO.db is required for this functionality")
  }
  aselect <- getFromNamespace('select', 'AnnotationDbi')
  p.all <- aselect(p.db,
                   AnnotationDbi::keys(p.db, keytype='GOSLIM_ID'),
                   columns=c('ENTREZ', 'GOSLIM_ID', 'GOSLIM_TERM'),
                   'GOSLIM_ID')
  p.all <- p.all[!is.na(p.all[['ENTREZ']]),,drop=FALSE]
  p.all <- p.all[order(p.all[['ENTREZ']]),]

  go <- aselect(GO.db::GO.db,
                unique(p.all[['GOSLIM_ID']]),
                c('GOID', 'TERM'),
                'GOID')
  go.missed <- go[is.na(go[['TERM']]),,drop=FALSE]
  ## 2015-09-02 (Bioc 3.1)
  ##       GOID TERM
  ## GO:0005083   NA
  ## GO:0006917   NA
  ## GO:0019204   NA
  go.add <- data.frame(
    GOID=c("GO:0005083", "GO:0006917", "GO:0019204"),
    TERM=c(
      "GTPase regulator activity",
      "apoptotic process",
      "nucleotide phosphatase activity (obsolete)"
    ), stringsAsFactors=FALSE)
  go <- rbind(
    go[!is.na(go[['TERM']]),,drop=FALSE],
    go.add
  )

  GO <- merge(p.all, go, by.x='GOSLIM_ID', by.y='GOID', all.x=TRUE)

  missed <- is.na(GO$TERM)
  mids <- unique(GO$GOSLIM_ID[missed])
  if (length(mids)) {
    warning(length(missed), " GOSLIM terms were not found in GO.db",
            immediate.=TRUE)
    GO$TERM <- ifelse(missed, GO$GOSLIM_ID, GO$TERM)
  }

  lol <- split(GO$ENTREZ, GO$TERM)
  url.fn <- function(x, y) 'http://www.pantherdb.org/panther/ontologies.jsp'

  gdb <- GeneSetDb(lol, collectionName='GOSLIM')
  xref <- match(gdb@table$name, GO$TERM)
  gdb@table$GOID <- GO$GOSLIM_ID[xref]
  gdb@table$ontology <- GO$GOSLIM_TERM[xref]
  geneSetCollectionURLfunction(gdb, 'GOSLIM') <- url.fn
  featureIdType(gdb, 'GOSLIM') <- EntrezIdentifier()
  gdb
}
