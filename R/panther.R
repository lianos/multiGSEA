##' Get pathways/GOslim information from PANTHER.db Biocondcutor package.
##'
##' @export
##'
##' @param type "pathway" or, "goslim"
##' @param species "human" or "mouse"
##'
##' @return A wired up GeneSetDb
getPanther <- function(type=c('pathway', 'goslim'),
                       species=c('human', 'mouse')) {
  if (!require('PANTHER.db')) {
    stop("The PANTHER.db bioconductor package is required")
  }
  species <- match.arg(species)
  type <- match.arg(type)

  if (species == 'human') {
    org.pkg <- 'org.Hs.eg.db'
    xorg <- 'Homo_sapiens'
  } else {
    org.pkg <- 'org.Mm.eg.db'
    xorg <- 'Mus_musculus'
  }
  if (!require(org.pkg, character.only=TRUE)) {
    stop(org.pkg, " bioconductor package required for this species query")
  }

  species(PANTHER.db) <- toupper(species)
  org.db <- get(org.pkg)

  out <- switch(type,
                pathway=getPantherPathways(PANTHER.db, org.db),
                goslim=getPantherGOSLIM(PANTHER.db, org.db))
  org(out) <- xorg
  out
}

getPantherPathways <- function(p.db, org.db) {
  p.all <- select(p.db, keys(p.db, keytype="PATHWAY_ID"),
                  columns=c("PATHWAY_ID", "PATHWAY_TERM", "UNIPROT"),
                  'PATHWAY_ID')
  ## Map uniprot to entrez
  umap <- select(org.db, p.all$UNIPROT, c('UNIPROT', 'ENTREZID'), 'UNIPROT')
  m <- merge(p.all, umap, by='UNIPROT')
  m <- subset(m, !is.na(ENTREZID))
  lol <- list(`panther pathway`=split(m$ENTREZID, m$PATHWAY_TERM))

  idxref <- unique(as.data.table(p.all)[, list(PATHWAY_TERM, PATHWAY_ID)])
  setkeyv(idxref, 'PATHWAY_TERM')

  url.fn <- function(coll, name) {
    ## Captures the lookup table for future use from parent.frame(!)
    ## Not sure if this will cause any serious memory leak, but ...
    pid <- idxref[J(name)]$PATHWAY_ID
    if (is.na(pid)) {
      return("http://pantherdb.org/panther/prowler.jsp?reset=1&selectedView=5")
    }
    paste0('http://pantherdb.org/pathway/pathwayDiagram.jsp?catAccession=', pid)
  }

  gdb <- GeneSetDb(lol)
  collectionUrlFunction(gdb, names(lol)) <- url.fn
  featureIdType(gdb, names(lol)) <- EntrezIdentifier()
  gdb
}

getPantherGOSLIM <- function(p.db, org.db) {
  if (!require("GO.db")) {
    stop("GO.db is required for this functionality")
  }
  p.all <- select(p.db,
                  keys(p.db, keytype='GOSLIM_ID'),
                  columns=c('ENTREZ', 'GOSLIM_ID', 'GOSLIM_TERM'),
                  'GOSLIM_ID')
  p.all <- p.all[order(p.all$ENTREZ),]

  library(GO.db)
  go <- select(GO.db,
               unique(p.all$GOSLIM_ID),
               c('GOID', 'TERM'),
               'GOID')
  go.missed <- subset(go, is.na(TERM))
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
  go <- rbind(subset(go, !is.na(TERM)), go.add)

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
  collectionUrlFunction(gdb, 'GOSLIM') <- url.fn
  featureIdType(gdb, 'GOSLIM') <- EntrezIdentifier()
  gdb
}
