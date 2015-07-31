##' Get pathways/GOslim information from PANTHER.db Biocondcutor package.
##'
##' Currently we only support the pathway collections from panther
##'
##' @export
##'
##' @param type "pathway" (or, eventually, "GOslim")
##' @param species "human" or "mouse"
##'
##' @return A wired up GeneSetDb
getPanther <- function(type=c('pathway', 'GOslim'),
                       species=c('human', 'mouse')) {
  if (!require('PANTHER.db')) {
    stop("The PANTHER.db bioconductor package is required")
  }
  species <- match.arg(species)
  type <- match.arg(type)
  if (type != 'pathway') {
    stop("only pathway supported now")
  }

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

  ## p.db <- PANTHER.db
  species(PANTHER.db) <- toupper(species)
  org.db <- get(org.pkg)

  out <- getPantherPathways(PANTHER.db, org.db)
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

## Scratch ---------------------------------------------------------------------
if (FALSE) {
  ## Just exploring what's in here
  pathways <- select(PANTHER.db, keys(PANTHER.db, keytype="PATHWAY_ID"),
                     columns=c("PATHWAY_ID", "PATHWAY_TERM",
                               "CLASS_ID", "CLASS_TERM",
                               "COMPONENT_ID", "COMPONENT_TERM"),
                     'PATHWAY_ID')
  p <- as.data.table(pathways)
}
