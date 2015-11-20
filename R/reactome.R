##' @importMethodsFrom AnnotationDbi select
##' @importMethodsFrom AnnotationDbi keys
getReactomeGeneSetDb <- function(species='human', rm.species.prefix=TRUE) {
  if (!requireNamespace('reactome.db')) {
    stop("reactome.db package required for this functionality")
  }
  species <- resolve.species(species)

  ## Find all KEGG pathways for the given species.
  ## Pathways are prefixed with the organism name like so:
  ##  `<GENUS> <SPECIES>: <PATHWAY NAME>`
  org.prefix <- paste0('^', sub('_', ' ', species), ': *')
  pathnames <- keys(reactome.db::reactome.db, keytype='PATHNAME')
  org.keep <- grepl(org.prefix, pathnames)
  org.pathnames <- pathnames[org.keep]

  info <- suppressWarnings({
    ## Generates 1:many mapping because of ENTREZID
    select(reactome.db::reactome.db,
           columns=c('PATHID', 'PATHNAME', 'ENTREZID'),
           keys=org.pathnames,
           keytype='PATHNAME')
  })
  stopifnot(setequal(org.pathnames, info$PATHNAME))
  info <- subset(info, !is.na(ENTREZID))
  setDT(info)

  ## Are there multiple mappings for PATHNAME:PATHID combo?
  ## maybe from different organisms?
  u.id2name <- unique(info[, list(PATHID, PATHNAME)])
  u.id2name[, N.pathid := .N, by='PATHID']
  if (nrow(dups <- subset(u.id2name, N.pathid > 1))) {
    warning("Multiple PATHID to PATHNAME not resolved", immediate.=TRUE)
  }

  if (rm.species.prefix) {
    info[, PATHNAME := sub(org.prefix, '', PATHNAME)]
  }

  ## Somehow only one of the 'name's of the reatome genesets is encoded in
  ## UTF-8:
  ##   Loss of proteins required for interphase microtubule organization
  ##   from the centrosome
  ## All of the rest are "unknown". This causes annoying warnings when
  ## data.table tries to join on collection,name -- so I'm just nuking the
  ## encoding here.
  info <- info[, list(collection='reactome', name=PATHNAME, featureId=ENTREZID,
                      pathId=PATHID)]
  Encoding(info$name) <- 'unknown'
  gdb <- GeneSetDb(info)

  ## fids <- split(info$ENTREZID, info$PATHNAME)
  ## lol <- list(reactome=fids)
  ## gdb <- GeneSetDb(lol)
  ## Add reactome PATHID to geneSets()
  ## gdb@table[, PATHID := info$PATHID[match(name, info$name)]]

  gdb
}
