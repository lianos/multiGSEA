#' Retrieve annotations from reactome.db
#'
#' @importFrom GSEABase EntrezIdentifier
#' @export
#'
#' @param species the species to get pathay information for
#' @param rm.species.prefix pathways are provided with species prefixes from
#'   `reactome.db`, when `TRUE` (default), these are stripped from the gene set
#'   names.
#' @param id.type `"entrez"` or `"ensembl"`
#' @return a reactome GeneSetDb object
getReactomeGeneSetDb <- function(species = 'human',
                                 id.type = c("entrez", "ensembl"),
                                 rm.species.prefix = TRUE) {
  id.type <- match.arg(id.type)
  id.col <- if (id.type == "entrez") "ENTREZID" else "ENSG"
  rdb <- tryCatch(loadNamespace('reactome.db'), error=function(e) NULL)
  if (is.null(rdb)) {
    stop("reactome.db package required for this functionality")
  }
  dbi <- tryCatch(loadNamespace('AnnotationDbi'), error=function(e) NULL)
  if (is.null(dbi)) {
    stop("AnnotationDbi required")
  }

  # species <- resolve.species(species)
  si <- species_info(species)
  # The species info in reactome needs to look like:
  # Homo_sapiens, Mus_musculus, etc.
  species <- sub(" ", "_", si$species) #

  ## Find all KEGG pathways for the given species.
  ## Pathways are prefixed with the organism name like so:
  ##  `<GENUS> <SPECIES>: <PATHWAY NAME>`
  org.prefix <- paste0('^', sub('_', ' ', species), ': *')
  pathnames <- dbi$keys(rdb$reactome.db, keytype='PATHNAME')
  org.keep <- grepl(org.prefix, pathnames)
  org.pathnames <- pathnames[org.keep]

  info <- suppressMessages({
    ## Generates 1:many mapping because of ENTREZID
    dbi$select(rdb$reactome.db,
               columns=c('PATHID', 'PATHNAME', 'ENTREZID'),
               keys=org.pathnames,
               keytype='PATHNAME')
  })
  stopifnot(setequal(org.pathnames, info$PATHNAME))
  info <- info[!is.na(info[['ENTREZID']]),,drop=FALSE]
  setDT(info)

  N.pathid <- PATHNAME <- ENTREZID <- PATHID <- NULL # silence R CMD check NOTEs

  ## Are there multiple mappings for PATHNAME:PATHID combo?
  ## maybe from different organisms?
  u.id2name <- unique(info[, c('PATHID', 'PATHNAME'), with=FALSE])
  u.id2name[, N.pathid := .N, by='PATHID']
  if (nrow(dups <- u.id2name[N.pathid > 1])) {
    warning("Multiple PATHID to PATHNAME not resolved", immediate.=TRUE)
  }

  if (rm.species.prefix) {
    info[, PATHNAME := sub(org.prefix, '', PATHNAME)]
  }

  ## Somehow only one of the 'name's of the reactome genesets is encoded in
  ## UTF-8:
  ##   Loss of proteins required for interphase microtubule organization
  ##   from the centrosome
  ## All of the rest are "unknown". This causes annoying warnings when
  ## data.table tries to join on collection,name -- so I'm just nuking the
  ## encoding here.
  info <- info[, list(collection = "Reactome", name = PATHNAME,
                      feature_id = ENTREZID, pathId = PATHID)]
  Encoding(info$name) <- "unknown"
  gdb <- GeneSetDb(info)

  fn <- function(collection, name, gdb, ...) {
    if (missing(gdb) || !is(gdb, "GeneSetDb")) {
      return("https://reactome.org/")
    }
    gs <- geneSets(gdb)
    gset <- gs[gs[["collection"]] == "Reactome" & gs[["name"]] == name,]
    pathway_id <- gset[["pathId"]]
    sprintf("https://reactome.org/content/detail/%s", pathway_id)
  }
  geneSetCollectionURLfunction(gdb, "Reactome") <- fn

  featureIdType(gdb, "Reactome") <- EntrezIdentifier()

  org(gdb, "Reactome") <- species
  gdb
}
