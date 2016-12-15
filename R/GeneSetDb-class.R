##' Creates a GeneSetDb from a variety of different types of inputs.
##'
##' @description
##' The GeneSetDb class serves the same purpose as the
##' \code{GSEABase::GeneSetCollection} class does: it acts as a centralized
##' object to hold collections of Gene Sets. The reason for its existence is
##' because there are things that I wanted to know about my gene set
##' collections that weren't easily inferred from what is essentially a
##' "list of GeneSets" that is the \code{GeneSetCollection} class.
##'
##' Gene Sets are internally represented by a \code{data.table} in "a tidy"
##' format, where we minimally require non \code{NA} values for the following
##' three \code{character} columns:
##'
##' \enumerate{
##'   \item{collection}
##'   \item{name}
##'   \item{featureId}
##' }
##'
##' The (\code{collection}, \code{name}) compound key is the primary key of
##' a gene set. There will be as many entries with the same
##' (\code{collection}, \code{name}) as there are genes/features in that set.
##'
##' The \code{GeneSetDb} tracks metadata about genesets at \emph{the collection}
##' level. This means that we assume that all of the \code{featureId}'s used
##' within a collection use the same type of feature identifier (such as
##' a \code{GSEABase::EntrezIdentifier()}), were defined in the same organism,
##' etc.
##'
##' \strong{Please refer to the "GeneSetDb" section of the vignette for more
##' details regarding the construction and querying of a \code{GeneSetDb}
##' object.}
##'
##' @section GeneSetDb Construction:
##'
##' The \code{GeneSetDb} constructor is sufficiently flexible enough to create
##' a \code{GeneSetDb} object from a variety of formats that are commonly used
##' in the bioconductor echosystem, such as:
##'
##' \enumerate{
##'   \item{GSEABase::GeneSetCollection}{
##'     If you've already got a \code{GeneSetCollection} on your hands, you can
##'     simply pass it to the \code{GeneSetDb} constructor.
##'   }
##'   \item{list of ids}{
##'     This format is commonly used to define gene sets in the edgeR/limma
##'     universe for testing with camera, roast, romer, etc. The names of
##'     the list items are the gene set names, and their values are a character
##'     vector of gene identifiers. When it's a single list of lists, you must
##'     provide a value for \code{collectionName}. You can embed multiple
##'     collections of gene sets by having a three-deep list-of-lists-of-ids.
##'     The top level list define the different collections, the second level
##'     are the genesets, and the third level are the feature identifiers for
##'     each gene set. See the examples for clarification.
##'   }
##'   \item{data.frame-like object}{
##'     To keep track of your own custom gene sets, you have probably realized
##'     the importance of maintaing your own sanity, and likely have gene sets
##'     organized in a table like object that has something like
##'     the \code{collection}, \code{name}, and \code{featureId} required for
##'     a \code{GeneSetDb}. Simply rename the appropriate columns to the ones
##'     prescribed here, and pass that into the constructor. Any other
##'     additional columns (symbol, perhaps(?)) in the table will be copied
##'     into the \code{GeneSetDb}.
##'   }
##' }
##'
##' @section Interrogating a GeneSetDb:
##'
##' You might wonder what gene sets are defined in a \code{GeneSetDb}: see
##' the \code{\link{geneSets}} function.
##'
##' Curious about what features are defined in your \code{GeneSetDb}? See
##' the \code{\link{featureIds}} function.
##'
##' Want the details of a particular gene set? Try the \code{\link{geneSet}}
##' function. This will return a \code{data.frame} of the gene set definition.
##' Calling \code{\link{geneSet}} on a \code{\link{MultiGSEAResult}} will
##' return the same \code{data.frame} along with the differential expression
##' statistics for the individual members of the geneSet across the contrast
##' that was tested in the \code{\link{multiGSEA}} call that created the
##' \code{\link{MultiGSEAResult}}.
##'
##' @rdname GeneSetDb-class
##' @aliases GeneSetDb
##' @export
##'
##' @param x A \code{GeneSetCollection}, a "two deep" list of either
##'   \code{GeneSetCollection}s or lists of character vectors, which are
##'   the gene identifers. The "two deep" list represents the different
##'   collections (top level) at the top level, and each such list is a named
##'   list itself, which represents the gene sets in the given collection.
##' @param featureIdMap A data.frame with  2 character columns. The first
##'   column is the ids of the genes (features) used to identify the genes in
##'   \code{gene.sets}, the second second column are IDs that this should be
##'   mapped to. Useful for testing probelevel microarray data to gene level
##'   geneset information.
##' @param collectionName If \code{x} represents a singular collection, ie.
##'   a single \code{GeneSetCollection} or a "one deep" (named (by geneset))
##'   list of genesets, then this parameter provides the name for the
##'   collection. If \code{x} is multiple collections, this can be character
##'   vector of same length with the names. In all cases, if a collection name
##'   can't be defined from this, then collections will be named anonymously.
##'   If a value is passed here, it will overide any names stored in the list of
##'   \code{x}.
##'
##' @examples
##' ## Some Gene Sets from Rooney et al. Note that the feature id vectors
##' ## do not have to be named, they just are here to add clarity.
##'
##' ## list of ids
##' gs.list <- list(
##'   'cytolytic activity'=c("GZMA"="3001", "PRF1"="5551"),
##'   'NK cells'=c("KLRC1"="3821", "KLRF1"="51348"),
##'   'MHC Class I'=c("B2M"="567", "HLA-A"="3105", "TAP1"="6890"))
##' gdb1 <- GeneSetDb(gs.list, collectionName='rooney_hacohen')
##'
##' ## list of lists of ids
##' gs.lol <- list(
##'   rooney_hacohen=gs.list,
##'   angelova=list(
##'     'Activated B Cells'=c('AKNA'='80709', 'ARHGAP25'='9938', 'CCL21'='6366'),
##'     'DC'=c('C1QC'='714', 'CCDC88A'='55704', 'CCL13'='6357')))
##' gdb2 <- GeneSetDb(gs.lol)
##'
##' ## Using a data.frame for input
##' library('dplyr')
##' gs.df <- lapply(names(gs.lol), function(coll.name) {
##'   lapply(names(gs.lol[[coll.name]]), function(gs.name) {
##'     tibble(collection=coll.name, name=gs.name,
##'            featureId=gs.lol[[coll.name]][[gs.name]],
##'            ## symbol is optional, but helpful
##'            symbol=names(gs.lol[[coll.name]][[gs.name]]))
##'   }) %>% bind_rows
##' }) %>% bind_rows
##' gdb3 <- GeneSetDb(gs.df)
##'
##' ## GeneSetDb Interrogation
##' (gsets <- geneSets(gdb3))
##' (nkcells <- geneSet(gdb3, 'rooney_hacohen', 'NK cells'))
##' (fids <- featureIds(gdb3))
GeneSetDb <- function(x, featureIdMap=NULL, collectionName=NULL) {
  gdb <- if (is(x, 'GeneSetCollection')) {
    GeneSetDb.GeneSetCollection(x, featureIdMap, collectionName)
  } else if (is(x, 'data.frame')) {
    GeneSetDb.data.frame(x, featureIdMap, collectionName)
  } else if (is(x, 'list')) {
    GeneSetDb.list(x, featureIdMap, collectionName)
  } else {
    stop("No GeneSetDb constructor defined for: ", class(x)[1L])
  }

  proto <- new("GeneSetDb")
  setkeyv(gdb@db, key(proto@db))
  gdb
}

GeneSetDb.data.frame <- function(x, featureIdMap=NULL, collectionName=NULL) {
  stopifnot(is.data.frame(x) && nrow(x) > 0)
  x <- setDT(as.data.frame(copy(x)))
  if (!'collection' %in% names(x)) {
    if (!is.character(collectionName) &&
        !length(collectionName) %in% c(1L, nrow(x))) {
      stop("If no `collection` column is provided in `x`, ",
           "collectionName must be well defined")
    }
    x[, collection := collectionName]
  }
  req.cols <- c('collection', 'name', 'featureId')
  cols.missed <- setdiff(req.cols, names(x))
  if (length(cols.missed)) {
    stop("The following columns are missing from `x`:\n ",
         paste(cols.missed, collapse=", "))
  }

  x <- unique(x, by=req.cols)
  lol <- sapply(unique(x[['collection']]), function(col) {
    with(x[collection == col], split(featureId, name))
  }, simplify=FALSE)
  gdb <- GeneSetDb(lol, featureIdMap=featureIdMap)

  add.cols <- setdiff(names(x), req.cols)
  if (length(add.cols)) {
    db <- merge(gdb@db, x, by=req.cols, all.x=TRUE)
    db0 <- setkeyv(copy(gdb@db), req.cols)
    setkeyv(db, req.cols)
    if (!all.equal(db0, db[, req.cols, with=FALSE])) {
      warning("Something unexpected happened merging more feature metadata",
              immediate.=TRUE)
    }
  }
  gdb
}

GeneSetDb.list <- function(x, featureIdMap=NULL, collectionName=NULL) {
  if (!is.list(x) || length(x) == 0L) {
    stop("A non-empty list is required for this function")
  }
  if (is(x[[1]], 'GeneSetCollection')) {
    if (is.null(collectionName)) {
      collectionName <- names(x)
    }
    out <- GeneSetDb.list.of.GeneSetCollections(x, featureIdMap, collectionName)
    return(out)
  }

  proto <- new("GeneSetDb")

  ## Is this just a "one deep" list of genesets? If so, let's wrap it in
  ## a list
  if (is.single.list.of.feature.vectors(x)) {
    x <- list(x)
  }
  if (is.null(collectionName)) {
    collectionName <- names(x)
  }
  if (is.null(collectionName)) {
    collectionName <- sprintf('anon_collection_%d', seq(x))
  }
  if (!is.character(collectionName)) {
    stop("Character vector expected for `collectionName`")
  }
  if (length(collectionName) != length(x)) {
    stop("length(collectionName) != length(x)")
  }
  names(x) <- collectionName

  db <- init.gsd.db.from.list.of.lists(x)
  tbl <- init.gsd.table.from.db(db)

  meta <- tbl[, {
    list(name=c('count', "url_function"),
         value=list(.N, function(x, y) NA_character_))
  }, by='collection']
  setkeyv(meta, c('collection', 'name'))

  if (is.null(featureIdMap)) {
    .ids <- unique(db$featureId)
    featureIdMap <- data.table(featureId=.ids, x.id=.ids, x.idx=NA_integer_)
    setkeyv(featureIdMap, key(featureIdMap(proto, .external=FALSE)))
  }

  out <- .GeneSetDb(table=tbl,
                    db=db,
                    featureIdMap=featureIdMap,
                    collectionMetadata=meta)
  out
}

##' @importFrom GSEABase setName geneIds geneIdType
GeneSetDb.list.of.GeneSetCollections <- function(x, featureIdMap=NULL,
                                                 collectionName=names(x)) {
  stopifnot(is.list(x))
  stopifnot(length(x) > 0)
  stopifnot(all(sapply(x, is, 'GeneSetCollection')))
  if (is.null(collectionName)) {
    collectionName <- sprintf('anon_collection_%d', seq(x))
  }
  if (!is.character(collectionName) && length(collectionName) != length(x)) {
    stop("Invalid value for `collectionName`")
  }

  lol <- lapply(1:length(x), function(i) {
    gsc.name <- collectionName[i]
    gsc <- x[[i]]
    id.list <- lapply(gsc, geneIds)
    org <- unique(sapply(gsc, GSEABase::organism))
    id.type <- unique(sapply(gsc, function(x) class(geneIdType(x))))
    if (length(org) > 1) {
      warning("multiple organisms defined in geneset collection: ", gsc.name,
              immediate.=TRUE)
    }
    if (length(id.type) > 1) {
      stop("different idtypes used in genesets: ", paste(id.type, collapse=','))
    }
    setNames(id.list, sapply(gsc, setName))
  })
  names(lol) <- collectionName
  GeneSetDb.list(lol, featureIdMap, collectionName)
}


GeneSetDb.GeneSetCollection <- function(x, featureIdMap=NULL,
                                        collectionName='anon_collection') {
  stopifnot(is.character(collectionName) && length(collectionName) == 1)
  gsc.list <- setNames(list(x), collectionName)
  GeneSetDb.list.of.GeneSetCollections(gsc.list, featureIdMap, collectionName)
}


##' @importFrom GSEABase GeneSetCollection GeneSet
setAs("GeneSetDb", "GeneSetCollection", function(from) {
  gs <- geneSets(from, .external=FALSE)
  n.coll <- length(unique(gs$collection))

  ## We explicitly set the key type after subsetting here in the event that
  ## a 0 row data.table is returned -- this isn't keyed in 1.9.4, which seems
  ## like a bug
  id.type <- subset(collectionMetadata(from, .external=FALSE),
                    name == 'id_type')
  setkeyv(id.type, 'collection')
  org <- subset(collectionMetadata(from, .external=FALSE),
                name == 'organism')
  setkeyv(org, 'collection')

  gsl <- lapply(1:nrow(gs), function(i) {
    name <- gs$name[i]
    coll <- gs$collection[i]
    idt <- id.type[coll]$value[[1]]
    ids <- featureIds(from, coll, name, 'featureId')
    xorg <- org[coll]$value[[1]]
    if (is.null(xorg)) {
      xorg <- ""
    }
    set.name <- name
    if (n.coll > 1L) {
      set.name <- paste0(coll, ';', set.name)
    }
    GeneSet(ids, setName=set.name, geneIdType=idt, organism=xorg)
  })
  gsc <- GeneSetCollection(gsl)
  gsc
})

setAs("GeneSetCollection", "GeneSetDb", function(from) {
  GeneSetDb(from)
})


## Constructor Helper Functions ------------------------------------------------
init.gsd.db.from.list.of.lists <- function(x) {
  proto <-.GeneSetDb()
  ## Ensure x is list of list of geneset features
  x <- validate.gene.sets.input(x)

  ## "melt" the x list-of-lists Create internal db data.table that
  ## stores the "pristine" geneset membership information passed into this
  ## function via `x`.
  db <- local({
    groups <- lapply(names(x), function(g.name) {
      members <- lapply(names(x[[g.name]]), function(id) {
        ids <- unique(x[[g.name]][[id]])
        data.table(collection=g.name, name=id, featureId=ids)
      })
      rbindlist(members)
    })
    setkeyv(rbindlist(groups), key(proto@db))
  })

  db
}

init.gsd.table.from.db <- function(db) {
  proto <- .GeneSetDb()
  out <- db[, list(active=FALSE, N=.N, n=NA_integer_), by=c('collection', 'name')]
  setkeyv(out, key(proto@table))
}


setMethod("show", "GeneSetDb", function(object) {
  proto <- .GeneSetDb()
  msg <- paste("GeneSetDb with %d defined genesets across %d collections",
               "(%d gene sets are active)")
  msg <- sprintf(msg,
                 nrow(unique(object@db, by=key(proto@db))),
                 length(unique(object@db$collection)),
                 sum(object@table$active))
  is.conf <- paste("  Conformed:", ifelse(is.conformed(object), 'yes', 'no'))
  hr <- paste(rep("=", nchar(msg)), collapse='')
  hr.sub <- gsub('=', '-', hr)
  cat(hr, "\n", msg, "\n", is.conf, "\n", hr.sub, "\n", sep="")
  data.table:::print.data.table(geneSets(object, .external=FALSE))
  cat(hr.sub, "\n", msg, "\n", is.conf, "\n", hr, "\n", sep="")
})

## -----------------------------------------------------------------------------
## Functions to check validity of GeneSetDb
setValidity("GeneSetDb", function(object) {
  proto <- .GeneSetDb()

  ## Get classes of the slots
  sclass <- sapply(slotNames(proto), function(x) class(slot(proto, x))[1L])

  ## Check all data.tables to ensure that they have a superset of the columns
  ## to the comparable prototype versions *and* they share the same defined
  ## keys.
  dt.errs <- sapply(names(sclass)[sclass == 'data.table'], function(s) {
    check.dt(slot(object, s), slot(proto, s))
  }, simplify=FALSE)

  u <- unlist(dt.errs)
  if (is.character(u)) {
    return(u)
  }

  cm.errs <- .validateCollectionMetadata(object)
  if (!isTRUE(cm.errs)) {
    return(cm.errs)
  }

  ## ---------------------------------------------------------------------------
  ## Further check @db slot:
  ## 1. ensure all features in @db have a row in the @featureIdMap
  if (!all(object@db$featureId %in% featureIdMap(object, .external=FALSE)$featureId)) {
    return("Some @db$featureId's are not in featureIdMap(object)$featureId")
  }
  if (any(is.na(object@db$featureId))) {
    return("NA's not permitted in @db$featureId")
  }
  ## Ensure that the collection,id combination is unique in @table
  if (any(duplicated(object@table, by=key(proto@table)))) {
    return("Duplicated gene set entries in @table")
  }

  ## Ensure that that collection,id,featureId are unique in @db
  if (any(duplicated(object@db, by=c('collection', 'name', 'featureId')))) {
    return("Duplicated collection,id,featureId in @db")
  }

  TRUE
})

##' Checks validaty of collectionMetadata of a GeneSetDb
##'
##' This function is for internal use only
##'
##' The following assertions are tested:
##' \itemize{
##'   \item{All collection,name entries are unique}
##'   \item{All collections have a url_function}
##'   \item{
##'     The collections that are listed in collectionMetadata have >= 1 defined
##'     genesets in \code{geneSets(object)}
##'   }
##' }
##'
##' @param object A \code{GeneSetDb}
##' @return TRUE if the collectionMetadata is kosher, otherwise a character
##'   vector of errors.
.validateCollectionMetadata <- function(object) {
  ## ---------------------------------------------------------------------------
  ## Check the collectionMetadata bits
  ## ---------------------------------------------------------------------------
  ##
  ## 1. Ensure that all collection,name entries are unique
  dupd <- duplicated(object@collectionMetadata, by=c('collection', 'name'))
  if (any(dupd)) {
    return('Duplicated (collection,name) entries in @collectionMetadata')
  }
  ## 2. Collect information about required metadata entries for each collection,
  ##    ie. the count of genesets in the collection and their url_function
  cm.info <- object@collectionMetadata[, {
    is.url.fn <- which(name == 'url_function')
    is.count <- which(name == 'count')
    if (length(is.url.fn) == 0) {
      url.fn.status <- 'not-defined'
    } else {
      fn <- value[[is.url.fn]]
      if (!is.function(fn)) {
        url.fn.status <- 'not-a-function'
      } else if (length(formalArgs(fn)) != 2) {
        url.fn.status <= 'function-not-2-args'
      } else {
        url.fn.status <- 'ok'
      }
    }
    if (length(is.count) == 0) {
      count <- NA_integer_
    } else {
      count <- value[[is.count]]
    }
    list(count=count, url.fn.stauts=url.fn.status)
  }, keyby='collection']
  ## 3. Ensure url functions are kosher
  bad.fns <- cm.info$url.fn.status != 'ok'
  if (any(bad.fns)) {
    msg <- paste('bad url fns:\n',
                 paste(sprintf('  %s:%s', cm.info$collection[bad.fns],
                               cm.info$url.fn.status[bad.fns]),
                       collapse='\n'))
    return(msg)
  }
  ## 4. Ensure gene set counts match
  gs.info <- geneSets(object, active.only=FALSE, .external=FALSE)[, {
    list(count=.N)
  }, by='collection']
  setkeyv(gs.info, 'collection')
  if (nrow(gs.info) != nrow(cm.info)) {
    msg <- paste('Number of defined collections in geneSets() does not match',
                 'defined collections in collectionMetadata')
    return(msg)
  }
  if (!all.equal(gs.info[, list(collection, count)],
                 cm.info[, list(collection, count)])) {
    msg <- paste('gene set counts per collection do not match')
    return(msg)
  }

  TRUE
}

## -----------------------------------------------------------------------------
## Helper functions to enable setting up of the GeneSetDb


is.single.list.of.feature.vectors <- function(x) {
  is.list(x) && all(sapply(x, is.character))
}

validate.gene.sets.input <- function(gene.sets) {
  ## Did the user only enter a single list of character vectors? We need to
  ## change this into a list of lists
  if (is.single.list.of.feature.vectors(gene.sets)) {
    ## make this into a list of lists
    gene.sets <- list(undef=gene.sets)
  }

  ## is each top level entry a list?
  top.are.lists <- sapply(gene.sets, is.list)
  if (!all(top.are.lists)) {
    stop("The gene.sets input should be a list of lists")
  }

  ## Are these specified as characters?
  xxx <- unlist(gene.sets)
  if (!is.character(xxx)) {
    stop("Identifiers used in gene.set list must be characters")
  }

  groups <- names(gene.sets)
  if (!is.character(groups) || any(duplicated(groups))) {
    stop("names() of gene.sets list must be set and unique")
  }

  bad.gs <- !sapply(gene.sets, is.single.list.of.feature.vectors)
  if (any(bad.gs)) {
    report <- paste(head(which(bad.gs), 10), collapse=',')
    if (sum(bad.gs) > 10) {
      report <- paste0(report, ',...')
    }
    stop("These gene.set lists are bad. Are IDs characters(?): ", report)
  }
  gene.sets
}

setMethod("updateObject", "GeneSetDb",
function(object, ..., verbose=FALSE) {
  return(object)
  ## Unfortunately we've internally generated GeneSetDb objects that weren't
  ## entirely properly keyed.
  proto <- new("GeneSetDb")
  ekeys <- key(proto@db)
  xkeys <- key(object@db)
  if (!setequal(ekeys, xkeys) || !all.equal(ekeys, xkeys)) {

    setkeyv(object@db, ekeys)
  }
  object
})
