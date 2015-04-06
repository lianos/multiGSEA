##' Create a GeneSetDb to match an input expression object.
##'
##' @section Easily mapping rownames(x) to IDs used in gene.sets:
##'
##' It may be that the IDs used in the provided \code{gene.sets} are different
##' than the ones used in \code{rownames(x)}, for instance the IDs in
##' \code{gene.sets} might (are likely to) be Entrez IDs, and the rownames of
##' \code{x} might by affy probe set IDs, what then?
##'
##' This is where the \code{mapping} parameter becomes useful.
##'
##' In all cases, the \code{featureIds} stored in the \code{GeneSetDb} are
##' "in the space" of the expression object \code{x} the GeneSetDb has been
##' "conformed" to. To recover the original IDs used, you could reference the
##' \code{featureIdMap} table.
##'
##' TODO: document the GeneSetDb ID mapping more thoroughly.
##'
##' @import GSEABase
##' @importFrom matrixStats rowVars
##' @export
##'
##' @param x A list (of lists) of character vectors that specify
##'   the features in a particular geneset. The top level of the list (of lists)
##'   are the collection ids. The names of the inner lists define the id for the
##'   particular geneset.
##' @param featureIdMap A data.frame with  2 character columns. The first
##'   column is the ids of the genes (features) used to identify the genes in
##'   \code{gene.sets}, the second second column are IDs that this should be
##'   mapped to. Useful for testing probelevel microarray data to gene level
##'   geneset information.
##' @param x The expression values measured. Could be a matrix, ExpressionSet,
##'   etc. rownames are required on this "thing" for it to work.
##' @param min.gs.size The minimum number of features in a geneset required
##'   for testing (this is calculated after gene.sets to expression feature
##'   mapping).
##' @param max.gs.size Same as above, but is upper limit of geneset size for
##'   testing.
GeneSetDb <- function(x, featureIdMap=NULL) {
  if (!is.list(x)) {
    stop("Do not know how to handle gene.set object of class: ", class(x1))
  }

  proto <- new("GeneSetDb")

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
    setkeyv(featureIdMap, key(featureIdMap(proto)))
  }

  out <- .GeneSetDb(table=tbl,
                    db=db,
                    featureIdMap=featureIdMap,
                    collectionMetadata=meta)
  out
}

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
  msg <- paste("GeneSetDb with %d defined genesets across %d categories",
               "(%d gene sets are active)")
  msg <- sprintf(msg,
                 nrow(unique(object@db, by=key(proto@db))),
                 length(unique(object@db$collection)),
                 sum(object@table$active))
  is.conf <- paste("  Conformed:", ifelse(is.conformed(object), 'yes', 'no'))
  hr <- paste(rep("=", nchar(msg)), collapse='')
  hr.sub <- gsub('=', '-', hr)
  cat(hr, "\n", msg, "\n", is.conf, "\n", hr.sub, "\n", sep="")
  data.table:::print.data.table(geneSets(object))
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
  if (!all(object@db$featureId %in% featureIdMap(object)$featureId)) {
    return("Some @db$featureId's are not in featureIdMap(object)$featureId")
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
  gs.info <- geneSets(object, active.only=FALSE)[, {
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
