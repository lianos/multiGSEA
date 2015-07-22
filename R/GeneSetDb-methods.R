##' (Re)-map geneset IDs to the rows in an expression object.
##'
##' This uses the already-set featureIdMap for the GeneSetDb
##'
##' @rdname conform
##'
##' @param x The GeneSetDb
##' @param y The expression object/matrix to conform to
##' @param unique.by If there are multiple rows that map to the identifiers
##'   used in the genesets, this is a means to pick the single row for that ID
##' @param min.gs.size Ensure that the genesets that make their way to the
##'   \code{GeneSetDb@@table} are of a minimum size
##' @param max.gs.size Ensure that the genesets that make their way to the
##'   \code{GeneSetDb@@table} are smaller than this size
##'
##' @return A \code{GeneSetDb} that has been matched/conformed to an
##'   expression object target (\code{y})
setMethod("conform", c(x="GeneSetDb"),
function(x, y, unique.by=c('none', 'mean', 'var'),
         min.gs.size=3L, max.gs.size=Inf, match.tolerance=0.25, ...) {
  unique.by <- match.arg(unique.by)
  if (unique.by != 'none') {
    stop("`unique.by` must be 'none' for now")
  }
  if (min.gs.size == 1) {
    stop("Doing GSEA with 1 gene doesn't really make sense, does it?")
  }
  if (min.gs.size < 3) {
    warning("Think twice before using genesets of size 2 in your GSEA",
            immediate.=TRUE)
  }
  if (max.gs.size < min.gs.size) {
    stop("max.gs.size must be larger than min.gs.size")
  }
  if (!any(sapply(.valid.x, function(claz) is(y, claz)))) {
    stop("Illegal type of expression object to conform to")
  }

  fm <- featureIdMap(x)
  fm$x.idx <- match(fm$x.id, rownames(y))
  fraction.match <- mean(!is.na(fm$x.idx))
  if (fraction.match <= match.tolerance) {
    warning("Fraction of gene set IDs that match rownames in the expression ",
            "object are low: ", sprintf("%.2f%% ", fraction.match * 100),
            immediate.=TRUE)
  }
  if (fraction.match == 0) {
    stop("None of the rownames of the expression object match the featureIds ",
         "for the genesets")
  }

  x@table <- x@db[, {
    f <- featureId
    xref <- fm[J(f)]
    n <- sum(!is.na(xref$x.idx))
    active <- n >= min.gs.size && n <= max.gs.size
    list(active=active, N=.N, n=n)
  }, by=c('collection', 'name')]

  inactive <- !x@table$active
  if (any(inactive)) {
    msg <- paste("Deactivating %d gene sets because conformation of GeneSetDb",
                 "to the target creates gene sets smaller than %s or greater",
                 "than %s")
    msg <- sprintf(msg, sum(inactive), as.character(min.gs.size),
                   as.character(max.gs.size))
    warning(msg, immediate.=TRUE)
  }

  x@featureIdMap <- fm
  x
})

setMethod("unconform", "GeneSetDb", function(x, ...) {
  x@table <- transform(x@table, active=FALSE)
  x@featureIdMap <- transform(x@featureIdMap, x.idx=NA_integer_)
  x
})


##' Check if GeneSetDb has been conformed (optionally to a specific
##' expression object)
##'
##' @export
##' @param x a \code{GeneSetDb}
##' @param to an optional expression object
##' @return logical indicating "conform"-ation status in general, or
##'   specifically to an expression object \code{to}
is.conformed <- function(x, to) {
  if (!is(x, 'GeneSetDb')) {
    stop("Only works on the GeneSetDb")
  }
  if (missing(to)) {
    ans <- any(!is.na(featureIdMap(x)$x.idx))
  } else {
    ## Verify that gsd is properly conformed to x
    fm <- subset(featureIdMap(x), !is.na(x.idx))
    ans <- nrow(fm) > 0 && all(rownames(to)[fm$x.idx] == fm$x.id)
  }
  ans
}

##' Creates a 1/0 matrix to indicate geneset membership to target object.
##'
##' rows are number of genesetes in \code{x}, columns are number of rows in
##' \code{y}
##'
##' @export
##'
##' @param x A \code{GeneSetDb}
##' @param y (optiona) A target (expression) object \code{x} is (or can be)
##'   conformed to
##' @return incidence matrix with nrows = number of genesets and columns are
##'   featureIDs. If \code{y} is passed in, the columns of the returned value
##'   match the rows of \code{y}.
incidenceMatrix <- function(x, y) {
  stopifnot(is(x, 'GeneSetDb'))
  gs <- NULL
  if (missing(y)) {
    val <- 'featureId'
    gs <- geneSets(x)
    all.fids <- unique(x@db$featureId)
    out <- matrix(0L, nrow=nrow(gs), ncol=length(all.fids),
                  dimnames=list(
                    paste(gs$collection, gs$name, sep=';'),
                    all.fids))
  } else {
    val <- 'x.idx'
    if (!inherits(y, .valid.x)) {
      stop("Invalid expression object (x) type: ", class(x)[1L])
    }
    if (!is.conformed(x, y)) {
      x <- conform(x, y)
    }
    gs <- geneSets(x)
    out <- matrix(0L, nrow(gs), nrow(y),
                  dimnames=list(
                    paste(gs$collection, gs$name, sep=';'),
                    rownames(y)))
  }

  for (i in 1:nrow(gs)) {
    fids <- featureIds(x, gs$collection[i], gs$name[i], value=val)
    out[i, fids] <- 1L
  }

  if (val == 'featureId') {
    out <- out[, colSums(out) > 0]
  }

  out
}


##' Interrogate "active" status of a given geneset.
##'
##' This method only works on one geneset at a time (ie. not vectorized)
##'
##' @export
##' @param x \code{GeneSetDb}
##' @param i collection of geneset
##' @param j name of geneset
##' @return logical indicating if geneset is active. throws an error if geneset
##'   does not exist in geneSets(x)
is.active <- function(x, i, j) {
  stopifnot(is(x, 'GeneSetDb'))
  idx <- .gsd.row.index(x, i, j)
  if (is.na(idx)) {
    stop(sprintf("Unknown geneset: (%s, %s)", i, j))
  }
  x@table$active[idx]
}

##' Fetch the IDs for the a given gene set.
##'
##' If the GeneSetDb \code{x} has been conformed to an expression object this
##' will default to return the featureId's as they are used/matched to the
##' expression object, otherwise it will return the featureIds used in the
##' definition of the gene set database.
##'
##' @rdname featureIds
##'
##' @param x The GeneSetDb
##' @param i The collection level identifier for the geneset
##' @param j The name level identifier for the genest
##' @param value What form do you want the id's in?
##'   \describe{
##'     \item{featureId}{the IDs used in the original geneset definitions}
##'     \item{x.id}{the ids of the features as they are used in the expression
##'           object}
##'     \item{x.idx}{The integer index into the expresion object \code{x} that
##'           the GeneSetDb has been conformed to.}
##'   }
##' @param fetch.all By default, this function only returns the IDs for the
##'   geneset that are matched in the target expression object if this
##'   \code{GeneSetDb} has been "conformed" (ie. \code{is.conformed(x) == TRUE}.
##'   Set this to \code{TRUE} if you want to get all featureIds irrespective of
##'   this constraint.
##'
##' @return A vector of identifiers (or indexes into an expression object) for
##'   the features in the given geneset. NA is returned if the geneset is not
##'   "active" (ie. listed in geneSets(x))
setMethod("featureIds", c(x="GeneSetDb"),
function(x, i, j, value, fetch.all=FALSE, active.only=is.conformed(x), ...) {
  if (missing(value)) {
    value <- if (is.conformed(x)) 'x.id' else 'featureId'
  }
  value <- match.arg(value, c('featureId', 'x.id', 'x.idx'))
  if (!isSingleCharacter(i)) {
    stop("collection (i) must be length 1 character vectors")
  }
  if (missing(j)) {
    whole.collection <- TRUE
  } else {
    if (!isSingleCharacter(j)) {
      stop("gene set name (j) must be length 1 character vectors")
    }
    whole.collection <- FALSE
  }

  gs <- geneSets(x, active.only=active.only)
  gs <- gs[, key(gs), with=FALSE]
  gs <- gs[J(i)]

  if (nrow(gs) == 0L) {
    stop("There are no ", if (active.only) "active " else NULL,
         "genesets in collection: ", i)
  }

  if (whole.collection) {
    db <- unique(x@db[gs], by='featureId')
  } else {
    ## I am purposefully not using `hasGeneSet` here for performance reasons
    ## hasGeneSet(x, i, j, as.error=TRUE)
    db <- x@db[J(i, j)]
    if (is.na(db$featureId[1L])) {
      msg <- sprintf("collection=%s, name=%s does not exist in GeneSetDb db",
                     i, j)
      stop(msg)
    }
  }

  fid.map <- merge(db, featureIdMap(x), by='featureId')
  if (is.conformed(x) && !fetch.all) {
    fid.map <- subset(fid.map, !is.na(x.idx))
  }

  fid.map[[value]]
})

setMethod("featureIdMap", c(x="GeneSetDb"), function(x) x@featureIdMap)

##' Replacing the featureIdMap blows away the "conform"-ation status of x
##'
##' This method ensures that there is only one featureId <-> x.id mapping value.
##' Note that this does not mean that this enforces a 1:1 mapping, it's just
##' that the same mapping is not listed more than once.
setReplaceMethod('featureIdMap', 'GeneSetDb', function(x, value) {
  if (!is.data.frame(value)) {
    stop("data.frame/table required for featureIdMap<-")
  }
  if (!ncol(value) == 2) {
    stop("featureIdMap must be a 2 column data.frame")
  }
  if (!all(x@db$featureId %in% value[[1]])) {
    stop("Some @db$featureId's are not in first column of new featureIdMap")
  }

  value <- as.data.table(value)
  setnames(value, c('featureId', 'x.id'))
  value <- unique(value, by=c('featureId', 'x.id'))
  setkeyv(value, 'featureId')

  x@featureIdMap <- value
  unconform(x)
})

setMethod("geneSets", c(x="GeneSetDb"),
function(x, active.only=is.conformed(x), ...) {
  if (active.only[1L]) x@table[active == TRUE] else x@table
})

##' @exportMethod geneSet
setMethod("geneSet", c(x="GeneSetDb"),
  function(x, i, j, active.only=is.conformed(x), fetch.all=FALSE, ...) {
  fids <- featureIds(x, i, j, value='featureId', active.only=active.only,
                     fetch.all=fetch.all, ...)
  info <- geneSets(x, active.only=FALSE)[J(i, j)]
  finfo <- featureIdMap(x)[J(fids)]
  cbind(info[rep(1, nrow(finfo))], finfo)
})

##' Subset GeneSetDb to only include specified genesets.
##'
##' This isn't exported yet because I don't like its implementation
##'
##' @param x \code{GeneSetDb}
##' @param keep logical vector as long as
##'   \code{nrow(geneSets(x, active.only=FALSE))}
##' @return a \code{GeneSetDb} that has only the results for the specified
##'   genesets.
subset.GeneSetDb <- function(x, keep) {
  stopifnot(is(x, 'GeneSetDb'))
  nr <- nrow(geneSets(x, active.only=FALSE))

  if (!is.logical(keep) && length(keep) != nr) {
    stop("The `keep` vector is FUBAR'd")
  }

  ## 1. Remove rows from x@table
  ## 2. Remove rows in x@db that belong to collection,name that do not exist
  ##    due to (1)
  ## 3. remove entries in x@featureIdMap for features that no longer exist in
  ##    updated db from (2)
  ## 4. Update x@collectionMetadata to:
  ##    a. remove all metadata for collections that are completely gone
  ##    b. update remaining collection,count entries for remaining collections

  ## 1
  keep.table <- x@table[keep]

  ## 2
  gs.keys <- keep.table[, key(keep.table), with=FALSE]
  setkeyv(gs.keys, key(keep.table))
  keep.db <- x@db[gs.keys, nomatch=0] ## only keep entries in db in gs.keys

  ## 3
  keep.featureIdMap <- subset(x@featureIdMap, featureId %in% keep.db$featureId)

  ## 4a
  keep.cm <- subset(x@collectionMetadata, collection %in% keep.db$collection)
  ## 4b
  cc <- keep.table[, list(name='count', value=.N), by='collection']
  setkeyv(cc, key(keep.cm))
  ## Currently (data.table v1.9.4( there's nothing I can do to make i.value a
  ## list element and this `set` mojog doesn't work either
  suppressWarnings(keep.cm[cc, value := list(i.value)])
  ## update.idxs <- keep.cm[cc, which=TRUE]
  ## val.idx <- which(colnames(keep.cm) == 'value')
  ## for (idx in seq_along(update.idxs)) {
  ##   set(keep.cm, update.idxs[idx], val.idx, list(cc$value[idx]))
  ## }

  out <- .GeneSetDb(table=keep.table,
                    db=keep.db,
                    featureIdMap=keep.featureIdMap,
                    collectionMetadata=keep.cm)
  out
}

## -----------------------------------------------------------------------------
## Functions over collections

##' Check if a collection exists in the \code{GeneSetDb}
##'
##' @export
##' @param x A \code{GeneSetDb}
##' @param collections character vector of name(s) of the collections to query
##' @param as.error logical if TRUE, this will error instead of returning FALSE
##' @return logical indicating if this collection exists
hasGeneSetCollection <- function(x, collections, as.error=FALSE) {
  stopifnot(is(x, 'GeneSetDb'))
  stopifnot(is.character(collections))

  ## I'm being paranoid
  ## if (!isTRUE(.validateCollectionMetadata(x))) {
  ##   stop("This GeneSetDb is not valid")
  ## }

  meta.idxs <- match(collections, collectionMetadata(x)$collection)
  gsc.exists <- !is.na(meta.idxs)
  if (!all(gsc.exists) && as.error) {
    bad <- paste("    * ", collections[!gsc.exists], collapse='\n', sep='')
    stop("The following collections to not exist:\n", bad)
  }

  gsc.exists
}

##' Check to see if the GeneSetDb has a collection,name GeneSet defined
##'
##' @export
##' @param x GeneSetDb
##' @param collection character indicating the collection
##' @param name character indicating the name of the geneset
##' @param as.error If \code{TRUE}, a test for the existance of the geneset
##'   will throw an error if the geneset does not exist
##' @return logical indicating whether or not the geneset is defined.
hasGeneSet <- function(x, collection, name, as.error=FALSE) {
  stopifnot(isSingleCharacter(collection) && isSingleCharacter(name))
  gs.exists <- !is.na(.gsd.row.index(x, collection, name))
  return(gs.exists)

  has.gsc <- hasGeneSetCollection(x, collection, as.error=as.error)
  if (!gs.exists && as.error) {
    stop(sprintf("geneset ('%s', '%s') does not exist", collection, name))
  }
  gs.exists
}

setMethod("collectionMetadata",
  c(x="GeneSetDb", collection="missing", name="missing"),
  function(x, collection, name) x@collectionMetadata)

setMethod("collectionMetadata",
  c(x="GeneSetDb", collection="character", name="missing"),
  function(x, collection, name) {
    stopifnot(isSingleCharacter(collection))
    hasGeneSetCollection(x, collection, as.error=TRUE)
    x@collectionMetadata[collection]
  })

setMethod("collectionMetadata",
  c(x="GeneSetDb", collection="character", name="character"),
  function(x, collection, name) {
    stopifnot(isSingleCharacter(collection))
    stopifnot(isSingleCharacter(name))
    hasGeneSetCollection(x, collection, as.error=TRUE)
    .col <- collection
    .name <- name
    cmd <- x@collectionMetadata
    idx <- cmd[J(.col, .name), which=TRUE]
    if (is.na(idx[1L])) {
      msg <- sprintf("metadata not defined for collection:%s, varname:%s",
                     .col, .name)
      stop(msg)
    }
    cmd$value[[idx]]
  })

setMethod("geneSetURL", c(x="GeneSetDb"), function(x, i, j, ...) {
  ## hasGeneSet(x, i, j, as.error=TRUE)
  ## fn <- geneSetCollectionURLfunction(x, i)
  ## fn(i, j)
  stopifnot(is.character(i) && is.character(j))
  stopifnot(length(i) == length(j))
  collections <- unique(i)
  col.exists <- hasGeneSetCollection(x, collections)
  url.fns <- Map(collections, col.exists, f=function(col, exists) {
    if (exists) {
      geneSetCollectionURLfunction(x, col)
    } else {
      function(x, y) NA_character_
    }
  })
  mapply(i, j, FUN=function(col, name) url.fns[[col]](col, name))
})

setMethod("geneSetCollectionURLfunction", "GeneSetDb", function(x, i, ...) {
  fn <- collectionMetadata(x, i, 'url_function')
  stopifnot(isSingleCharacter(i))
  fn.dt <- x@collectionMetadata[J(i, 'url_function'), nomatch=0]
  if (nrow(fn.dt) == 0) {
    ## validObject(x) : this call is slow and should never be FALSE anyway
    stop(sprintf("No url_function for collection '%s' found", i))
  }
  if (nrow(fn.dt) > 1) {
    ## validObject(x) : this call is slow and should never be FALSE anyway
    stop(sprintf("Multiple url_function defined for collection '%s' found", i))
  }

  fn <- fn.dt$value[[1L]]
  if (!is.function(fn)) {
    ## validObject(x) : this call is slow and should never be FALSE anyway
    stop(sprintf("The URL function for collection '%s' is missing", i))
  }

  fn
})

setReplaceMethod("geneSetCollectionURLfunction", "GeneSetDb",
function(x, i, value) {
  stopifnot(is.character(i))
  if (!is.function(value)) {
    stop("Function required")
  }
  if (length(formalArgs(value)) != 2L) {
    stop("URL function needs to take two arguments")
  }
  idx <- x@collectionMetadata[J(i, 'url_function'), which=TRUE]
  if (is.na(idx)) {
    ## This should never have happened, since x would fail validObject(x), but
    ## I'm being generous here
    add.me <- x@collectionMetadata[NA]
    add.me$collection[1L] <- i
    add.me$name[1L] <- 'url_function'
    add.me$value[[1L]] <- value
    x@collectionMetadata <- rbind(x@collectionMetadata, add.me)
  } else {
    x@collectionMetadata$value[[idx]] <- value
  }
  x
})


##' @importFrom BiocGenerics append
##' @export append
##' @export
setMethod("append", c(x='GeneSetDb'), function(x, values, after=NA) {
  if (!missing(after)) {
    warning("`after` argument is ignored in append,GeneSetDb")
  }
  if (!is(values, 'GeneSetDb')) {
    values <- GeneSetDb(values)
  }
  if (!is(values, 'GeneSetDb')) {
    stop("GeneSetDb expected by now")
  }

  ## Combine the db and featureIdMap(s)
  db <- rbindlist(list(x@db, values@db), use.names=TRUE, fill=TRUE)
  db <- unique(db, by=c('collection', 'name', 'featureId'))
  db <- setkeyv(db, key(x@db))

  fms <- list(featureIdMap(x), featureIdMap(values))
  fm <- rbindlist(fms, use.names=TRUE, fill=TRUE)
  ## ensure that a featureId entry maps to only one x.id entry
  ## DEBUG: Is this uniquification necessary?
  fm <- unique(fm, by=c('featureId', 'x.id'))
  fm[, x.idx := NA_integer_]  ## blow out any `conform`-ation information
  setkeyv(fm, 'featureId')

  cmeta <- rbind(x@collectionMetadata, values@collectionMetadata)
  cmeta <- unique(cmeta, by=key(x@collectionMetadata))
  setkeyv(cmeta, key(x@collectionMetadata))

  out <- .GeneSetDb(db=db, featureIdMap=fm, table=init.gsd.table.from.db(db),
                    collectionMetadata=cmeta)

  ## Transfer over any extra metadata (columns) of the @table slots from
  ## the two inputs incase the user stored extra data at the geneset level
  ## in them.
  gs <- rbindlist(list(x@table, values@table), use.names=TRUE, fill=TRUE)
  add.gs.cols <- setdiff(names(gs), names(out@table))
  if (length(add.gs.cols) > 0) {
    gs.keys <- key(out@table)
    new.table <- merge(out@table, gs[, c(gs.keys, add.gs.cols), with=FALSE],
                       by=gs.keys, all.x=TRUE)
    if (!all.equal(out@table, new.table[, names(out@table), with=FALSE])) {
      stop("There was a problem adding additional columns to geneset@table ",
           "during `append`")
    }
    out@table <- new.table
  }

  out
})

setMethod("nrow", "GeneSetDb", function(x) nrow(geneSets(x)))

##' Checks equality (feature parity) between GeneSetDb objects
##'
##' @method all.equal GeneSetDb
##' @export
all.equal.GeneSetDb <- function(target, current, features.only=FALSE, ...) {
  msg <- TRUE

  dbt <- setkeyv(copy(target@db), c('collection', 'name', 'featureId'))
  gst <- geneSets(target, active.only=FALSE)

  dbc <- setkeyv(copy(current@db), key(dbt))
  gsc <- geneSets(current, active.only=FALSE)

  proto <- new("GeneSetDb")
  if (features.only) {
    dbt <- dbt[, names(proto@db), with=FALSE]
    dbc <- dbc[, names(proto@db), with=FALSE]
    gst <- gst[, names(geneSets(proto)), with=FALSE]
    gsc <- gsc[, names(gst), with=FALSE]
  }

  msg <- all.equal(dbt, dbc)
  if (!isTRUE(msg)) {
    return(msg)
  }

  msg <- all.equal(gst, gsc)
  if (!isTRUE(msg) || features.only) {
    return(msg)
  }

  all.equal(collectionMetadata(target), collectionMetadata(current))
}

##' \code{as.list.GeneSetDb} and \code{as.expression.indexes} intentionally
##' have different default values for the \code{active.only} parameter.
##' \code{as.list} returns the ids of features, and \code{as.expression.indexes}
##' returns the integer index into the expression matrix that the
##' \code{GeneSetDb} is conformed to.
##'
##' @aliases as.list.GeneSetDb
##' @rdname GeneSetDb-conversion
##' @method as.list GeneSetDb
##' @export
as.list.GeneSetDb <- function(x, nested=FALSE, value=c('x.id', 'x.idx'),
                              active.only=is.conformed(x),
                              ...) {
  value <- match.arg(value)
  as.expression.indexes(x, value, active.only, nested)
}

##' Unrolls the GeneSetDb into a list of index vectors per "active" gene set
##'
##' @aliases as.expression.indexes
##' @rdname GeneSetDb-conversion
##'
##' @param x A \code{GeneSetDb}
##' @param value \code{x.idx} returns indexes into the conformed expression
##'   object as integers (ie. row index numbers) and \code{x.id} returns them
##'   as their featureId's as used in the target expression object.
##' @param active.only Only include "active" gene sets? Defaults to conformed
##'   status of the \code{x}
as.expression.indexes <- function(x, value=c('x.idx', 'x.id'),
                                  active.only=is.conformed(x),
                                  nested=FALSE) {
  value <- match.arg(value)
  if (!is(x, 'GeneSetDb')) {
    stop('GeneSetDb required')
  }
  if (value == 'x.idx' && !is.conformed(x)) {
    stop("GeneSetDb has not been conformed to an expression object")
  }
  gs <- geneSets(x, active.only=active.only)
  cats <- gs$collection
  nms <- gs$name
  out <- lapply(seq(cats), function(idx) {
    featureIds(x, cats[idx], nms[idx], value=value)
  })
  if (nested) {
    names(out) <- nms
    out <- split(out, cats)
  } else {
    names(out) <- paste(cats, nms, sep=';;')
  }
  out
}

## -----------------------------------------------------------------------------
## indexing

##' Identify the row number of the geneset being queried.
##'
##' The method only allows querying by single length character vectors
##'
##' @param x A GeneSetDb
##' @param i the collection of the geneset
##' @param j the name of the geneset
##'
##' @return The row index of the geneset under question. If i,j do not match
##' a given geneset, then NA is returned.
.gsd.row.index <- function(x, i, j) {
  stopifnot(is(x, 'GeneSetDb'))
  if (!(isSingleCharacter(i) && isSingleCharacter(j))) {
    stop("i (collection) and j (name) must both be length 1 character vectors")
  }
  ## if (!isTRUE(validObject(x))) {
  ##   stop("Invalid GeneSetDb")
  ## }
  geneSets(x, active.only=FALSE)[J(i, j), which=TRUE]
}

## setMethod("[", "GeneSetDb", function(x, i, j, ..., drop) {
##   if (length(list(...)) > 0L)
##     stop("invalid subsetting")
## })
