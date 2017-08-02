##' Annotates rows of a data.frame with geneset membership from a GeneSetDb
##'
##' This is helpful when you don't have a monsterly sized GeneSetDb. There will
##' be as many new columns added to \code{x} as there are active genesets in
##' \code{gdb}.
##'
##' @export
##'
##' @param x A data.frame with genes/features in rows
##' @param gdb A \code{GeneSetDb} object with geneset membership
##' @param x.ids The name of the column in \code{x} that holds the feautre
##'   id's in \code{x} that match the featureId's in \code{gdb}, or a vector
##'   of id's to use for each row in \code{x} that represent these.
##' @return Returns the original \code{x} with additional columns: each is a
##'   logical vector that indicates membership for genesets defined in
##'   \code{gdb}.
##'
##' @examples
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- getMSigGeneSetDb('h', 'human')
##' mg <- multiGSEA(gdb, vm, vm$design, 'tumor', methods=NULL)
##' lfc <- logFC(mg)
##' annotated <- annotateGeneSetMembership(lfc, gdb, 'featureId')
##'
##' ## Show only genes that are part of 'HALLMARK_ANGIOGENESIS' geneset
##' (angio <- subset(annotated, `h;;HALLMARK_ANGIOGENESIS`))
annotateGeneSetMembership <- function(x, gdb, x.ids=NULL, ...) {
  stopifnot(is(x, 'data.frame'))
  stopifnot(is(gdb, 'GeneSetDb'))

  if (is.null(x.ids)) {
    ## Guess the column of x that has featureIds by identifying the column of
    ## x that has the highest percent membership in gdb@featureIdMap
    membership <- sapply(x, function(col) {
      if (!is.character(col)) 0 else mean(col %in% featureIdMap(gdb)$x.id)
    })
    x.ids <- names(x)[which.max(membership)]
  }
  if (!is.character(x.ids)) {
    stop("character expected for x.ids")
  }
  if (length(x.ids) == 1L && is.character(x[[x.ids]])) {
    ## This is a column in df that has the gene IDs that `x` expects
    x.ids <- x[[x.ids]]
  }
  if (nrow(x) != length(x.ids)) {
    stop("x.ids must be a vector of gene ids, or a column name in x of these")
  }

  im <- incidenceMatrix(gdb, x.ids, ...)
  storage.mode(im) <- 'logical'
  cbind(x, t(im))
}

setMethod("length", "GeneSetDb", function(x) nrow(geneSets(x)))

##' (Re)-map geneset IDs to the rows in an expression object.
##'
##' @description \code{conform}-ing, a \code{GeneSetDb} to a target expression
##' object is an important step required prior to perform any type of GSEA. This
##' function maps the featureIds used in the GeneSetDb to the elements of a
##' target expression object (ie. the rows of an expression matrix, or the
##' elements of a vector of gene-level statistics).
##'
##' After \code{conform}-ation, each geneset in the \code{GeneSetDb} is flagged
##' as active (or inactive) given the number of its features that are
##' successfully mapped to \code{target} and the minimum and maximum number of
##' genes per geneset required as specified by the \code{min.gs.size} and
##' \code{max.gs.size} parameters, respectively.
##'
##' Only genesets that are marked with \code{active=TRUE} will be used in any
##' downstream gene set operations.
##'
##' @section Related Functions:
##' \enumerate{
##'   \item{\code{unconform}}{Resets the conformation mapping.}
##'   \item{\code{is.conformed}}{
##'     If \code{to} is missing, looks for evidence if \code{coform} has
##'     been called (at all) on \code{x}. If \code{to} is provided, specifically
##'     checks that \code{x} has been conformed to the target object \code{to}.
##'   }
##' }
##'
##' @rdname conform
##'
##' @param x The GeneSetDb
##' @param target The expression object/matrix to conform to. This could also
##'   just be a character vector of IDs.
##' @param unique.by If there are multiple rows that map to the identifiers used
##'   in the genesets, this is a means to pick the single row for that ID
##' @param min.gs.size Ensure that the genesets that make their way to the
##'   \code{GeneSetDb@@table} are of a minimum size
##' @param max.gs.size Ensure that the genesets that make their way to the
##'   \code{GeneSetDb@@table} are smaller than this size
##' @param match.tolerance Numeric value between [0,1]. If the fraction of
##'   \code{featureId}s used in \code{x} that match \code{rownames(y)} is below
##'   this number, a warning will be fired.
##' @param ... moar args
##'
##' @return A \code{GeneSetDb} that has been matched/conformed to an expression
##'   object target (\code{y})
##'
##' @examples
##' es <- exampleExpressionSet()
##' gdb <- exampleGeneSetDb()
##' head(geneSets(gdb))
##' gdb <- conform(gdb, es)
##' ## Note the updated values `active` flag, and n (the number of features
##' ## mapped per gene set)
##' head(geneSets(gdb))
setMethod("conform", c(x="GeneSetDb"),
function(x, target, unique.by=c('none', 'mean', 'var'),
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
  ## We are allowing y to be a character vector of featureIds here
  if (is.vector(target) && is.character(target)) {
    target <- matrix(1L, nrow=length(target), ncol=1L,
                     dimnames=list(target, NULL))
  }
  if (!any(sapply(.valid.x, function(claz) is(target, claz)))) {
    stop("Illegal type of expression object to conform to")
  }
  # x <- copy(x)
  fm <- featureIdMap(x, .external=FALSE)
  fm$x.idx <- match(fm$x.id, rownames(target))
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
    xref <- fm[list(f)]
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

##' @export
##' @rdname conform
is.conformed <- function(x, to) {
  if (!is(x, 'GeneSetDb')) {
    stop("Only works on the GeneSetDb")
  }
  if (missing(to)) {
    ans <- any(!is.na(featureIdMap(x, .external=FALSE)$x.idx))
  } else {
    ## Verify that gsd is properly conformed to x
    fm <- subset(featureIdMap(x, .external=FALSE), !is.na(x.idx))
    to.ids <- if (is.character(to))
      to
    else if (is.vector(to))
      names(to)
    else
      rownames(to)
    if (!is.character(to.ids)) {
      stop("featureIds unsuccessfully extracted from `to`")
    }
    ans <- nrow(fm) > 0 && all(to.ids[fm$x.idx] == fm$x.id)
  }
  ans
}

##' Creates a 1/0 matrix to indicate geneset membership to target object.
##'
##' Generates an inidcator matrix to indicate membership of genes (columns)
##' to gene sets (rows). If \code{y} is provided, then the incidence is mapped
##' across the entire feature-space of \code{y}.
##'
##' @export
##'
##' @param x A \code{GeneSetDb}
##' @param y (optional) A target (expression) object \code{x} is (or can be)
##'   conformed to
##' @return incidence matrix with nrows = number of genesets and columns are
##'   featureIDs. If \code{y} is passed in, the columns of the returned value
##'   match the rows of \code{y}.
##'
##' @examples
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- getMSigGeneSetDb('h', 'human')
##' im <- incidenceMatrix(gdb)
##' imv <- incidenceMatrix(gdb, vm)
incidenceMatrix <- function(x, y, ...) {
  stopifnot(is(x, 'GeneSetDb'))
  gs <- NULL
  if (missing(y)) {
    val <- 'featureId'
    ynames <- unique(x@db$featureId)
    ncol <- length(ynames)
  } else {
    val <- 'x.idx'
    x <- conform(x, y, ...)
    if (is.vector(y)) {
      ynames <- y
      ncol <- length(y)
    } else {
      ynames <- rownames(y)
      ncol <- nrow(y)
    }
  }

  gs <- geneSets(x, .external=FALSE)
  dimnames <- list(paste(gs$collection, gs$name, sep=';;'), ynames)
  out <- matrix(0L, nrow(gs), ncol, dimnames=dimnames)

  for (i in 1:nrow(gs)) {
    fids <- featureIds(x, gs$collection[i], gs$name[i], value=val)
    out[i, fids] <- 1L
  }

  out
}


##' Interrogate "active" status of a given geneset.
##'
##' Returns the \code{active} status of genesets, which are specified by
##' their collection,name compound keys. This function is vectorized and
##' supports query of multiple gene sets at a time. If a requested
##' collection,name gene set doesn't exist, this throws an error.
##'
##' @export
##' @param x \code{GeneSetDb}
##' @param i collection of geneset(s)
##' @param j name of geneset(s) (must be same length as \code{i})
##' @return logical indicating if geneset is active. throws an error if
##'   any requested geneset does not exist in \code{x}.
is.active <- function(x, i, j) {
  stopifnot(is(x, 'GeneSetDb'))
  stopifnot(is.character(i), is.character(j), length(i) == length(j))
  gsx <- geneSets(x, active.only=FALSE, .external=FALSE)[list(i,j), nomatch=NA]
  res <- gsx$active
  isna <- is.na(res)
  if (any(isna)) {
    unk <- paste(i[isna], j[isna], sep=":", collapse=",")
    stop(sprintf("Unknown genesets: %s", unk))
  }
  res
}

setMethod("subsetByFeatures", c(x="GeneSetDb"),
function(x, features, value=c('featureId', 'x.id', 'x.idx'), ...) {
  value <- match.arg(value)
  ## some good old data.table voodoo going on inside here
  unk.f <- setdiff(features, featureIds(x, value=value))
  if (length(unk.f)) {
    warning(length(unk.f), "/", length(features), " do not exist in GeneSetDb")
    features <- setdiff(features, unk.f)
  }

  dat <- merge(x@db, featureIdMap(x, .external=FALSE), by='featureId')
  hits <- unique(dat[dat[[value]] %in% features, list(collection, name)])
  gs.all <- geneSets(x, active.only=FALSE, .external=FALSE)
  keep <- rep(FALSE, nrow(gs.all))
  gs.idx <- gs.all[hits, which=TRUE]
  if (length(gs.idx)) {
    keep[gs.idx] <- TRUE
  }
  x[keep]
})

##' @rdname featureIds
setMethod("featureIds", c(x="GeneSetDb"),
function(x, i, j, value=c('featureId', 'x.id', 'x.idx'),
         active.only=is.conformed(x), ...) {
  if (missing(value)) {
    value <- if (is.conformed(x)) 'x.id' else 'featureId'
  }
  value <- match.arg(value, c('featureId', 'x.id', 'x.idx'))

  if (missing(i) && missing(j)) {
    ## User isn't asking about any particular collection, but just wants all
    ## features in the GeneSetDb as a whole ... OK(?)
    fmap <- x@featureIdMap
    if (is.conformed(x) && active.only) {
      fmap <- fmap[!is.na(x.idx)]
    }
    out <- unique(fmap, by=value)[[value]]
    return(out[!is.na(out)])
  }

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

  gs <- geneSets(x, active.only=active.only, .external=FALSE)
  gs <- gs[, key(gs), with=FALSE]
  gs <- gs[list(i)]

  if (nrow(gs) == 0L) {
    stop("There are no ", if (active.only) "active " else NULL,
         "genesets in collection: ", i)
  }

  if (whole.collection) {
    db <- unique(x@db[gs], by='featureId')
  } else {
    ## I am purposefully not using `hasGeneSet` here for performance reasons
    ## hasGeneSet(x, i, j, as.error=TRUE)
    db <- x@db[list(i, j)]
    if (is.na(db$featureId[1L])) {
      msg <- sprintf("collection=%s, name=%s does not exist in GeneSetDb db",
                     i, j)
      stop(msg)
    }
  }

  fid.map <- featureIdMap(x, .external=FALSE)[db$featureId]
  if (is.conformed(x) && active.only) {
    fid.map <- fid.map[!is.na(x.idx)]
  }

  fid.map[[value]]
})

setMethod("featureIdMap", c(x="GeneSetDb"), function(x, .external=TRUE) {
  ret.df(x@featureIdMap, .external=.external)
})

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
  value[, x.idx := NA_integer_]
  x@featureIdMap <- value
  unconform(x)
})

##' @rdname geneSets
setMethod("geneSets", c(x="GeneSetDb"),
function(x, active.only=is.conformed(x), ... , .external=TRUE) {
  out <- if (active.only[1L]) x@table[active == TRUE] else x@table
  ret.df(out, .external=.external)
})

##' @rdname geneSet
setMethod("geneSet", c(x="GeneSetDb"),
function(x, i, j, active.only=is.conformed(x), with.feature.map=FALSE, ...,
         .external=TRUE) {
  stopifnot(isSingleCharacter(i), isSingleCharacter(j))
  fids <- featureIds(x, i, j, value='featureId', active.only=active.only, ...)
  info <- geneSets(x, active.only=FALSE, .external=FALSE)[list(i, j)]
  info <- info[, c("collection", "name", "active", "N", "n"), with=FALSE]

  ## Fetch information from x@db. Extra information per feature are stored here
  dbx <- x@db[list(i,j,fids)]
  out <- cbind(info[rep(1L, nrow(dbx))], dbx[, -(1:2), with=FALSE])

  ## Add featureIdMap info
  if (with.feature.map) {
    fminfo <- featureIdMap(x, .external=FALSE)[list(fids)]
    out <- cbind(out, fminfo[, -1L, with=FALSE])
  }

  ret.df(out, .external=.external)
})

##' Subset GeneSetDb to only include specified genesets.
##'
##' This is a utility function that is called by \code{[.GeneSetDb} and is not
##' exported because it is not meant for external use.
##'
##' DEBUG: If `keep` is all FALSE, this will explode. What does an empty
##' GeneSetDb look like, anyway? Something ...
##'
##' @param x \code{GeneSetDb}
##' @param keep logical vector as long as
##'   \code{nrow(geneSets(x, active.only=FALSE))}
##' @return a \code{GeneSetDb} that has only the results for the specified
##'   genesets.
subset.GeneSetDb <- function(x, keep) {
  stopifnot(is(x, 'GeneSetDb'))
  if (all(keep == FALSE)) {
    stop("Cannot subset GeneSetDb down to empty (`keep` is all FALSE)")
  }
  nr <- nrow(geneSets(x, active.only=FALSE, .external=FALSE))

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
  ## list element and this `set` mojo doesn't work either
  value <- NULL # silence R CMD check NOTEs
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

##' Subset whole genesets from a GeneSetDb
##' @exportMethod [
setMethod("[", "GeneSetDb", function(x, i, j, ..., drop=FALSE) {
  if (length(list(...)) > 0) {
    stop("parameters in '...' not supported")
  }
  if (!missing(drop)) {
    warning("drop argument ignored", immediate.=TRUE)
  }
  if (!missing(j)) {
    warning("j is ignored", immediate.=TRUE)
  }
  if (!is.logical(i) && length(i) != nrow(geneSets(x))) {
    stop("i must be a logical vector as long as nrow(geneSets(x))")
  }
  subset.GeneSetDb(x, i)
})


if (FALSE) {
## @exportMethod subset
## @importFrom BiocGenerics subset
setMethod("subset", "GeneSetDb", function(x, subject, select, drop=FALSE, ...) {
  data.table:::subset.data.table(geneSets(x, .external=FALSE), substitute(subject))
})
}

## -----------------------------------------------------------------------------
## Functions over collections

##' Check if a collection exists in the \code{GeneSetDb}
##'
##' @export
##' @param x A \code{GeneSetDb}
##' @param collection character vector of name(s) of the collections to query
##' @param as.error logical if TRUE, this will error instead of returning FALSE
##' @return logical indicating if this collection exists
hasGeneSetCollection <- function(x, collection, as.error=FALSE) {
  stopifnot(is(x, 'GeneSetDb'))
  stopifnot(is.character(collection))
  meta.idxs <- match(collection,
                     collectionMetadata(x, .external=FALSE)$collection)
  gsc.exists <- !is.na(meta.idxs)
  if (!all(gsc.exists) && as.error) {
    bad <- paste("    * ", collection[!gsc.exists], collapse='\n', sep='')
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
##'
##' @examples
##' gdb <- exampleGeneSetDb()
##' hasGeneSet(gdb, c('c2', 'c7'), c('BIOCARTA_AGPCR_PATHWAY', 'something'))
hasGeneSet <- function(x, collection, name, as.error=FALSE) {
  gs.exists <- !is.na(.gsd.row.index(x, collection, name))
  if (as.error && !all(gs.exists)) {
    msg <- "The follwing %d gene sets do not exist: \n%s"
    not <- paste(collection[!gs.exists], name[!gs.exists],
                 sep=":", collapse="\n")
    stop(sprintf(msg, sum(!gs.exists), not))
  }
  gs.exists
}

setMethod("collectionMetadata",
  c(x="GeneSetDb", collection="missing", name="missing"),
  function(x, collection, name, .external=TRUE) {
    ret.df(x@collectionMetadata, .external=.external)
  })

setMethod("collectionMetadata",
  c(x="GeneSetDb", collection="character", name="missing"),
  function(x, collection, name, .external=TRUE) {
    stopifnot(isSingleCharacter(collection))
    hasGeneSetCollection(x, collection, as.error=TRUE)
    ret.df(x@collectionMetadata[collection], .external=.external)
  })

setMethod("collectionMetadata",
  c(x="GeneSetDb", collection="character", name="character"),
  function(x, collection, name, .external=TRUE) {
    stopifnot(isSingleCharacter(collection))
    stopifnot(isSingleCharacter(name))
    hasGeneSetCollection(x, collection, as.error=TRUE)
    .col <- collection
    .name <- name
    cmd <- x@collectionMetadata
    idx <- cmd[list(.col, .name), which=TRUE]
    if (is.na(idx[1L])) {
      msg <- sprintf("metadata not defined for collection:%s, varname:%s",
                     .col, .name)
      stop(msg)
    }
    cmd$value[[idx]]
  })

setMethod("geneSetURL", c(x="GeneSetDb"), function(x, i, j, ...) {
  stopifnot(is.character(i), is.character(j), length(i) == length(j))
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
  fn.dt <- x@collectionMetadata[list(i, 'url_function'), nomatch=0]
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
  valid <- function(v) {
    if (!isTRUE(is.function(v))) return(FALSE)
    if (length(formalArgs(v)) != 2L) {
      ## "URL function needs to take two arguments"
      return(FALSE)
    }
    TRUE
  }
  addCollectionMetadata(x, i, 'url_function', value, valid)
})

# setReplaceMethod("collectionUrlFunction", "GeneSetDb", function(x, i, value) {
#   valid <- function(v) {
#     if (!isTRUE(is.function(v))) return(FALSE)
#     if (length(formalArgs(v)) != 2L) {
#       ## "URL function needs to take two arguments"
#       return(FALSE)
#     }
#     TRUE
#   }
#   addCollectionMetadata(x, i, 'url_function', value, valid)
# })

setReplaceMethod("featureIdType", "GeneSetDb", function(x, i, value) {
  valid <- function(v) is(v, 'GeneIdentifierType')
  addCollectionMetadata(x, i, 'id_type', value, valid)
})

setMethod("featureIdType", "GeneSetDb", function(x, i, ...) {
  x@collectionMetadata[list(i, 'id_type')]$value[[1L]]
})

setReplaceMethod("org", "GeneSetDb", function(x, i, value) {
  valid <- function(v) {
    if (!is.character(v) && length(v) == 1L) {
      return(FALSE)
    }
    info <- strsplit(v, '_')[[1L]]
    if (length(info) != 2L) {
      return(FALSE)
    }
    known <- c('Mus_musculus', 'Homo_sapiens')
    if (!v %in% known) {
      warning("Value '",v,"' for organism is unrecognized (but will be used).",
              v, immediate.=TRUE)
    }
    TRUE
  }
  if (missing(i)) {
    colls <- unique(x@collectionMetadata$collection)
    for (coll in colls) {
      org(x, coll) <- value
    }
    return(x)
  }
  addCollectionMetadata(x, i, 'organism', value, valid)
})

setMethod("org", "GeneSetDb", function(x, i, ...) {
  if (missing(i)) {
    x@collectionMetadata[name == 'organism']
  } else {
    x@collectionMetadata[list(i, 'organism')]$value[[1L]]
  }
})

##' @section Adding arbitrary collectionMetadata:
##'
##' Adds arbitrary metadata to a gene set collection of a GeneSetDb
##'
##' Note that this is not a replacement method! You must catch the returned
##' object to keep the one with the updated `collectionMetadata`. Although this
##' function is exported, I imagine this being used mostly through predefined
##' replace methods that use this as a utility function, such as the replacement
##' methods for \code{\link{org}}, and \code{\link{featureIdType}}.
##'
##' \preformatted{
##'   gdb <- getMSigGeneSetDb('h')
##'   gdb <- addCollectionMetadata(gdb, 'h', 'foo', 'bar')
##' }
##'
##' @export
##' @rdname collectionMetadata
##'
##' @param x \code{GeneSetDb}
##' @param xcoll The collection name
##' @param xname The name of the metadata variable
##' @param value The value of the metadata variable
##' @param validate.value.fn If a function is provided, it is run on
##'   \code{value} and msut return \code{TRUE} for addition to be made
##' @param allow.add If \code{FALSE}, this xcoll,xname should be in the
##'   \code{GeneSetDb} already, and this will fail because something is
##'   deeply wrong with the world
##' @return The updated \code{GeneSetDb}
addCollectionMetadata <- function(x, xcoll, xname, value,
                                  validate.value.fn=NULL, allow.add=TRUE) {
  stopifnot(is(x, 'GeneSetDb'))
  stopifnot(is.character(xcoll))
  if (!hasGeneSetCollection(x, xcoll)) {
    stop("GeneSetDb does not have collection: ", xcoll)
  }
  if (is.function(validate.value.fn)) {
    if (!isTRUE(validate.value.fn(value))) {
      stop(sprintf("Invalid value used to update %s,%s", xcoll, xname))
    }
  }

  ## update or replace
  if (!is.list(value)) {
    value <- list(value)
  }
  idx <- x@collectionMetadata[list(xcoll, xname), which=TRUE]
  if (is.na(idx)) {
    if (!allow.add) {
      msg <- sprintf("%s,%s metadata does not exist in the GeneSetDb, but it",
                     "should be there. Your GeneSetDb is hosed")
      stop(msg)
    }
    ## the variable you want to ender here is not there yet, so let's put it in
    add.me <- x@collectionMetadata[NA]
    add.me$collection[1L] <- xcoll
    add.me$name[1L] <- xname
    add.me$value[[1L]] <- value
    x@collectionMetadata <- setkeyv(rbind(x@collectionMetadata, add.me),
                                    key(x@collectionMetadata))
  } else {
    x@collectionMetadata$value[idx] <- value
  }
  x
}

##' @importMethodsFrom BiocGenerics append
##' @exportMethod append
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

  fms <- list(featureIdMap(x, .external=FALSE),
              featureIdMap(values, .external=FALSE))
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
    new.table <- merge(out@table,
                       gs[, c(gs.keys, add.gs.cols), with=FALSE],
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
##' @export
##' @method all.equal GeneSetDb
##' @param target The reference \code{GeneSetDb} to compare against
##' @param current The \code{GeneSetDb} you wan to compare
##' @param features.only Only compare the "core" columns of \code{target@@db}
##'   and \code{target@@table}. It is possible that you added additional columns
##'   (to keep track of symbols in \code{target@@db}, for instance) that you
##'   want to ignore for the purposes of the equality test.
##' @param ... moar args.
##' @return \code{TRUE} if equal, or \code{character} vector of messages if not.
all.equal.GeneSetDb <- function(target, current, features.only=FALSE, ...) {
  stopifnot(is(target, 'GeneSetDb'))
  stopifnot(is(current, 'GeneSetDb'))

  msg <- TRUE

  dbt <- setkeyv(copy(target@db), c('collection', 'name', 'featureId'))
  gst <- geneSets(target, active.only=FALSE, .external=FALSE)

  dbc <- setkeyv(copy(current@db), key(dbt))
  gsc <- geneSets(current, active.only=FALSE, .external=FALSE)

  proto <- new("GeneSetDb")
  if (features.only) {
    dbt <- dbt[, names(proto@db), with=FALSE]
    dbc <- dbc[, names(proto@db), with=FALSE]
    gst <- gst[, names(geneSets(proto, .external=FALSE)), with=FALSE]
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

  all.equal(collectionMetadata(target, .external=FALSE),
            collectionMetadata(current, .external=FALSE))
}

##' Convert a GeneSetDb to other formats.
##'
##' @description
##' As awesome as a GeneSetDb is, you might find a time when you'll need your
##' gene set information in an other format. To do that, we provide the
##' following functions:
##'
##' \enumerate{
##'   \item{\code{as.data.frame}}{
##'     Perhaps the most natural format to convert to in order to save locally
##'     an examine outside of Bioconductor's GSEA universe, but not many other
##'     tools accet gene set definitions in this way
##'   }
##'   \item{\code{as.list}}{
##'     A named list of feature identifiers. This is the format that many of
##'     the limma gene set testing methods use
##'   }
##'   \item{\code{as(gdb, 'GeneSetCollection')}}{
##'     The \code{\link[GSEABase]{GeneSetCollection}} class.
##'   }
##' }
##'
##' The \code{as.*} functions accept a \code{value} parameter which indicates
##' the type of IDs you want to export in the conversion. The following choices
##' are available:
##'
##' \enumerate{
##'   \item{\code{"featureId"}}{
##'     The ID used as originally entered into the \code{GeneSetDb}
##'   }
##'   \item{\code{"x.idx"}}{
##'     Only valid if the GeneSetDb \code{x} has been \code{conform}-ed to
##'     an expession container. This option will export the features as the
##'     integer rows of the expression container.
##'   }
##'   \item{\code{"x.id"}}{
##'     Only valid if the GeneSetDb \code{x} has been \code{conform}-ed. The
##'     target expression container might use feature identifiers that are
##'     different than what is in the GeneSetDb. If an active featureMap is
##'     set on the GeneSetDb, this will convert the original feature identifiers
##'     into a different target space (entrez to ensembl, for instance). Using
##'     this option, the features will be provided in the target space.
##'   }
##' }
##'
##' @export
##' @rdname conversion
##' @name conversion
##' @aliases as.data.frame as.list
##' @method as.data.frame GeneSetDb
##'
##' @param x A \code{GeneSetDb} object
##' @param value The value type to export for the feature ids
##' @param active.only If the \code{GeneSetDb} is conformed, do you want to only
##'   return the features and genests that match target and are "active"?
##' @param ... nothing
##' @return a converted \code{GeneSetDb}
##'
##' @examples
##' es <- exampleExpressionSet()
##' gdb <- conform(exampleGeneSetDb(), es)
##' gdf <- as.data.frame(gdb)
##' gdfi <- as.data.frame(gdb, 'x.idx')
##' gdl <- as.list(gdb)
as.data.frame.GeneSetDb <- function(x, row.names=NULL, optional=FALSE,
                                    value=c('featureId', 'x.id', 'x.idx'),
                                    active.only=is.conformed(x), ...) {
  stopifnot(is(x, 'GeneSetDb'))
  value <- match.arg(value)
  if (!is.conformed(x) && value %in% c('x.id', 'x.idx')) {
    stop("must use value='featureId' for non-conformed GeneSetDb'")
  }

  fid.map <- featureIdMap(x, .external=FALSE)
  if (is.conformed(x) && active.only) {
    # x.idx <- NULL # silence R CMD check NOTEs
    fid.map <- fid.map[!is.na(x.idx)]
  }

  gs <- copy(geneSets(x, active.only=active.only, .external=FALSE))

  gene2cat <- merge(x@db, fid.map, by='featureId')
  gene2cat <- gene2cat[!is.na(gene2cat[[value]]),]
  gene2cat$finalId <- gene2cat[[value]]

  # collection <- name <- finalId <- featureId <- NULL # silence R CMD check NOTEs
  out <- gene2cat[, list(collection, name, featureId=finalId)]

  more.cols <- setdiff(names(gene2cat),
                       c(names(out), names(fid.map), 'finalId'))
  if (length(more.cols)) {
    for (col in more.cols) {
      out[, (col) := gene2cat[[col]]]
    }
  }

  setkeyv(out, key(x@db))
  trimmed <- out[gs[, list(collection, name)]]
  setDF(trimmed)
}

## Split and conserve ordering
##
## Using base::split will turn f into a factor and won't preserve the ordering
## of the elements in f. This function preserves the split order in f
csplit <- function(x, f) {
  f <- as.character(f)
  ff <- factor(f, levels=unique(f))
  split(x, ff)
}

##' @rdname conversion
##' @method as.list GeneSetDb
##' @export
as.list.GeneSetDb <- function(x, value=c('featureId', 'x.id', 'x.idx'),
                              active.only=is.conformed(x), nested=FALSE,
                              ...) {
  stopifnot(is(x, 'GeneSetDb'))
  value <- match.arg(value)
  df <- as.data.frame(x, value=value, active.only=active.only)
  if (nested) {
    colls <- unique(df$collection)
    ## Using the "naive split" call converts xdf$name to a factor and doesn't
    ## preserve ordering
    out <- sapply(colls, function(coll) {
      xdf <- df[df[['collection']] == coll,,drop=FALSE]
      csplit(xdf$featureId, xdf$name)
    }, simplify=FALSE)
  } else {
    df$key <- paste(df$collection, df$name, sep=";;")
    out <- csplit(df$featureId, df$key)
  }
  out
}

##' @importFrom GSEABase GeneSetCollection GeneSet NullIdentifier
setAs("GeneSetDb", "GeneSetCollection", function(from) {
  gs <- geneSets(from, .external=FALSE)
  n.coll <- length(unique(gs$collection))

  ## We explicitly set the key type after subsetting here in the event that
  ## a 0 row data.table is returned -- this isn't keyed in 1.9.4, which seems
  ## like a bug
  cmd <- collectionMetadata(from, .external=FALSE)
  org <- subset(cmd, name == 'organism')
  id.type <- subset(cmd, name == 'id_type')
  setkeyv(id.type, 'collection')
  setkeyv(org, 'collection')

  gsl <- lapply(1:nrow(gs), function(i) {
    name <- gs$name[i]
    coll <- gs$collection[i]
    idt <- id.type[coll]$value[[1]]
    if (!is(idt, 'GeneIdentifierType')) {
      idt <- NullIdentifier()
    }
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

## -----------------------------------------------------------------------------
## indexing

##' Identify the row number of the geneset being queried.
##'
##' The method only allows querying by single length character vectors
##'
##' @rdname gsd.row.index
##' @param x A GeneSetDb
##' @param i the collection of the geneset
##' @param j the name of the geneset
##'
##' @return The row index of the geneset under question. If i,j do not match
##' a given geneset, then NA is returned.
.gsd.row.index <- function(x, i, j) {
  stopifnot(is(x, 'GeneSetDb'), is.character(i), is.character(j))
  geneSets(x, active.only=FALSE, .external=FALSE)[list(i, j), which=TRUE]
}
