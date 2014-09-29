##' A class that holds a data.table with geneset membership.
##'
##' @exportClass GeneSetTable
##'
##' @slot table The data.table with geneset information
##' @slot feature.lookup Maps the ids used in the geneset lists to the ids
##' (rows) over the expression data the GSEA is run on
##' @slot species Either "human" or "mouse" -- only (really) used when the object
##' was built from a call to MSigDB
setClass("GeneSetTable",
         slots=c(
           table="data.table",
           feature.lookup="data.table",
           species="character"),
         prototype=prototype(
           table=data.table(
             group=character(), id=character(), N=integer(),
             n=integer(), membership=as.list(logical())),
           feature.lookup=data.table(
             gset.id=character(),
             x.id=character(),
             x.index=integer()),
           species=character()))

##' Create a GeneSetTable to match an input expression object.
##'
##' @importFrom matrixStats rowVars
##' @export
##'
##' @param gene.sets This can be several things: (a) a character vector that
##' lists the MSigDB set id's (in conjunction with the \code{species}
##' parameter); (b) a list of logical vectors into the rows of \code{x}. This
##' should be a list of lists, where the top-level lists identified the "group"
##' its child elements (genesets) belong to (like the c1, c2, etc. grouping)
##' in the MSigDB lists; ie. the provenance of the gene set (c) An already
##' rigged up GeneSetTable;
##' @param x The expression values measured. Could be a matrix, ExpressionSet,
##' etc. rownames are required.
##' @param mapping A data.frame with two columns. The first column is the ids
##' used to identify the genes in \code{gene.sets}, the second column are the
##' rows of x that match the particular gene.set id. This is minimally meant to
##' support GSEA over microarray datsets, where the geneset IDs are entrez IDs,
##' but \code{rownames(x)} are probeset IDs.
GeneSetTable <- function(gene.sets, x, mapping=NULL, min.gs.size=5,
                         stop.on.duplicate.x.xref=TRUE,
                         unique.by=c('mean', 'var'),
                         species='human',
                         version=.wehi.msigdb.current) {
  if (is(gene.sets, 'GeneSetTable')) {
    return(conform(gene.sets, x))
  }

  species <- match.species(species)
  if (is.character(gene.sets)) {
    gene.sets <- getMSigDBset(gene.sets, species, version)
  }

  if (!is.list(gene.sets)) {
    stop("We need a list of named gene.set vectors to create a GeneSetTable")
  }

  is.msigdb.lookup <- is.character(gene.sets)
  ## gene.sets can be either a list of index vectors
  ## a list of lists of index vectors
  ## a character vector indicating the MSigDB gene set to use
  if (is.msigdb.lookup) {
    gene.sets <- getMSigDBset(gene.sets, species=species)
  }

  if (is.null(mapping)) {
    mapping <- data.frame(rownames(x), rownames(x), stringsAsFactors=FALSE)
  }
  lookup <- buildMappingTable(mapping, x)


  if (!is(gene.sets, 'list')) {
    stop("Do not know how to handle gene.set object of class: ",
         class(gene.sets[1L]))
  }
  if (is.list.of.index.vectors(gene.sets)) {
    ## make this into a list of lists
    gene.sets <- list(undef=gene.sets)
  }
  bad.gs <- !sapply(gene.sets, is.list.of.index.vectors)
  if (sum(bad.gs)) {
    stop("Bad gene.sets elements: ", paste(which(bad.gs), collapse=','))
  }
  dt <- list.of.geneset.lists.to.data.table(x, gene.sets, lookup)
  out <- new("GeneSetTable",
             table=dt,
             feature.lookup=lookup,
             species=if (is.msigdb.lookup) species else character())

  if (min.gs.size > 0) {
    out@table <- out@table[n >= min.gs.size]
  }

  if (any(duplicated(out@feature.lookup$x.id))) {
    out <- resolveDuplicateIdMapping(out, x, unique.by=unique.by)
  }

  out
}

resolveDuplicateIdMapping <- function(gst, x, unique.by=c('mean', 'var')) {
  ## TODO: deal with multiplie gene set IDs mapping to rows in x
  warning("Need to resolve duplicate mapping from geneset IDs to rows in x",
          immediate.=TRUE)
  ## This will just run through the @table$membership vectors and set the
  ## elements that will be removed to FALSE
}

##' Builds a table that maps gene set IDs to the data in the expression object x
##'
##' @param mapping A two-column data.frame, with the first column being the
##' id's of the things used in the genesets, and the second column being the
##' id's that these things are mapped to in \code{x}
##' @param x The "rectangular thing" that holds the epxression data.
##'
##' @return a data.table with the following columns, keyed on "gset.id"
##' \enumerate{
##'   \item gset.id The ids used in the genesets
##'   \item x.id The ids of the same things in \code{x}
##'   \item x.index The row index of "the thing" in \code{x}
##' }
buildMappingTable <- function(mapping, x, expr.id.fn=rownames) {
  if (!is.data.frame(mapping) || ncol(mapping) != 2) {
    stop("`mapping` must be a two-column data.frame of characters, see ",
         "?GeneSetTable")
  }
  if (!all(sapply(mapping, is.character))) {
    stop("Both columns of `mapping` must be character")
  }
  if (!is.character(rownames(x))) {
    stop("x must be an object with rownames, this should not have gotten here")
  }
  mapping <- as.data.table(mapping)
  setnames(mapping, c('gset.id', 'x.id'))
  mapping[, x.index := match(x.id, expr.id.fn(x))]
  mapping <- mapping[is.na(x.index) == FALSE]
  setkeyv(mapping, 'gset.id')
}

##' @importFrom Biobase featureNames
setMethod("featureNames", c(object='GeneSetTable'), function(object) {
  fn <- character(max(object@feature.lookup$x.index))
  fn[object@feature.lookup$x.index] <- object@feature.lookup$x.id
  fn
})

##' Conforms x to y.
##'
##' @exportMethod conform
setGeneric("conform", function(x, ...) standardGeneric("conform"))

setMethod("conform", c(x='GeneSetTable'),
function(x, y, mapping=x@feature.lookup, unique.by=c('mean', 'var'),
         expr.id.fn=rownames, ...) {
  y <- validateInputs(y)$x

  ## Will update x@feature.lookup to match to the rows in y and return a
  ## GeneSetTable that's ready to be analyzed against `data`

  ## Build a list of lists out of this thing, and rip it through GeneSetTable
  x.ids <- featureNames(x)
  groups <- unique(x@table$group)
  lol <- sapply(groups, function(g) {
    these <- x@table[group == g]
    out <- lapply(1:nrow(these), function(idx) {
      x.ids[these$membership[[idx]]]
    })
    names(out) <- these$id
    out
  }, simplify=FALSE)

  ## mapping <- data.frame(x@feature.lookup$gset.id, x@feature.lookup$x.id)
  ## lookup <- buildMappingTable(mapping, y, expr.id.fn=expr.id.fn)

  gst <- GeneSetTable(lol, y, unique.by=unique.by, min.gs.size=1)
  xref <- match(paste(gst@table$group, gst@table$id, sep='.'),
                paste(x@table$group, x@table$id, sep='.'))
  gst@table[, N := x@table$N[xref]]
  gst
})

setMethod('show', 'GeneSetTable', function(object) {
  cat("GeneSetTable with XX genesets across YY groups\n")
  cat("----------------------------------------------\n")
  data.table:::print.data.table(object@table)
})

setValidity("GeneSetTable", function(object) {
  proto <- new('GeneSetTable')

  ## Check lookup table
  kosher.lookup <- .check.dt.columns(object@feature.lookup,
                                     proto@feature.lookup)
  if (!isTRUE(kosher.lookup)) {
    return(kosher.lookup)
  }
  na.xrows <- is.na(object@feature.lookup$x.index)
  if (any(na.xrows)) {
    return("There are NA row indices into the eXpression object")
  }

  kosher.data <- .check.dt.columns(object@table, proto@table)
  if (!isTRUE(kosher.data)) {
    return(kosher.data)
  }
  not.logical.idx <- !sapply(object@table$membership, is, 'logical')
  if (n.bad <- sum(not.logical.idx)) {
    msg <- "%d%s list elements are not logical"
    some.bad <- head(which(not.logical.idx, 5))
    msg <- sprintf(msg, some.bad, if (n.bad > 5) '...' else '')
    return(msg)
  }

  lens <- sapply(object@table$membership, length)
  ulens <- unique(lens)
  if (length(ulens) != 1) {
    return("Not all membership index vectors are the same length")
  }

  TRUE
})

## -----------------------------------------------------------------------------
## Non exported utility functions for GeneSetTable stuff.
is.list.of.index.vectors <- function(x) {
  if (!is(x, 'list')) {
    return(FALSE)
  }
  all(sapply(x, function(xx) is.vector(xx) && !is.list(xx)))
}

##' Helper function to create data.table out of list-of-list genesets
##'
##' @param x An "expression thing"
##' @param gene.sets A list of lists of genesets into \code{x}
list.of.geneset.lists.to.data.table <- function(x, gene.sets, lookup) {
  groups <- lapply(names(gene.sets), function(g.name) {
    indexes <- lapply(gene.sets[[g.name]], function(v) {
      index.vector.to.xrow(x, v, lookup)
    })
    data.table(group=g.name,
               id=names(gene.sets[[g.name]]),
               N=sapply(indexes, '[[', 'N'),
               n=sapply(indexes, '[[', 'n'),
               membership=lapply(indexes, '[[', 'membership'))
  })
  out <- rbindlist(groups)
  out
}

index.vector.to.xrow <- function(x, index, lookup) {
  if (!is.vector(index)) {
    stop("index needs to be a vector")
  }
  if (is.logical(index)) {
    if (length(index) != nrow(x)) {
      stop("Length of logical index vector does not match nrows in x")
    }
    bool <- index
    N <- sum(index)
  } else {
    N <- length(index)
    if (is.character(index)) {
      idx.dt <- data.table(gset.id=index, key='gset.id')
      xref <- lookup[idx.dt]$x.index
      nomatch <- is.na(xref)
      if (sum(nomatch)) {
        ## msg <- sprintf("%d (%.2f%%) of named genes do not match x",
        ##                sum(nomatch), mean(nomatch) * 100)
        ## warning(msg, immediate.=TRUE)
      }
      index <- xref[!nomatch]
    }
    if (is.numeric(index)) {
      if (all(as.integer(index) != index)) {
        stop("Only integer index vectors are allowed")
      }
      index <- as.integer(index)
      if (min(index) < 1 || max(index) > nrow(x)) {
        stop("Numeric index vector has out of bounds vectors")
      }
    } else {
      stop("Unknown type for index vector: ", class(index)[1L])
    }

    bool <- logical(nrow(x))
    bool[index] <- TRUE
  }

  list(N=N, n=sum(bool), membership=bool)
}


.check.dt.columns <- function(x, prototype) {
  if (!is.data.table(x)) {
    return('Object is not a data.table')
  }
  if (!setequal(names(x), names(prototype))) {
    bad.cols <- setdiff(names(x), names(prototype))
    msg <- sprintf("Illegal columns in x: ", paste(bad.cols, collapse=','))
    return(msg)
  }

  invalid.cols <- sapply(names(prototype), function(col) {
    proto.class <- class(prototype[[col]])[1L]
    dat <- x[[col]]
    is.null(dat) || !is(dat, proto.class)
  })

  if (any(invalid.cols)) {
    bad.cols <- names(prototype)[invalid.cols]
    msg <- "Invalid columns in data.table:\n    %s"
    sprintf(msg, paste(bad.cols, collapse="\n    - "))
  } else {
    TRUE
  }
}
