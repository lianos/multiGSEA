##' @exportClass
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
##' @param xref A named character vector which maps the rownames of \code{x}
##' to the gene ids in the \code{gene.sets}. The values of this vector are
##' the IDs used in \code{gene.sets} (typically entrez.id) and the \code{names}
##' are the \code{rownames} of \code{x} (the features in x)
GeneSetTable <- function(x, gene.sets, xref=NULL, min.gs.size=5,
                         stop.on.duplicate.x.xref=TRUE,
                         unique.by=c('mean', 'var'),
                         species=.wehi.msigdb.species,
                         version=.wehi.msigdb.current) {
  species <- match.arg(species)
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

  if (is.character(xref)) {
    if (any(duplicated(names(xref)))) {
      stop("Duplicated names in `xref`: this is a no go!")
    }
    if (any(duplicated(xref))) {
      if (stop.on.duplicate.x.xref) {
        stop("Duplicate geneset ID matches to rows in x")
      }
      ## TODO: Need to take multiple/duplicated ids in \code{x} into account,
      ##       for instance when `x` are rows from a microarray experiment, and
      ##       each row is a separate probeset -- many rows will be measuring
      ##       the same gene. In this case you probably want to keep the row
      ##       that is unique in some way, ie: highest average expression or
      ##       variability.
      warning("Duplicate matches from geneset IDs to rownames(x) not fully ",
              "implemented yet",
              immediate.=TRUE)
    }
    lookup <- data.table(gset.id=xref,
                         x.id=names(x.id),
                         x.index=match(names(x.id), rownames(x)),
                         key='gset.id')
    lookup <- lookup[is.na(x.index) == FALSE]
  } else {
    lookup <- data.table(gset.id=rownames(x),
                         x.id=rownames(x),
                         x.index=1:nrow(x),
                         key='gset.id')
  }

  if (is(gene.sets, 'GeneSetTable')) {
    out <- gene.sets
    validObject(out)
  } else {
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
  }

  if (min.gs.size > 0) {
    out@table <- out@table[n >= min.gs.size]
  }

  out
}

##' Conforms x to y.
setGeneric("conform", function(x, y, ...), standardGeneric("conform"))

setMethod("conform", c(x='GeneSetTable', y='ANY'),
function(x, y, x.id.fn=rownames, lookup=x@feature.lookup, ...) {
  ## Will update x@feature.lookup to match to the rows in x and return a
  ## GeneSetTable that's ready to be analyzed against `data`
  x.id <- x.id.fn(x)
  ## TODO: Finish up `conform,(GeneSetTable,ANY)`
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

##' @param x An "expression thing"
##' @param gene.sets A list of lists of genesets into \code{x}
##' @param xref A named character vector which maps the rownames of \code{x}
##' to the gene ids in the \code{gene.sets}. The values of this vector are
##' the IDs used in \code{gene.sets} (typically entrez.id) and the \code{names}
##' are the \code{rownames} of \code{x} (the features in x)
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

