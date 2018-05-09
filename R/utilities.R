extract_preranked_stats <- function(x, design, contrast, robust.fit=FALSE,
                                    robust.eBayes=FALSE, logFC=NULL,
                                    score.by=c('t', 'logFC', 'pval'), ...) {
  score.by <- match.arg(score.by)
  if (ncol(x) > 1) {
    if (is.null(logFC)) {
      logFC <- calculateIndividualLogFC(x, design, contrast, robust.fit,
                                        robust.eBayes, ..., as.dt=TRUE)
    } else {
      is.logFC.like(logFC, x, as.error=TRUE)
    }
    ## t will be NA if statistics were computed using edgeR from a DGEList
    if (score.by == 't' && any(is.na(logFC[['t']]))) {
      warning("t statistics not found in dge results, switching to logFC",
              immediate.=TRUE)
      score.by <- 'logFC'
    }
    stats <- setNames(logFC[[score.by]], logFC$featureId)
  } else {
    ## This is already a column matrix of precomputed things (logFC, perhaps)
    ## to rank
    stats <- setNames(x[, 1L], rownames(x))
  }

  if (any(is.na(stats))) {
    stop("NA values are found in the stats vector for geneSetTest")
  }

  if (!setequal(rownames(x), names(stats))) {
    stop("Identifiers are not setequal among stats vector and x matrix")
  }

  stats[rownames(x)]
}

##' Converts collection,name combination to key for geneset
##'
##' The "key" form often comes out as rownames to matrices and such, or
##' particularly for sending genesets down into wrapped methods, like do.camera.
##'
##' @export
##' @rdname gskey
##' @return a character vector
encode_gskey <- function(x, y, sep=";;") {
  if (is.data.frame(x)) {
    stopifnot(
      "collection" %in% colnames(x),
      "name" %in% colnames(x))
    y <- x[['name']]
    x <- x[['collection']]
  }
  if (is.factor(x)) x <- as.character(x)
  if (is.factor(y)) y <- as.character(y)
  stopifnot(is.character(x), is.character(y))
  paste(x, y, sep=sep)
}

##' Splits collection,name combinations to collection,name data.frames
##'
##' @export
##' @rdname gskey
##' @return a data.frame with (collection,name) columns
split_gskey <- function(x, sep=";;") {
  stopifnot(all(grepl(sep, x)))
  data.frame(
    collection=sub(sprintf("%s.*$", sep), "", x),
    name=sub(sprintf(".*?%s", sep), "", x),
    stringsAsFactors=FALSE)
}

## Helper function to extract a vector of "pre-ranked" stats for a GSEA. This
## can come from: (1) a user provided vector; (2) the logFCs or t-stats of
## an internal call to calculalateIndividualLogFCs from a "full design"ed
## matrix
generate.preranked.stats <- function(x, design, contrast, logFC=NULL,
                                     score.by=c('t', 'logFC', 'pval')) {
  if (!is.null(logFC)) {
    is.logFC.like(logFC, x, as.error=TRUE)
    score.by <- match.arg(score.by)
    out <- setNames(logFC[[score.by]], logFC[['featureId']])
  } else {
    ## If multiGSEA was called with a preranked vector, the validateInputs function
    ## would have converted it into a column matrix with rownames, but most
    ## preranked functions want a named vector
    out <- setNames(as.vector(x), rownames(x))
  }
  out
}


##' Converts an expression container like object to a matrix for analysis
##'
##' This function is not exported. It will convert various expression-like
##' containers into a a matrix of values for use in GSEA or single-sample based
##' scoring methods.
##'
##' There's nothing too fancy here. Keep in mind, however, that if \code{y}
##' is something that typically stores counts (a \code{DGEList} or
##' \code{DESeqDataSet}), then it is transformed into a matrix of values
##' on the log scale.
##'
##' @param y an object to convert into an expression matrix for use in various
##'   internal gene set based methods
##' @return a matrix of values to use downstream of internal gene set based
##'   methods.
as_matrix <- function(y) {
  if (is.vector(y)) {
    y <- t(t(y)) ## column vectorization that sets names to rownames
  } else if (is(y, 'EList')) {
    y <- y$E
  } else if (is(y, 'DGEList')) {
    y <- cpm(y, prior.count=5, log=TRUE)
  } else if (is(y, 'eSet')) {
    y <- Biobase::exprs(y)
  } else if (is(y, 'DESeqDataSet')) {
    y <- SummarizedExperiment::assay(DESeq2::normTransform(y, pc=5))
  } else if (is.data.frame(y)) {
    y <- as.matrix(y)
  } else if (is(y, "SingleCellExperiment")) {
    if (!"logcounts" %in% assayNames(y)) {
      stop("`logcounts` assay missing from SingleCellExperiment")
    }
    y <- SummarizedExperiment::assay(y, "logcounts")
  }
  stopifnot(is.matrix(y) && is.numeric(y))
  y
}

##' Reads in a semi-annotated genelist (one symbol per line)
##'
##' Often we are given a list of gene names, and the symbols provided are not
##' the official HGNC ones. In these cases (when small enough) I will replace
##' the symbol provided by the official one, and leave the submitted symbol
##' there after a comment character ("#"). This just strips the stuff after
##' the comment character to provide only the offiical symbols.
##'
##' @export
##' @param fn the path to the gene list file
##' @return character vector of gene names.
readGeneSymbols <- function(fn) {
  out <- readLines(fn)
  sub(' +#.*', '', out)
}

isSingleCharacter <- function(x, allow.na=FALSE) {
  is.character(x) && length(x) == 1L && (!is.na(x) || allow.na)
}

isSingleInteger <- function(x, allow.na=FALSE) {
  is.integer(x) && length(x) == 1L && (!is.na(x) || allow.na)
}

isSingleNumeric <- function(x, allow.na=FALSE) {
  is.numeric(x) && length(x) == 1L && (!is.na(x) || allow.na)
}

isSingleLogical <- function(x, allow.na=FALSE) {
  is.logical(x) && length(x) == 1L && (!is.na(x) || allow.na)
}

##' Check a data.table vs a reference data.table
##'
##' This function ensures that a \code{data.table} \code{x} has a superset of
##' columns of a reference table \code{ref}, and that both tables are keyed by
##' the same columns.
##'
##' @param x A \code{data.table} you want to be validated
##' @param ref The \code{data.table} to use as the reference/model data.table
##'
##' @return \code{TRUE} if all things check out, otherwise a character vector
##'   indicating what the problems were.
check.dt <- function(x, ref) {
  if (!is.data.table(x)) {
    stop("Input is not a data.table")
  }
  missed.cols <- setdiff(names(ref), names(x))
  if (length(missed.cols)) {
    msg <- paste('columns missing:', paste(missed.cols, collapse=','))
    return(msg)
  }

  pk <- key(ref)
  if (!is.null(pk)) {
    xk <- key(x)
    if (length(pk) != length(xk) || !all(pk == xk)) {
      return('illegal key set')
    }
  }

  TRUE
}

## -----------------------------------------------------------------------------
## File Utilities

##' Utility function to determine if paths listed are real directories
##'
##' @param ... character vectors containing file paths. Tilde-expansion is
##' done: see 'path.expand'.
##'
##' @return A logical vector indicating whether or not the queried directories
##' exist.
dir.exists <- function(...) {
  sapply(file.info(...)$isdir, isTRUE)
}

##' Checks that the directory provided is writable by the current user.
##'
##' This works by testing to put a temporary file into an already existing
##' directory
##'
##' @param path The path to a directory to check.
##'
##' @return \code{logical}, \code{TRUE} if \code{path} is writable by the
##' current user, otherise \code{FALSE}
dir.writable <- function(path) {
  if (!dir.exists(path)) {
    stop("The directory provided does not exist: ", path)
  }
  tmp.fn <- tempfile(tmpdir=path, fileext='.test.tmp')
  on.exit({
    if (file.exists(tmp.fn)) unlink(tmp.fn)
  })

  tryCatch({
    suppressWarnings(writeLines('test', tmp.fn))
    TRUE
  }, error=function(e) FALSE)
}

## Random Utilitis -------------------------------------------------------------
##' Utility function to cat a message to stderr (by default)
##'
##' @export
##' @param ... pieces of the message
##' @param file where to send the message. Defaults to \code{stderr()}
msg <- function(..., file=stderr()) {
  cat(paste(rep('-', 80), collapse=''), '\n', file=file)
  cat(..., '\n', file=file)
  cat(paste(rep('-', 80), collapse=''), '\n', file=file)
}

##' Utility function to try and fail with grace.
##'
##' @export
##' @param default the value to return if \code{expr} fails
##' @param expr the expression to take a shot at
##' @param frame the frame to evaluate the expression in
##' @param message the error message to display if \code{expr} fails. Deafults
##'   to \code{geterrmessage()}
##' @param silent if \code{TRUE}, sends the error message to \code{msg}
##' @param file where msg sends the message
##' @return the result of \code{expr} if successful, otherwis \code{default}
failWith <- function(default=NULL, expr, frame=parent.frame(),
                     message=geterrmessage(), silent=FALSE, file=stderr()) {
  tryCatch(eval(expr, frame), error=function(e) {
    if (!silent) msg(message, file)
    default
  })
}
