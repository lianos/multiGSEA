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
