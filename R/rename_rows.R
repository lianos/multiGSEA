#' Smartly/easily rename the rows of an object.
#'
#' The most common usecase for this is when you have a SummarizedExperiment,
#' DGEList, matrix, etc. that is "rownamed" by some gene idnetifiers (ensembl,
#' entrez, etc) that you want to "easily" convert to be rownamed by symbols.
#' And perhaps the most common use-case for this, again, would be able to
#' easily change rownames of a heatmap to symbols.
#'
#' The rownames that can't successfully remapped will keep their old names.
#' This function should also guarantee that the rows of the incoming matrix
#' are the same as the outgoing one.
#'
#' @export
#'
#' @param x an object to whose rows need renaming
#' @param xref an object to help with the renaming.
#'   * A character vector where length(xref) == nrow(x). Every row in
#'     x should correspond to the renamed value in the same position in
#'     xref
#'   * If x is a DGEList, SummarizedExperiment, etc. this can be a string.
#'     In this case, the string must name a column in the data container's
#'     fData-like data.frame. The values in that column will be the new
#'     candidate rownames for the object.
#'   * A two column data.frame. The first column has entries in rownames(x),
#'     and the second column is the value to rename it to.
#' @examples
#' eset <- exampleExpressionSet(do.voom = FALSE)
#' ess <- rename_rows(eset, "symbol")
#'
#' vm <- exampleExpressionSet(do.voom = TRUE)
#' vms <- rename_rows(vm, "symbol")
rename_rows <- function(x, xref, ...) {
  UseMethod("rename_rows", x)
}

#' Returns a two-column data.frame with rownames. The rownames are entriex from
#' x, first should be the same, and the second column is the value that x
#' should be renamed to.
#'
#' @noRd
.rename_rows.df <- function(x, xref = NULL, rowmeta.df = NULL, ...) {
  stopifnot(is.character(x))
  if (!is.data.frame(xref)) {
    stopifnot(
      is.character(xref),
      length(xref) == length(x))
    xref <- data.frame(from = x, to = xref, stringsAsFactors = FALSE)
  }
  if (is(xref, "tbl") || is(xref, "data.table")) {
    xref <- as.data.frame(xref, stringsAsFactors = FALSE)
  }
  stopifnot(
    is.data.frame(xref),
    ncol(xref) == 2L,
    is.character(xref[[1L]]), is.character(xref[[2L]]))

  # If there is NA in rename_to column, use the value from first column
  xref[[2L]] <- ifelse(is.na(xref[[2L]]), xref[[1L]], xref[[2L]])

  # Are there entries in x that don't appear in first colum of xref? If so,
  # we expand `xref` to include these entries and have them "remap" to identity
  missed.x <- setdiff(x, xref[[1L]])
  if (length(missed.x)) {
    add.me <- data.frame(old = missed.x, new = missed.x)
    colnames(add.me) <- colnames(xref)
    xref <- rbind(xref, add.me)
  }

  # Remove ambiguity in remapping process. If the same original ID can be
  # remapped to several other ones, then only one will be picked.
  out <- xref[!duplicated(xref[[1L]]),,drop = FALSE]
  # If there are duplicated values in the entries that x can be translated to,
  # then those renamed entries will remap x to itself
  out[[2L]] <- ifelse(duplicated(out[[2L]]), out[[1L]], out[[2L]])
  out
}

#' @export
rename_rows.default <- function(x, xref = NULL, ...) {
  if (is.null(xref)) {
    warning("No `xref` provided, returning object unchanged", immediate. = TRUE)
    return(x)
  }
  rn <- rownames(x)
  if (is.null(rn)) {
    warning("`x` has no rownames, returning as is", immediate. = TRUE)
    return(x)
  }
  if (length(dim(x)) != 2L) {
    stop("The input object isn't 2d-subsetable")
  }
  xref <- .rename_rows.df(rownames(x), xref, ...)
  nomatch <- setdiff(rn, xref[[1L]])
  if (length(nomatch)) {
    stop(length(nomatch), " rownames do not have a lookup to use in renaming")
  }
  lookup <- match(rn, xref[[1L]])
  rownames(x) <- xref[[2L]][lookup]
  x
}

#' @noRd
.rename_rows.bioc <- function(x, xref = NULL, ...) {
  if (is.null(xref)) return(x)
  if (is.character(xref) && length(xref) == 1L) {
    xref <- data.frame(from = rownames(x),
                       to = fdata(x)[[xref]],
                       stringsAsFactors = FALSE)
  }
  out <- rename_rows.default(x, xref = xref, ...)
  rownames(fdata(out)) <- rownames(out)
  out
}


#' @export
#' @noRd
rename_rows.EList <- function(x, xref = NULL, ...) {
  .rename_rows.bioc(x, xref, ...)
}

#' @export
#' @noRd
rename_rows.DGEList <- function(x, xref = NULL, ...) {
  .rename_rows.bioc(x, xref, ...)
}

#' @export
#' @noRd
rename_rows.SummarizedExperiment <- function(x, xref, ...) {
  .rename_rows.bioc(x, xref, ...)
}

#' @export
#' @noRd
rename_rows.eSet <- function(x, xref, ...) {
  .rename_rows.bioc(x, xref, ...)
}
