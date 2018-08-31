#' When you have an object with rownames, and you want to swap them by lookup.
#'
#' For instance, have a SummarizedExperiment, DGEList, matrix, etc. that have
#' ensembl id rowname`single cell experiment with ENSGID rownames that you
#' want to convert to symbols? Now you can.
#'
#' The rownames that can't successfully remapped will keep their old names.
#' This function should also guarantee that the rows of the incoming matrix
#' are the same as the outgoing one.
#'
#' @export
#'
#' @param x an object to whose rows need renaming
#' @param rename.rows an objec to help with the renaming.
#'   * A character vector where length(rename.rows) == nrow(x). Every row in
#'     x should correspond to the renamed value in the same position in
#'     rename.rows
#'   * If x is a DGEList, SummarizedExperiment, etc. this can be a string.
#'     In this case, this mus be a column in the fData/rowData/etc. The
#'     values in that column will be the new
#'   * A two column data.frame. The first column has entries in rownames(x),
#'     and the second column is the value to rename it to.
#'
rename_rows <- function(x, rename.rows, ...) {
  UseMethod("rename_rows", x)
}

.rename_rows.df <- function(x, rename.rows, rowmeta.df = NULL, ...) {
  stopifnot(is.character(x))
  if (is.data.frame(rowmeta.df)) {
    stopifnot(nrow(rowmeta.df) >= length(x))
  }

  if (is.character(rename.rows)) {
    if (length(rename.rows) == 1L && is.data.frame(rowmeta.df)) {
      rename.rows <- rowmeta.df[[rename.rows]]
    }
    if (length(rename.rows) != length(x)) {
      stop("rename.rows needs to be a character vector as long as nrow(x)")
    }
    rename.rows <- data.frame(oldnames = x, newnames = rename.rows,
                              stringsAsFactors = FALSE)
  }

  stopifnot(
    is.data.frame(rename.rows),
    ncol(rename.rows) == 2,
    is.character(rename.rows[[1]]), is.character(rename.rows[[2]]))
  # in case we got a tibble
  rename.rows <- as.data.frame(rename.rows, stringsAsFactors = FALSE)

  # If there is NA in rename_to column, use the value from first column
  rename.rows[[2]] <- ifelse(is.na(rename.rows[[2]]),
                             rename.rows[[1]], rename.rows[[2]])
  missed.x <- setdiff(x, rename.rows[[1]])
  if (length(missed.x)) {
    add.me <- data.frame(old = missed.x, new = missed.x)
    colnames(add.me) <- colnames(rename.rows)
    rename.rows <- rbind(rename.rows, add.me)
  }
  out <- rename.rows[!duplicated(rename.rows[[1]]),,drop = FALSE]
  out[[2]] <- ifelse(duplicated(out[[2]]), out[[1]], out[[2]])
  rownames(out) <- out[[1]]
  out[x,,drop=FALSE]
}

#' @export
rename_rows.default <- function(x, rename.rows = NULL, rowmeta.df = NULL, ...) {
  if (is.null(rename.rows)) return(x)
  if (length(dim(x)) != 2) stop("This thing isn't 2d-subsetable")
  rename.rows <- .rename_rows.df(rownames(x), rename.rows, rowmeta.df, ...)
  if (is.matrix(x)) {
    out <- x[rownames(rename.rows),,drop = FALSE]
  } else {
    out <- x[rownames(rename.rows),]
  }
  rownames(out) <- rename.rows[[2L]]
  out
}


#' @export
rename_rows.EList <- function(x, rename.rows = NULL, ...) {
  if (is.null(rename.rows)) return(x)
  # rename.rows <- .rename_rows.df(rownames(x), rename.rows, ...)
  #
  # out <- x[rownames(rename.rows),,drop = FALSE]
  # rownames(out) <- rename.rows$newnames
  # out
  out <- NextMethod(rowmeta.df = x$genes, ...)
  # rownames of genes isn't updated
  rownames(out$genes) <- rownames(out)
  out
}

#' @export
rename_rows.DGEList <- function(x, rename.rows = NULL, ...) {
  if (is.null(rename.rows)) return(x)
  # rename.rows <- .rename_rows.df(rownames(x), rename.rows, ...)
  #
  # out <- x[rownames(rename.rows),,drop = FALSE]
  # rownames(out) <- rename.rows$newnames
  # out
  out <- NextMethod(rowmeta.df = x$genes, ...)
  # rownames of genes isn't updated
  rownames(out$genes) <- rownames(out)
  out
}

#' @export
rename_rows.SummarizedExperiment <- function(x, rename.rows, ...) {
  if (is.null(rename.rows)) return(x)
  ns <- requireNamespace("SummarizedExperiment")
  # rename.rows <- .rename_rows.df(rownames(x), rename.rows, ...)
  # out <- x[rownames(rename.rows),,drop = FALSE]
  # rownames(out) <- rename.rows$newnames
  # out
  NextMethod(rowmeta.df = as.data.frame(rowData(x)), ...)
}

#' @export
rename_rows.eSet <- function(x, rename.rows, ...) {
  if (is.null(rename.rows)) return(x)
  ns <- requireNamespace("Biobaset")
  # rename.rows <- .rename_rows.df(rownames(x), rename.rows, ...)
  # out <- x[rownames(rename.rows),,drop = FALSE]
  # rownames(out) <- rename.rows$newnames
  # out
  NextMethod(rowmeta.df = Biobase::fData(x), ...)
}
