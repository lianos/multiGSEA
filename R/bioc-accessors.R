#' Wrapper for easy fetching of rowData/fData, etc.
#'
#' @noRd
fdata <- function(x, as.df = FALSE, ...) {
  UseMethod("fdata", x)
}

#' @noRd
#' @export
fdata.default <- function(x, ...) NULL

`fdata<-` <- function(x, value) {
  UseMethod("fdata<-", x)
}

#' @noRd
pdata <- function(x, as.df = FALSE, ...) {
  UseMethod("pdata", x)
}

#' @noRd
#' @export
pdata.default <- function(x, ...) NULL

`pdata<-` <- function(x, value) {
  UseMethod("pdata<-", x)
}

# DGEList ======================================================================
fdata.DGEList <- function(x, as.df = FALSE, ...) {
  out <- x$genes
  if (!is.data.frame(out)) {
    warning("No `genes` data.frame found in (DG)EList", immediate. = TRUE)
  }
  out
}

`fdata<-.DGEList` <- function(x, value) {
  x$genes <- value
  x
}

pdata.DGEList <- function(x, ...) {
  x$samples
}

`pdata<-.DGEList` <- function(x, value) {
  x$samples <- value
  x
}

# EList ========================================================================
fdata.EList <- fdata.DGEList
`fdata<-.EList` <- `fdata<-.DGEList`

pdata.EList <- function(x, ...) {
  x$targets
}

`pdata<-.EList` <- function(x, value) {
  x$targets <- value
  x
}

# eSet =========================================================================
fdata.eSet <- function(x, ...) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  ns$fData(x)
}

`fdata<-.eSet` <- function(x, value) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  x <- ns$`fData<-`(x, value)
  x
}

pdata.eSet <- function(x, ...) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  ns$pData(x)
}

`pdata<-.eSet` <- function(x, value) {
  ns <- tryCatch(loadNamespace("Biobase"), error = function(e) NULL)
  if (is.null(ns)) stop("Biobase required")
  x <- ns$`pData<-`(x, value)
  x
}

# SummarizedExperiment =========================================================

fdata.SummarizedExperiment <- function(x, as.df = FALSE, ...) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("SummarizedExperiment package required")
  }
  out <- SummarizedExperiment::rowData(x)
  if (as.df) out <- SummarizedExperiment::as.data.frame(out)
  out
}

`fdata<-.SummarizedExperiment` <- function(x, value) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("SummarizedExperiment package required")
  }
  SummarizedExperiment::rowData(x) <- value
  x
}

pdata.SummarizedExperiment <- function(x, as.df = FALSE, ...) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("SummarizedExperiment package required")
  }
  out <- SummarizedExperiment::colData(x)
  if (as.df) out <- SummarizedExperiment::as.data.frame(out)
  out
}

`pdata<-.SummarizedExperiment` <- function(x, value) {
  if (!requireNamespace("SummarizedExperiment")) {
    stop("SummarizedExperiment package required")
  }
  SummarizedExperiment::colData(x) <- value
  x
}
