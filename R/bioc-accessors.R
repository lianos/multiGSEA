#' Wrapper for easy fetching of rowData/fData, etc.
#'
#' @noRd
fdata <- function(x, as.df = FALSE, ...) {
  UseMethod("fdata", x)
}

`fdata<-` <- function(x, value) {
  UseMethod("fdata<-", x)
}

#' @noRd
pdata <- function(x, as.df = FALSE, ...) {
  UseMethod("pdata", x)
}

`pdata<-` <- function(x, value) {
  UseMethod("pdata<-", x)
}

# DGEList ======================================================================
fdata.DGEList <- function(x, as.df = FALSE, ...) {
  out <- x$genes
  if (!is.data.frame(out)) stop("No `genes` data.frame found in (DG)EList")
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
  if (!requireNamespace("Biobase")) stop("Biobase package required")
  Biobase::fData(x)
}

`fdata<-.eSet` <- function(x, value) {
  if (!requireNamespace("Biobase")) stop("Biobase package required")
  # `Biobase::fData<-`(x, value)
  Biobase::fData(x) <- value
  x
}

pdata.eSet <- function(x, ...) {
  if (!requireNamespace("Biobase")) stop("Biobase package required")
  Biobase::pData(x)
}

`pdata<-.eSet` <- function(x, value) {
  if (!requireNamespace("Biobase")) stop("Biobase package required")
  Biobase::pData(x) <- value
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
