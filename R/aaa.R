# Helps internal data.table functions work when devoloping this package with
# devtools
.datatable.aware <- TRUE

# valid types of objects that can be used for "Expression" (x)'s
.valid.x <- c('matrix', 'eSet', 'EList', 'DGEList', 'SummarizedExperiment',
              'Matrix')

#' Lists the supported GSEA methods by multiGSEA
#'
#' @export
#' @param names.only if `TRUE` (default), only returns the names of the methods,
#'   otherwise also returns meta information about the methods, such as the
#'   package it is found in.
#' @return a character vector of GSEA names, or a list of metadata for each
#'   method.
multiGSEA_methods <- function() {
  methods <- list(
    camera = list(package = "edgeR", type = "required"),
    cameraPR = list(package = "edgeR", type = "required"),
    fgsea = list(package = "fgsea", type = "suggested"),
    enrichtest = list(package = "limma", type = "required"),
    fry = list(package = "edgeR", type = "required"),
    roast = list(package = "edgeR", type = "required"),
    romer = list(package = "edgeR", type = "required"),
    goseq = list(package = "goseq", type = "suggested"),
    geneSetTest = list(package = "edgeR", type = "required"),
    logFC = list(package="multiGSEA", type = "required"),
    svdGeneSetTest = list(package="multiGSEA", type = "required")
  )

  fn <- system.file("extdata", "gsea-methods.csv", package = "multiGSEA")
  methods <- read.csv(fn, stringsAsFactors = FALSE)
  methods
}

#' Helper function to check if a method is supported in multiGSEA
#'
#' @noRd
#' @param methods a character string of methods
check.gsea.methods <- function(methods) {
  if (!is.character(methods)) {
    stop("`methods` is not a character vector")
  }
  if (length(methods) == 0L) {
    stop("No `methods` are specified (length(methods) == 0)")
  }

  mg.methods <- multiGSEA_methods()
  bad.methods <- setdiff(methods, mg.methods[["method"]])
  if (length(bad.methods)) {
    stop("unknown GSEA methods: ", paste(bad.methods, collapse=', '))
  }

  for (method. in methods) {
    info <- subset(mg.methods, method == method.)
    if (info[["dependancy"]] == "suggested") {
      pkg <- info[["package"]]
      if (!requireNamespace(pkg, quietly = TRUE)) {
        stop("The '", method., "' GSEA method requires the '", pkg,
             "', which does not seem to be installed")
      }
    }
  }
}
