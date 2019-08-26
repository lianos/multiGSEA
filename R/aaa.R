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
multiGSEA.methods <- function(names.only = TRUE) {
  methods <- list(
    camera = list(package = "edgeR", type = "required"),
    cameraPR = list(package = "edgeR", type = "required"),
    fgsea = list(package = "fgsea", type = "suggested"),
    geneSetTest = list(package = "edgeR", type = "required"),
    goseq = list(package = "goseq", type = "suggested"),
    logFC = list(package="multiGSEA", type = "required"),
    enrichtest = list(package = "limma", type = "required"),
    fry = list(package = "edgeR", type = "required"),
    roast = list(package = "edgeR", type = "required"),
    romer = list(package = "edgeR", type = "required"),
    svdGeneSetTest = list(package="multiGSEA", type = "required")
  )

  if (names.only) {
    methods <- names(methods)
  }

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

  mg.methods <- multiGSEA.methods(names.only = FALSE)
  bad.methods <- setdiff(methods, names(mg.methods))
  if (length(bad.methods)) {
    stop("unknown GSEA methods: ", paste(bad.methods, collapse=', '))
  }

  for (method in methods) {
    desc <- mg.methods[[method]]
    if (desc$type == "suggested") {
      if (!requireNamespace(method, quietly = TRUE)) {
        stop("The '", method, "' GSEA method requires the '", desc$package,
             "', which does not seem to be installed")
      }
    }
  }
}
