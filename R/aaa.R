## Helps internal data.table functions work when devoloping this package with
## devtools
.datatable.aware <- TRUE

## valid types of objects that can be used for "Expression" (x)'s
.valid.x <- c('matrix', 'eSet', 'EList', 'DGEList', 'SummarizedExperiment')

##' Lists the supported GSEA methods by multiGSEA
##'
##' @export
multiGSEA.methods <- function() {
  c('camera', 'roast', 'fry', 'romer', 'geneSetTest',
    'goseq', 'hyperGeometricTest', 'fgsea',
    'svdGeneSetTest', 'logFC')
}
.unsupportedGSEAmethods <- function(x, throw.error=TRUE) {
  bad.methods <- setdiff(x, multiGSEA.methods())
  if (length(bad.methods)) {
    stop("unknown GSEA methods: ", paste(bad.methods, collapse=', '))
  }
  bad.methods
}

