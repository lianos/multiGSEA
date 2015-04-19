## Helps internal data.table functions work when devoloping this package with
## devtools
.datatable.aware <- TRUE

## .multi.gsea.methods <- c('camera', 'roast', 'gsd', 'npGSEA')
.multi.gsea.methods <- c('camera', 'roast', 'geneSetTest', 'hyperGeometricTest',
                         'logFC')

## valid types of objects that can be used for "Expression" (x)'s
.valid.x <- c('matrix', 'eSet', 'EList', 'DGEList', 'SummarizedExperiment')

.unsupportedGSEAmethods <- function(x, throw.error=TRUE) {
  bad.methods <- setdiff(x, .multi.gsea.methods)
  if (length(bad.methods)) {
    stop("unknown GSEA methods: ", paste(bad.methods, collapse=', '))
  }
  bad.methods
}
