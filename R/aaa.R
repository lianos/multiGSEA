## Helps internal data.table functions work when devoloping this package with
## devtools
.datatable.aware <- TRUE

.multi.gsea.methods <- c('camera', 'roast', 'gst', 'npGSEA')

.unsupportedGSEAmethods <- function(x, throw.error=TRUE) {
  bad.methods <- setdiff(x, .multi.gsea.methods)
  if (length(bad.methods)) {
    stop("unknown GSEA methods: ", paste(bad.methods, collapse=', '))
  }
  bad.methods
}
