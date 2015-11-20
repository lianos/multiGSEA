## Helps internal data.table functions work when devoloping this package with
## devtools
.datatable.aware <- TRUE

## .multi.gsea.methods <- c('camera', 'roast', 'gsd', 'npGSEA')
.multi.gsea.methods <- c('camera', 'roast', 'fry', 'romer', 'geneSetTest',
                         'goseq', 'hyperGeometricTest', 'logFC')

## valid types of objects that can be used for "Expression" (x)'s
.valid.x <- c('matrix', 'eSet', 'EList', 'DGEList', 'SummarizedExperiment')

.unsupportedGSEAmethods <- function(x, throw.error=TRUE) {
  bad.methods <- setdiff(x, .multi.gsea.methods)
  if (length(bad.methods)) {
    stop("unknown GSEA methods: ", paste(bad.methods, collapse=', '))
  }
  bad.methods
}

##' Returns data.table objects in formats people like.
##'
##' Apparently some people don't like getting data.table(s). Use
##' \code{options(multiGSEA.df.return='data.frame')} if you want to get
##' a \code{data.frame}
ret.df <- function(x, .external=TRUE) {
  if (!.external) {
    return(x)
  }
  if (!is(x, 'data.frame')) {
    return(x)
  }
  df.return <- getOption('multiGSEA.df.return')
  if (!df.return %in% c('data.table', 'data.frame')) {
    warning("Invalid value for options(multiGSEA.df.return) (",
            df.return, "), returning default object", immediate.=TRUE)
    return(x)
  }

  clazz <- class(x)[1L]
  switch(df.return,
         data.table=if (is.data.table(x)) x else setDT(copy(x)),
         data.frame=if (clazz == 'data.frame') {
           x
         } else if (is.data.table(x)) {
           setDF(copy(x))
          } else {
            warning("Unknown object type, returning as is: ", clazz)
            x
          })
}
