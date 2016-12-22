## Helps internal data.table functions work when devoloping this package with
## devtools
.datatable.aware <- TRUE

## .multi.gsea.methods <- c('camera', 'roast', 'gsd', 'npGSEA')
.multi.gsea.methods <- c('camera', 'roast', 'fry', 'romer', 'geneSetTest',
                         'gseap', 'goseq', 'hyperGeometricTest',
                         'svdGeneSetTest', 'logFC')

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
##' The data.table package is heavily used internally in this package, and I
##' want those objects to (always) be data.table objects while they are passed
##' around from one function to the next. When data.tables are passed back to
##' the user, they will be automatically converted into to data.frames, by
##' default.
##'
##' Internal multiGSEA functions that return a data.table are always
##' (1) parameterized with an \code{.external=TRUE} parameter; and
##' (2) pass the outgoing data.table through this function.
##'
##' If \code{.external=TRUE}, then the data.table will be transformed on the way
##' out back to whatever is set for \code{getOption('multiGSEA.df.return')}.
##'
##' \code{getOption('multiGSEA.df.return')} is automatically set to "data.frame"
##' when the package is loaded if this option isn't already set in the user's
##' enviroment (ie. in their \code{.Rprofile}, for instance). The assumption is
##' that user's are most accustomed to using data.frames, and some of the
##' idiosyncracies of indexing into a data.table won't trip them up.
##'
##' Does this extra monkey business make the internal code of multiGSEA a bit
##' more painful to write? Yes. So ... you're welcome.
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

  x <- copy(x)
  clazz <- class(x)[1L]
  switch(df.return,
         data.table=if (is.data.table(x)) x else setDT(x),
         data.frame=if (clazz == 'data.frame') {
           x
         } else if (is.data.table(x)) {
           setDF(x)
          } else {
            warning("Unknown object type, returning as is: ", clazz)
            x
          })
}
