
##' Returns the featureIds associated with a given GeneSetTable or one of the
##' individual gene sets within it.
##'
##' @rdname featureIds
##' @exportMethod
setGeneric("featureIds", function(object, ...) {
  standardGeneric("featureIds")
})

##' Conforms (GeneSetTable) x to (ExpressionSet) y, for instance.
##'
##' @exportMethod conform
setGeneric("conform", function(x, ...) standardGeneric("conform"))

