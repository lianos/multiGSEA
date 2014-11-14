##' @rdname featureIds
##' @export
setGeneric("featureIds", function(object, ...) {
  standardGeneric("featureIds")
})

##' Conforms (GeneSetTable) x to (ExpressionSet) y, for instance.
##'
##' @exportMethod conform
setGeneric("conform", function(x, ...) standardGeneric("conform"))

