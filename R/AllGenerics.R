
##' Returns the featureIds associated with a given GeneSetDb or one of the
##' individual gene sets within it.
##'
##' @rdname featureIds
##' @exportMethod featureIds
setGeneric("featureIds", function(x, ...) standardGeneric("featureIds"))

##' Fetch the featureIdMap for a \code{GeneSetDb}
##' @exportMethod featureIdMap
setGeneric("featureIdMap", function(x, ...) standardGeneric("featureIdMap"))

setGeneric("featureIdMap<-", function(x, value) {
  standardGeneric("featureIdMap<-")
})


setGeneric("collectionMetadata", signature=c("x", "collection", "name"),
function(x, collection, name) {
  standardGeneric("collectionMetadata")
})

setGeneric("collectionMetadata<-", signature=c("x", "collection", "name"),
function(x, collection, name, value) {
  standardGeneric("collectionMetadata<-")
})

##' Provides an URL that provides a detailed description for a geneset
##'
##' This function requires that the \code{GeneSetDb@@collectionMetadata} slot has
##' a valid \code{url.fn} column that provides a function to parse the
##' collection,name information for its genesets into a valid URL. The functions
##' in \code{url.fn} should take two arguments which correspond to the geneset
##' \code{collection} and \code{name}, respectively.
##'
##' @exportMethod geneSetURL
setGeneric("geneSetURL", signature="x", function(x, i, j, ...) {
  standardGeneric("geneSetURL")
})

##' exportMethod geneSetCollectionURLfunction
setGeneric("geneSetCollectionURLfunction", signature="x", function(x, i, ...) {
  standardGeneric("geneSetCollectionURLfunction")
})

setGeneric("geneSetCollectionURLfunction<-", signature="x", function(x, i, value) {
  standardGeneric("geneSetCollectionURLfunction<-")
})

##' @exportMethod geneSet
setGeneric("geneSet", signature="x", function(x, i, j, ...) {
  standardGeneric("geneSet")
})

##' Fetch the gene sets "active" (or otherwise) in the GeneSetDb or Result
##'
##' @export
##' @param x GeneSetDb
##' @return a data.table with geneset information.
##'
##' @rdname geneSets
##' @exportMethod geneSets
setGeneric("geneSets", function(x, ...) standardGeneric('geneSets'))

##' Map featureId's in a GeneSetDb to a target expression object.
##'
##' @exportMethod conform
setGeneric("conform", function(x, ...) standardGeneric("conform"))

##' @exportMethod unconform
setGeneric("unconform", function(x, ...) standardGeneric("unconform"))

##' Summarizes different results into tabular form
##'
##' @exportMethod summarized
setGeneric("summarized", function(x, ...) standardGeneric("summarized"))

##' @exportMethod plot
setGeneric("plot")
