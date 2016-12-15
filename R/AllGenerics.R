##' Returns the relevant featureIds for a given geneset.
##'
##' @description
##' Gene sets are defined by the unique compound key consisting of their
##' \code{collection} and \code{name}. To fetch the featureIds associated with
##' a specific geneset, you must provide values for \code{i} and \code{j}. If
##' these are missing, then a character vector of all the unique feature ids
##' within \code{x} are returned.
##'
##' If the GeneSetDb \code{x} has been conformed to an expression object this
##' will default to return the featureId's as they are used/matched to the
##' expression object, otherwise it will return the featureIds used in the
##' definition of the gene set database.
##'
##' \code{x} can be either a \code{GeneSetDb} or a \code{MultiGSEAResult}. If
##' its the latter, then this call simply delegates to the internal
##' \code{GeneSetDb}.
##'
##' @rdname featureIds
##' @exportMethod featureIds
##'
##' @return A vector of identifiers (or indexes into an expression object,
##'   depending on the \code{value} argument) for the features in the specified
##'   geneset. NA is returned if the geneset is not "active" (ie. listed in
##'   geneSets(x))
setGeneric("featureIds", function(x, ...) standardGeneric("featureIds"))

##' Fetch the featureIdMap for a \code{GeneSetDb}
##' @exportMethod featureIdMap
setGeneric("featureIdMap", function(x, ...) standardGeneric("featureIdMap"))

setGeneric("featureIdMap<-", function(x, value) {
  standardGeneric("featureIdMap<-")
})


##' Query for metadata associated with a collection.
##'
##' @exportMethod collectionMetadata
##' @param x Object to extract the collectionMetadata from
##' @param collection The geneset collection to to query
##' @param name The name of the metadata variable to get the value for
setGeneric("collectionMetadata", signature=c("x", "collection", "name"),
function(x, collection, name, ...) {
  standardGeneric("collectionMetadata")
})

##' @export
setGeneric("collectionMetadata<-", signature=c("x", "collection", "name"),
function(x, collection, name, value) {
  standardGeneric("collectionMetadata<-")
})

##' Constructs URL that provides a detailed description for a geneset
##'
##' This function requires that the \code{GeneSetDb@@collectionMetadata} slot has
##' a valid \code{url.fn} column that provides a function to parse the
##' collection,name information for its genesets into a valid URL. The functions
##' in \code{url.fn} should take two arguments which correspond to the geneset
##' \code{collection} and \code{name}, respectively.
##'
##' This function is vectorized over \code{i} and \code{j}
##'
##' @exportMethod geneSetURL
##'
##' @param x A \code{GeneSetDb} or \code{MultiGSEAResult}
##' @param i The names of the collections
##' @param j The names of the gene sets belonging to the corresponding
##'   collection identified in \code{i}
##'
##' @return A character vector of URLs for each of the genesets identified by
##'   \code{i, j}. \code{NA} is returned for genesets \code{i,j} that are not
##'   found in \code{x}.
setGeneric("geneSetURL", signature="x", function(x, i, j, ...) {
  standardGeneric("geneSetURL")
})

##' exportMethod geneSetCollectionURLfunction
setGeneric("geneSetCollectionURLfunction", signature="x", function(x, i, ...) {
  standardGeneric("geneSetCollectionURLfunction")
})

##' @export
setGeneric("geneSetCollectionURLfunction<-", signature="x", function(x, i, value) {
  standardGeneric("geneSetCollectionURLfunction<-")
})

##' @export
setGeneric("collectionUrlFunction<-", signature="x", function(x, i, value, ...) {
  standardGeneric("collectionUrlFunction<-")
})

##' Query/set the featureIdType for an collection
##'
##' EntrezIdIdentifer() or?
##' @exportMethod featureIdType
setGeneric("featureIdType", signature="x", function(x, i, ...) {
  standardGeneric("featureIdType")
})

##' @export
setGeneric("featureIdType<-", signature="x", function(x, i, value) {
  standardGeneric("featureIdType<-")
})

##' Query/set the organism for a given geneset (db or collection)
##' @exportMethod org
##' @rdname org
setGeneric("org", signature="x", function(x, i, ...) {
  standardGeneric("org")
})

##' @export
##' @rdname org
setGeneric("org<-", signature="x", function(x, i, value) {
  standardGeneric("org<-")
})


##' Fetches information about the genes in a geneset.
##'
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

##' Creates a geneset summary table for genesets defined by gene membership
##'
##' @exportMethod geneSetSummaryByGenes
##' @param x A GeneSetDb or a multiGSEA ResultmultiGSEA result
##' @param features The featureIds to use to select genesets from \code{x}
##' @param with.features Inlcude the features as columns of the outgoing result?
##'   (Default: \code{TRUE})
##' @param feature.rename If features are added as columns, then we might want
##'   to rename them from featureIds to something else. By default these columns
##'   are prefixed with "featureId_", unless \code{feature.rename=FALSE}
##' @param name If specified, one of \code{\link{resultNames}} from the
##'   \code{MultiGSEAResult} object \code{x} to use as a look up and filter
##'   geneset by.
##' @param max.p The FDR threshold to filter geneset results by in combination
##'   with the \code{name} parameter.
##' @return A summary data.frame
setGeneric("geneSetSummaryByGenes", signature=c("x"),
           function(x, features, with.features=TRUE, feature.rename=NULL, ...)
             standardGeneric("geneSetSummaryByGenes"))


##' Map featureId's in a GeneSetDb to a target expression object.
##'
##' @exportMethod conform
##' @param x \code{GeneSetDb}
setGeneric("conform", function(x, ...) standardGeneric("conform"))

##' Subset a GeneSetDb to only include geneSets with specified features.
##'
##' @param x \code{GeneSetDb}
##' @param featureIds Character vector of featureIds
##' @return A subset of \code{x} which contains only the geneSets that contain
##'   features found in \code{featureIds}
##' @exportMethod subsetByFeatures
setGeneric("subsetByFeatures", function(x, features, ...) {
  standardGeneric("subsetByFeatures")
})

##' @exportMethod unconform
setGeneric("unconform", function(x, ...) standardGeneric("unconform"))

##' Summarizes different results into tabular form
##'
##' @exportMethod summarized
setGeneric("summarized", function(x, ...) standardGeneric("summarized"))

##' @exportMethod plot
setGeneric("plot")
