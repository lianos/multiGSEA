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
##' will default to return only the featureId's that are matched to the target
##' expression object, and they will be returned using the same identifiers that
##' the target expression object uses. To change this behavior, tweak the values
##' for the \code{activee.only} and \code{value} parameters, respectively.
##'
##' \code{x} can be either a \code{GeneSetDb} or a \code{MultiGSEAResult}. If
##' its the latter, then this call simply delegates to the internal
##' \code{GeneSetDb}.
##'
##' @rdname featureIds
##' @exportMethod featureIds
##'
##' @param x Object to retrieve the gene set from, either a \code{GeneSetDb} or
##'   a \code{MultiGSEAResult}.
##' @param i,j The collection,name compound key identifier of the gene set
##' @param value What form do you want the id's in?
##'   \describe{
##'     \item{featureId}{the IDs used in the original geneset definitions}
##'     \item{x.id}{the ids of the features as they are used in the expression
##'           object}
##'     \item{x.idx}{The integer index into the expresion object \code{x} that
##'           the GeneSetDb has been conformed to.}
##'   }
##' @param active.only only look for gene sets that are "active"? Defaults to
##'   \code{TRUE} if \code{x} is conformed to a target expression object, else
##'   \code{FALSE}. See \code{\link{conform}} for further details.
##' @return A vector of identifiers (or indexes into an expression object,
##'   depending on the \code{value} argument) for the features in the specified
##'   geneset. NA is returned if the geneset is not "active" (ie. listed in
##'   geneSets(x))
##'
##' @examples
##' gdb <- exampleGeneSetDb()
##' fids.gs <- featureIds(gdb, 'c2', 'BIOCARTA_AGPCR_PATHWAY')
##' fids.c2 <- featureIds(gdb, 'c2')
##' fids.all <- featureIds(gdb)
##'
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- conform(gdb, vm)
##' ## fewer than before
##' fids.gs2 <- featureIds(gdb, 'c2', 'BIOCARTA_AGPCR_PATHWAY')
##' ## same as before
##' fids.gs3 <- featureIds(gdb, 'c2', 'BIOCARTA_AGPCR_PATHWAY', active.only=FALSE)
##' ## returned as row indices into vm
##' fids.idxs <- featureIds(gdb, 'c2', value='x.idx')
setGeneric("featureIds", signature="x",
function(x, i, j, value=c('featureId', 'x.id', 'x.idx'),
         active.only=is.conformed(x), ...) {
  standardGeneric("featureIds")
})

##' Fetch the featureIdMap for a \code{GeneSetDb}
##' @exportMethod featureIdMap
setGeneric("featureIdMap", function(x, ...) standardGeneric("featureIdMap"))

setGeneric("featureIdMap<-", function(x, value) {
  standardGeneric("featureIdMap<-")
})


##' Gene Set Collection Metadata
##'
##' @description
##' The design of the GeneSetDb is such that we assume that groups of gene sets
##' are usually defined together and will therefore share similar metadata.
##' These groups of gene sets will fall into the same "collection", and,
##' therefore, metadata for particular gene sets are tracked at the collection
##' level.
##'
##' Types of metadata being referred to could be things like the organism
##' that a batch of gene sets were defined in, the type of feature identifiers
##' that a collection of gene sets are using (ie. \code{EntrezIdentifier}),
##' or a URL pattern that combines the collection,name compound key that one
##' can browse to in order to find out more information about the gene set.
##'
##' There are explicit helper functions that set and get these aforementioned
##' metadata, namely \code{org}, \code{featureIdType},
##' \code{geneSetCollectionURLfunction}, and \code{geneSetURL}. Aribtrary
##' metadata can be stored at the collection level using the
##' \code{addCollectionMetadata} function. More details are provided below.
##'
##' @exportMethod collectionMetadata
##' @rdname collectionMetadata
##' @param x Object to extract the collectionMetadata from
##' @param collection The geneset collection to to query
##' @param name The name of the metadata variable to get the value for
##'
##' @examples
##' gdb <- getMSigGeneSetDb('h')
##'
##' ## Gene Set URLs
##' geneSetURL(gdb, 'h', 'HALLMARK_ADIPOGENESIS')
##' geneSetURL(gdb, c('h', 'h'),
##'            c('HALLMARK_ADIPOGENESIS', 'HALLMARK_ANGIOGENESIS'))
##'
##' ## FeatureId TYpe
##' featureIdType(gdb, 'h')
##'
##' ## Organism
##' org(gdb, 'h')
##'
##' ## Arbitrary metadata
##' gdb <- addCollectionMetadata(gdb, 'h', 'foo', 'bar')
##' cmh <- collectionMetadata(gdb, 'h') ## print this to see
setGeneric("collectionMetadata", signature=c("x", "collection", "name"),
function(x, collection, name, ...) {
  standardGeneric("collectionMetadata")
})

# ##' @export
# ##' @rdname collectionMetadata
# setGeneric("collectionMetadata<-", signature=c("x", "collection", "name"),
# function(x, collection, name, value) {
#   standardGeneric("collectionMetadata<-")
# })

##' @section Gene Set URLs:
##'
##' A URL function can be defined per collection that takes the collection,name
##' compound key and generates a URL for the gene set that the user can browse
##' to for futher information. For instance, the
##' \code{geneSetCollectionURLfunction} for the MSigDB collections are defined
##' like so:
##'
##' \preformatted{
##'   url.fn <- function(collection, name) {
##'     url <- 'http://www.broadinstitute.org/gsea/msigdb/cards/%s.html'
##'     sprintf(url, name)
##'   }
##'   gdb <- getMSigGeneSetDb('h')
##'   geneSetCollectionURLfunction(gdb, 'h') <- url.fn
##' }
##'
##' In this way, a call to \code{geneSetURL(gdb, 'h', 'HALLMARK_ANGIOGENESIS')}
##' will return
##' \url{http://www.broadinstitute.org/gsea/msigdb/cards/HALLMARK_ANGIOGENESIS.html}.
##'
##' This function is vectorized over \code{i} and \code{j}
##'
##' @exportMethod geneSetURL
##' @inheritParams featureIds
##' @rdname collectionMetadata
##'
##' @return A character vector of URLs for each of the genesets identified by
##'   \code{i, j}. \code{NA} is returned for genesets \code{i,j} that are not
##'   found in \code{x}.
setGeneric("geneSetURL", signature="x", function(x, i, j, ...) {
  standardGeneric("geneSetURL")
})

##' @exportMethod geneSetCollectionURLfunction
##' @rdname collectionMetadata
setGeneric("geneSetCollectionURLfunction", signature="x", function(x, i, ...) {
  standardGeneric("geneSetCollectionURLfunction")
})

##' @export
##' @rdname collectionMetadata
setGeneric("geneSetCollectionURLfunction<-", signature="x", function(x, i, value) {
  standardGeneric("geneSetCollectionURLfunction<-")
})

##' @section Feature ID Types:
##'
##' When defining a set of gene sets in a collection, the identifiers used must
##' be of the same type. Most often you'll probably be working with Entrez
##' identifiers, simply because that's what most of the annotations work with.
##'
##' As such, you'd define that your collection uses geneset identifiers like
##' so:
##'
##' \preformatted{
##'   gdb <- getMSigGeneSetDb('h')
##'   featureIdType(gdb, 'h') <- EntrezIdentifier()
##'   ## or, equivalently (but you don't want to use this)
##'   gdb <- addCollectionMetadata(gdb, 'h', 'id_type', EntrezIdentifier())
##' }
##'
##' @exportMethod featureIdType
##' @rdname collectionMetadata
##' @inheritParams featureIds
setGeneric("featureIdType", signature="x", function(x, i, ...) {
  standardGeneric("featureIdType")
})

##' @export
##' @rdname collectionMetadata
setGeneric("featureIdType<-", signature="x", function(x, i, value) {
  standardGeneric("featureIdType<-")
})

##' @section Organism:
##'
##' You're going to want to keep track of the organism the experiments were run
##' in that were used to define this collection of gene sets.
##'
##' \preformatted{
##'   gdb <- getMSigGeneSetDb('h')
##'   org(gdb, 'h') <- 'Homo_sapiens'
##' }
##'
##' @exportMethod org
##' @rdname collectionMetadata
##' @inheritParams featureIds
setGeneric("org", signature="x", function(x, i, ...) {
  standardGeneric("org")
})

##' @export
##' @rdname collectionMetadata
setGeneric("org<-", signature="x", function(x, i, value) {
  standardGeneric("org<-")
})


##' Fetches information for a gene set
##'
##' @description
##' Gene sets inside a \code{GeneSetDb} are indexed by their collection,name
##' compound key. There is no special class to represent an individual gene set.
##' Instead, gene sets are returned as a data.frame, the rows of which enumerate
##' the features that belong to them.
##'
##' When \code{x} is a \code{\link{MultiGSEAResult}}, this function will append
##' the differential expression statistics for the individual features generated
##' across the contrast that defined \code{x}.
##'
##' @exportMethod geneSet
##' @rdname geneSet
##'
##' @inheritParams featureIds
##' @param with.feature.map If \code{TRUE}, then deatils of the feature mapping
##'   from the original featureId space to the target feature space are included
##'   (default: \code{FALSE}).
##' @param ... passed down to inner functinos
##' @template asdt-param
##' @return a data.(frame|table) of gene set information. If \code{x} is a
##'   \code{MultiGSEAResult} object, then differential expression statistics
##'   are added as columns to this result.
setGeneric("geneSet", signature="x", function(x, i, j, ...) {
  standardGeneric("geneSet")
})

##' Fetch the active (or all) gene sets from a GeneSetDb or MultiGSEAResult
##'
##' @export
##' @inheritParams featureIds
##' @template asdt-param
##' @return a data.table with geneset information.
##'
##' @rdname geneSets
##' @exportMethod geneSets
##'
##' @examples
##' gdb <- exampleGeneSetDb()
##' gs <- geneSets(gdb)
setGeneric("geneSets", function(x, ...) standardGeneric('geneSets'))


##' Summarize geneset:feature relationships for specified set of features
##'
##' This function creates a geneset by feature table with geneset membership
##' information for a specified feature set. Only the gene sets that have
##' any of the \code{features} are included.
##'
##' @rdname geneSetSummaryByGenes
##' @exportMethod geneSetSummaryByGenes
##'
##' @param x \code{GeneSetDb} or \code{MultiGSEAResult}
##' @param features a character vector of featureIds
##' @param with.features Include columns for \code{features}? If \code{x} is
##'  is a \code{GeneSetDb}, these columns are \code{TRUE}/\code{FALSE}. If
##'  \code{x} is a \code{MultiGSEAResult} object, the values are the logFC of
##'  the feature if present in the gene set, otherwise its \code{NA}.
##' @param feature.rename if \code{NULL}, the feature columns are prefixed with
##'   \code{featureId_}, if \code{FALSE}, no renaming is done. If \code{x} is
##'   a \code{MultiGSEAResult}, then this can be the column name found in
##'   \code{logFC(x)}, in which case the value for the feature from the given
##'   column name would be used (setting this to \code{"symbol"}) would be a
##'   common thing to do, for instance.
##' @return a data.frame of geneset <-> feature incidence/feature matrix.
##'
##' @examples
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- conform(exampleGeneSetDb(), vm)
##' mg <- multiGSEA(gdb, vm, vm$design, 'tumor', methods=NULL)
##' features <- c("55839", "8522", "29087")
##' gsm.hit <- geneSetSummaryByGenes(gdb, features)
##' gsm.fid <- geneSetSummaryByGenes(mg, features, feature.rename=NULL)
##' gsm.sym <- geneSetSummaryByGenes(mg, features, feature.rename='symbol')
setGeneric("geneSetSummaryByGenes", signature=c("x"),
           function(x, features, with.features=TRUE, feature.rename=NULL, ...)
             standardGeneric("geneSetSummaryByGenes"))

##' @rdname conform
##' @exportMethod conform
setGeneric("conform", function(x, ...) standardGeneric("conform"))

##' Subset a GeneSetDb to only include geneSets with specified features.
##'
##' @exportMethod subsetByFeatures
##' @rdname subsetByFeatures
##'
##' @param x \code{GeneSetDb}
##' @param featureIds Character vector of featureIds
##' @return A subset of \code{x} which contains only the geneSets that contain
##'   features found in \code{featureIds}
##'
##' @examples
##' gdb <- exampleGeneSetDb()
##' features <- c("55839", "8522", "29087")
##' (gdb.sub <- subsetByFeatures(gdb, features))
setGeneric("subsetByFeatures", signature="x",
function(x, features, value=c('featureId', 'x.id', 'x.idx'), ...) {
  standardGeneric("subsetByFeatures")
})

##' @exportMethod unconform
##' @rdname conform
setGeneric("unconform", function(x, ...) standardGeneric("unconform"))

##' Summarizes different results into tabular form
##'
##' @exportMethod summarized
setGeneric("summarized", function(x, ...) standardGeneric("summarized"))
