##' Turns the provided/requested geneset into a GeneSetTable
##'
##' This will ensure that the geneset membersihp in \code{x} gels with the
##' expression data in \code{data}.
##'
##' @param x The geneset somthing
##' @param data The expression data
##' @param species If \code{x} are MSigDB identifiers (ie. c1, c2, etc) then
##' this is used to specify which organism to fetch these from.
##'
##' @return A GeneSetTable object.
preprocessGeneSets <- function(x, data, species=c('human', 'mouse')) {
  species <- match.arg(species)
  if (is.character(x)) {
    x <- getMSigDBset(x, species)
  }
}

