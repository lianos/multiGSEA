##' Trim an input expressionset and create matching gsea table to use in
##' downstream GSEA analyses.
##'
##' @importFrom matrixStats rowVars
##' @export
##'
##' @param x An "Expression Object" to run GSEA stuff over
##' @param gene.sets Either a data.table generated from AnalyzeThis::prepMsigDB
##' or a named list of index vectors into \code{x} that define gene sets.
##' @param x2gsid A named character vector which maps the rownames of \code{x}
##' to the gene ids in the \code{gene.sets}. The values of this vector are
##' the IDs used in \code{gene.sets} (typically entrez.id) and the \code{names}
##' are the \code{rownames} of \code{x} (the features in x)
##' @param min.gs.size The minimum number of genes required in a geneset to
##' include it in testing.
##' @param unique.by If we need to map the rownames of x to unique entrez.id
##' this is the criteria used to select the row that "wins" to be the
##' representative feature for that entrez.id. The feature with the maximum
##' mean or variance is picked.
##'
##' @return a list with a "trimmed" \code{x} object that matches the vectors
##' in the AnalyzeThis::gs.table stored in \code{gs.table}
init.gsea.data <- function(x, gene.sets, x2gsid=NULL, min.gs.size=5,
                           unique.by=c('mean', 'var'),
                           species=c('human', 'mouse')) {
  unique.by <- match.arg(unique.by)
  species <- match.arg(species)

  entrez.fn <- switch(class(x)[1L],
                      matrix=function(x) sub('GeneID:', '', rownames(x)),
                      ExpressionSet=function(x) fData(x)$entrez.id,
                      EList=function(x) x$genes$entrez.id,
                      DGEList=function(x) x$genes$entrez.id)

  y <- limma::getEAWP(x)
  xdebug <- list(x=x, gene.sets=gene.sets, x2gsid=x2gsid)

  ## Create a geneset data.table and a geneset list that camera will use
  if (is.character(gene.sets)) {
    valid.gs <- paste0('c', 1:7)
    if (length(bad.gs <- setdiff(gene.sets, valid.gs))) {
      assign('.rmd.camera.debug', xdebug, .GlobalEnv)
      stop("gene.sets defined by name must be one of the msigdb curated sets\n.",
           "The following gene.sets could not be found: ",
           paste(bad.gs, collapse=','))
    }

    if (!all(.is.entrez(rownames(x)))) {
      if (is.null(x2gsid)) {
        warning("rownames(x) do not look like entrez ids and x2gsid is missing, ",
                "trying to autogenerate it myself")
        entrez <- entrez.fn(x)
        if (is.null(entrez) ||
            !all(.is.entrez(entrez[!is.na(entrez)])) ||
            sum(!is.na(entrez)) < 100) {
          assign('.rmd.camera.debug', xdebug, .GlobalEnv)
          stop("Could not guess entrez ids for input `x` object")
        }
        x2gsid <- setNames(entrez, rownames(x))
      }
    }

    if (!is.null(x2gsid) && any(duplicated(x2gsid))) {
      ## If we have to map features to entrez.ids, we need to ensure that there
      ## is only one unique feature per entrez.id -- which means we may need to
      ## remove some rows from x
      xx <- data.table(id=rownames(y$exprs),
                       mean=rowMeans(y$exprs),
                       var=rowVars(y$exprs))
      xx[, entrez.id := x2gsid[id]]
      keep <- !is.na(xx$entrez.id)
      if  (!all(keep)) {
        warning(sum(!keep), " features are being removed from x due to NA ids")
      }
      xx <- xx[keep]
      setkeyv(xx, c('entrez.id', unique.by))
      take <- xx[, list(idx=.I[.N], id=id[.N]), by='entrez.id']

      message("Keeping ", nrow(take), " / ", nrow(x), " features")
      x <- x[take$id,]
      x2gsid <- x2gsid[rownames(x)]
      xdebug$x2gsid <- x2gsid
    }

    sig.db <- AnalyzeThis::prepMsigDB(x, gene.sets, x2gs=x2gsid, species=species)
    ## gene.sets <- setNames(sig.db$membership, sig.db$id)
  } else if (AnalyzeThis::is.gs.table(gene.sets)) {
    x <- x[attr(gene.sets, 'gene.id'),]
    sig.db <- gene.sets
    ## gene.sets <- setNames(sig.db$membership, sig.db$id)
  } else if (is.list(gene.sets)) {
    ids <- names(gene.sets)
    if (is.null(ids) || !all(ids == make.unique(ids))) {
      assign('.rmd.camera.debug', xdebug, .GlobalEnv)
      stop("A list `gene.sets` object needs unique `names()`")
    }

    ## Check each genest index vector in the `gene.sets` list to see if they
    ## are the right length or their index values are kosher wrt `x`
    ## kosher <- sapply(1:length(gene.sets), function(i) {
    for (i in 1:length(gene.sets)) {
      gs <- gene.sets[[i]]
      if (is.logical(gs)) {
        if (length(gs) != nrow(x)) {
          assign('.rmd.camera.debug', xdebug, .GlobalEnv)
          msg <- paste0("gene.set[[", i, "]] is wrong length",
                        "(see .rmd.camera.debug variable in .GlobalEnv)")
          stop(msg)
        }
        ## return(TRUE)
        next
      }
      if (is.numeric(gs)) {
        idx <- as.integer(gs)
        if (!all(idx == gs)) {
          assign('.rmd.camera.debug', xdebug, .GlobalEnv)
          msg <- paste0("gene.set[[", i, "]] indices are not integers ")
          stop(msg)
        }
        if (any(idx < 1) || any(idx > nrow(x))) {
          assign('.rmd.camera.debug', xdebug, .GlobalEnv)
          msg <- paste0("gene.set[[", i, "]] indices are out-of-bounds")
          stop(msg)
        }
        ## Convert this to a logical vector
        gs.logical <- logical(nrow(x))
        gs.logical[idx] <- TRUE
        gene.sets[[i]] <- gs.logical
        ## return(TRUE)
        next
      }
      if (is.character(gs)) {
        idx <- match(gs, rownames(x))
        if (any(is.na(idx))) {
          assign('.rmd.camera.debug', xdebug, .GlobalEnv)
          msg <- paste0("gene.set[[", i, "]] indices have names not in ",
                        "rownames(x)")
          stop(msg)
        }
        gs.logical <- logical(nrow(x))
        gs.logical[idx] <- TRUE
        gene.sets[[i]] <- gs.logical
        ## return(TRUE)
        next
      }
      msg <- paste0("gene.set[[", i, "]] has indices of unusable type: ",
                    class(gs)[1L])
      stop(msg)
    }
    sig.db <- data.table(id=names(gene.sets), N=sapply(gene.sets, sum))
    sig.db[, n := N]
    sig.db[, membership := unname(gene.sets)]
    attr(sig.db, 'gene.id') <- rownames(x)
  } else {
    assign('.rmd.camera.debug', xdebug, .GlobalEnv)
    stop("Don't know how to parse the gene.sets passed in of type: ",
         paste(class(gene.sets), collapse=','))
  }

  list(x=x, gs.table=sig.db)
}
