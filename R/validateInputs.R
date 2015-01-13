## .input.reqs <- c('full-design', 'logFC-only')
## .method.to.input.reqs <- c(camera='full-design',
##                            roast='full-design',
##                            geneSetTest='logFC-only',
##                            hyperGeometricTest='full-design')

##' Validate the input objects to a GSEA call.
##'
##' Checks to ensure that the values for \code{x}, \code{design}, and
##' \code{contrast} are appropriate for the GSEA \code{methods} being used. If
##' they are kosher, then "normalized" versions of these objects are returned in
##' an (aptly) named list, otheerwise an error is thrown.
##'
##' This function is strange in that we both want to verify the objects, and
##' return them in some canonical form, so it is normal for the caller to then
##' use the values for \code{x}, \code{design}, and \code{contrast} that are
##' returned from this call, and not the original values for these objects
##' themselves
##'
##' I know that the validation/checking logic is a bit painful (and repetitive)
##' here. I will (perhaps) clean that up some day.
##'
##' @export
##'
##' @param x The expression object to use
##' @param design A design matrix, if the GSEA method(s) require it
##' @param contrast A contrast vector (if the GSEA method(s) require it)
##' @param A character vector of the GSEA methods that these inputs will be used
##'   for
##' @param require.x.rownames Leave this alone, should always be \code{TRUE} but
##'   have it in this package for dev/testing purposes.
##'
##' @return A list with "normalized" versions of \code{$x}, \code{$design}, and
##'   \code{$contrast} for downstream use.
validateInputs <- function(x, design=NULL, contrast=NULL, methods=NULL,
                           require.x.rownames=TRUE) {
  if (is.character(methods)) {
    .unsupportedGSEAmethods(methods)
  } else if (!is.null(methods)) {
    stop("Illegal type for `methods`: ", class(methods)[1L])
  }

  if (is(x, 'DGEList')) {
    ## TODO: stop if estimateDisp(x, design) was called on input DGEList
    warning("Not checking if estimateDisp() was called", immediate.=TRUE)
  }

  if (is.vector(x)) {
    x <- matrix(x, ncol=1L, dimnames=list(names(x), NULL))
  }

  if (!inherits(x, .valid.x)) {
    stop("Invalid expression object (x) type: ", class(x)[1L])
  }
  if (!is.character(rownames(x)) && require.x.rownames) {
    stop("The expression object does not have rownames ...")
  }

  if (is.matrix(design)) {
    design.errs <- .validateDesign(x, design)
    if (length(design.errs)) {
      stop("Design matrix problems:\n    * ",
           paste(names(design.errs), collapse='\n    * '))
    }
    if (is.null(contrast)) {
      contrast <- ncol(design)
    } else {
      contrast <- .validateContrastVector(contrast, design)
      if (is.list(contrast)) {
        stop("Contrast vector problems:\n    * ",
             paste(names(contrast), collapse='\n    * '))
      }
    }
  }

  if (is.character(methods)) {
    errs.all <- sapply(methods, function(method) {
      fn <- getFunction(paste0('validate.inputs.', method))
      errs <- fn(x, design, contrast)
    }, simplify=FALSE)

    errs.un <- unlist(errs.all)
    if (length(errs.un)) {
      stop("Erros in inputs:\n    * ",
           paste(names(errs.un), collapse='\n    * '))
    }
  }

  list(x=x, design=design, contrast=contrast)
}

## Returns a 0-length list when there is no error
.validateDesign <- function(x, design) {
  errs <- list()
  if (!is.matrix(design)) {
    errs$design.not.matrix <- TRUE
  } else {
    if (nrow(design) != ncol(x)) {
      errs$design.discordant.dims <- TRUE
    }
    if (!is.character(colnames(design))) {
      errs$design.no.colnames <- TRUE
    }
  }
  errs
}

## When there is no error, returns a valid contrast vector, otherwise will
## return a named list of errors.
##
## I know this is a bad design!
.validateContrastVector <- function(contrast, design, ret.err.only=FALSE) {
  errs <- list()
  if (!is.vector(contrast) || !(is.character(contrast)||is.numeric(contrast))) {
    errs$illegal.contrast.type <- TRUE
    return(errs)
  }

  if (length(contrast) == 1L) {
    if (is.character(contrast)) {
      contrast <- which(colnames(design) == contrast)
      if (length(contrast) != 1L) {
        errs$contrast.name.not.found <- TRUE
        return(errs)
      }
    } else {
      if (as.integer(contrast) != contrast) {
        errs$numeric.contrast.not.integer <- TRUE
        return(errs)
      }
      if (contrast < 1 || contrast > ncol(design)) {
        errs$illegal.contrast.bound <- TRUE
        return(errs)
      }
    }
  } else {
    if (!is.numeric(contrast)) {
      errs$long.contrast.vector.not.numeric <- TRUE
      return(errs)
    }

    if (length(contrast) != ncol(design)) {
      errs$illegal.contrast.length <- TRUE
      return(errs)
    }

    if (sum(contrast) != 0) {
      errs$sum.contrast.not0 <- TRUE
      return(errs)
    }
  }

  if (ret.err.only || length(errs)) errs else contrast
}

.validate.inputs.full.design <- function(x, design, contrast) {
  errs <- list()
  if (!inherits(x, .valid.x)) {
    errs <- c(errs,
              sprintf("Invalid expression object (x) type: %s", class(x)[1L]))
  } else {
    if (!is.character(rownames(x)) && require.x.rownames) {
      errs <- c(errs, "The expression object does not have rownames ...")
    }
  }
  if (!is.matrix(design)) {
    errs$design.matrix.required <- TRUE
  } else {
    errs <- c(errs, .validateDesign(x, design))
  }
  if (!is.vector(contrast)) {
    errs$contrast.vector.required <- TRUE
  } else {
    errs <- c(errs, .validateContrastVector(contrast, design, TRUE))
  }
  errs
}

.validate.inputs.logFC.only <- function(x, design, contrast) {
  errs <- list()
  if (ncol(x) > 1) {
    errs <- .validate.inputs.full.design(x, design, contrast)
  } else {
    if (is(x, 'DGEList')) {
      errs$DGEList.not.supported.for.gsd <- TRUE
    }
  }
  errs
}
