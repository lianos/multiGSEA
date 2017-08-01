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
##' @importFrom Matrix rankMatrix
##' @export
##'
##' @param x The expression object to use
##' @param design A design matrix, if the GSEA method(s) require it
##' @param contrast A contrast vector (if the GSEA method(s) require it)
##' @param A character vector of the GSEA methods that these inputs will be used
##'   for
##' @param require.x.rownames Leave this alone, should always be \code{TRUE} but
##'   have it in this package for dev/testing purposes.
##' @param ... other variables that called methods can check if they want
##' @return A list with "normalized" versions of \code{$x}, \code{$design}, and
##'   \code{$contrast} for downstream use.
validateInputs <- function(x, design=NULL, contrast=NULL, methods=NULL,
                           require.x.rownames=TRUE, ...) {
  if (is.character(methods)) {
    .unsupportedGSEAmethods(methods)
  } else if (!is.null(methods)) {
    stop("Illegal type for `methods`: ", class(methods)[1L])
  }

  if (is(x, 'DGEList') && !disp.estimated(x)) {
    stop("It does not look like estimateDisp has been run on DGEList")
  }

  if (is.vector(x)) {
    x <- matrix(x, ncol=1L, dimnames=list(names(x), NULL))
  }

  ## Check that x is generally OK
  x.kosher <- validate.X(x)
  if (!isTRUE(x.kosher)) {
    stop("Bad expression object x provided: ", paste(x.errs, collapse=','))
  }

  ## Validate the input expression object separately (not sure why now)
  if (!is.null(methods)) {
    is.valid.x <- sapply(methods, function(meth) {
      fn <- getFunction(paste0('validate.x.', meth))
      fn(x)
    }, simplify=FALSE)
    bad <- which(sapply(is.valid.x, Negate(isTRUE)))
    if (length(bad)) {
      msg <- paste("Error validating x for methods:",
                   paste(methods[bad], collapse=', '))
      stop(msg)
    }
  } else {
    if (!inherits(x, .valid.x)) {
      stop("Invalid expression object (x) type: ", class(x)[1L])
    }
    if (!is.character(rownames(x)) && require.x.rownames) {
      stop("The expression object does not have rownames ...")
    }
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

  ## method specific validation checks
  if (is.character(methods)) {
    errs.all <- sapply(methods, function(method) {
      fn <- getFunction(paste0('validate.inputs.', method))
      errs <- fn(x, design, contrast, ...)
    }, simplify=FALSE)

    errs.un <- unlist(errs.all)
    if (length(errs.un)) {
      msg <- paste("Errors in inputs:\n    *",
                   paste(names(errs.un), collapse='\n    * '))
      msg <- paste(msg, '=======', unname(errs.un), sep='\n')
      stop(msg)
    }
  }

  list(x=x, design=design, contrast=contrast, is.full.design=is.matrix(design))
}

##' Checkes that there are no NAs in x
na.check <- function(x) {
  if (is(x, 'DGEList')) x <- x$counts
  if (is(x, 'EList')) x <- x$E
  if (is(x, 'ExpressionSet')) x <- exprs(x)
  if (is(x, 'SummarizedExperiment')) x <- assay(x)
  if (any(is.na(x))) {
    stop("No NA's allowed in x")
  }
}

##' Checks a DGEList to see if estimateDisp() was run on it
##' @param x Input DGEList
##' @return TRUE if yes, FALSE if no or if x is not a DGEList.
disp.estimated <- function(x) {
  ## check that estimateDisp has been run
  if (!is(x, "DGEList")) return(FALSE)
  reqd <- c('common.dispersion', 'trended.dispersion', 'tagwise.dispersion')
  for (wut in reqd) {
    vals <- x[[wut]]
    if (is.null(vals) || is.na(vals) || !is.numeric(vals)) {
      warning(sprintf('[[%s]] was not found, did you call estimateDisp?"', wut),
              immediate.=TRUE)
      return(FALSE)
    }
  }
  return(TRUE)
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
    if (rankMatrix(design) != ncol(design)) {
      errs$design.not.full.rank <- TRUE
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

    if (abs(sum(contrast)) > 1e-5) {
      warning("Sum of contrast vector != 0", immediate.=TRUE)
      # errs$sum.contrast.not0 <- TRUE
      # return(errs)
    }
  }

  if (ret.err.only || length(errs)) errs else contrast
}

.validate.inputs.full.design <- function(x, design, contrast,
                                         require.x.rownames=FALSE, ...) {
  errs <- list()
  if (!inherits(x, .valid.x)) {
    errs <- c(errs,
              sprintf("Invalid expression object (x) type: %s", class(x)[1L]))
  } else {
    if (!is.character(rownames(x)) && require.x.rownames) {
      errs <- c(errs, "The expression object does not have rownames ...")
    }
  }
  if (ncol(x) == 1L) {
    errs <- c(errs, 'expression matrix needs more than one column')
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

.validate.inputs.logFC.only <- function(x, design, contrast, ...) {
  errs <- list()
  if (ncol(x) > 1) {
    errs <- .validate.inputs.full.design(x, design, contrast, ...)
  } else {
    if (is(x, 'DGEList')) {
      errs$DGEList.not.supported.for.gsd <- TRUE
    }
  }
  errs
}

## Some GSEA functions can work on a simple (named) pre-ranked vector of
## logFC's or t-statistcs, like limma::geneSetTest or fgsea, for example.
## The usre shoudl be able to simply provide such a pre-ranked vector, or can
## provide a "full.design" set of inputs from which logFC's or t-statistics
## could be computed using multiGSEA's internal calculateIndividualLogFC
## function.
.validate.inputs.preranked <- function(x, design, contrast, ...) {
  if (is.vector(x) || is.matrix(x) && ncol(x) == 1L) {
    .validate.inputs.logFC.only(x, design, contrast, ...)
  } else {
    .validate.inputs.full.design(x, design, contrast, ...)
  }
}

## Validation Methods for Expression Objects -----------------------------------

##' importFrom edgeR DGEList
validate.DGEList <- function(x) {
  if (!isTRUE(is(x, 'DGEList'))) {
    return("x is not a DGEList")
  }
  if (!is.numeric(x$common.dispersion)) {
    return("dispersion is not estimated, minimally call estimateDisp on x")
  }
  TRUE
}

validate.XwithWeights <- function(x) {
  if (!isTRUE(is(x, 'EList') || is(x, 'eSet'))) {
    return("x must be an EList or ExpressionSet")
  }
  if (is(x, 'EList')) {
    if (!isTRUE(dim(x$weights), dim(x))) {
      return("EList does not have weights")
    }
  }
  if (is(x, 'eSet')) {
    if (!'weights' %in% assayDataElementNames(x)) {
      return("weights assay not in eSet x")
    }
  }

  TRUE
}

## We expect x to be matrix like now
validate.X <- function(x) {
  if (!inherits(x, .valid.x)) {
    return("Invalid expression object (x) type: ", class(x)[1L])
  }
  na.check(x)
  if (!is.character(rownames(x))) {
    return("The expression object does not have rownames ...")
  }
  if (any(is.na(rownames(x)))) {
    return("NAs in rownames of x")
  }
  if (any(duplicated(rownames(x)))) {
    return("Duplicated rownames in x")
  }
  if (is(x, 'DGEList')) {
    return(validate.DGEList(x))
  }
  TRUE
}
