.valid.x <- c('matrix', 'eSet', 'EList', 'DGEList', 'SummarizedExperiment')

## This is painful error checking logic.
## validateInputs will throw an error if there are problems found in the inputs
## or will return a "cleaned" version of the inputs on no error.
validateInputs <- function(x, design=NULL, contrast=NULL,
                           methods=if (is.matrix(design)) 'camera' else NULL,
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
    stop("The expression object has nor rownames ...")
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
      fn <- getFunction(paste0('.validateInputs.', method))
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

.validateInputs.camera <- function(x, design, contrast) {
  errs <- list()
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

.validateInputs.roast <- function(x, design, contrast) {
  .validateInputs.camera(x, design, contrast)
}

.validateInputs.gst <- function(x, design, contrast) {
  errs <- list()
  if (is(x, 'DGEList')) {
    errs$DGEList.not.supported.for.gst <- TRUE
  }
  if (ncol(x) > 1) {
    errs <- c(errs, .validateInputs.camera(x, design, contrast))
  }
  errs
}
