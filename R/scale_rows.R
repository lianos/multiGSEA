#' Centers and scales the rows of a numeric matrix.
#'
#' This was for two reasons: (1) to avoid the (more commonly used)
#' `t(scale(t(x), ...)` idiom; and (2) to specify what values, columns of x,
#' etc. to use to calculate means and sd's to use in the scaling function.
#'
#' For instance, you might want to subtract the mean of a subset of columns
#' from each row in the matrix (like the columns that come from control samples)
#'
#' @section Transformation based on specific columns:
#' `center` and `scale` can be a logical, character, or numeric-like vector.
#' The flexibility enables the following scenarios:
#'
#' 1. The user can set it to `TRUE` to center all values on the mean of their
#'    row. (`FALSE` does no centering)
#' 2. A (named) vector of values that is a superset of rownames(x). These will
#'    be the values that are subtracted from each row.
#' 3. A logical vector as long as ncol(x). Each value will be centered to the
#'    mean of the values of the columns specified as TRUE.
#' 4. An integer vector, the is the analog of 3 but specifies the columns to
#'    use for centering.
#'
#' @export
#' @param x the matrix-like object
#' @param center Either a logical, character, or numeric-like value that
#'   specifies what to center
#' @param scale Either a logical, characeter, or numeric-like value that
#'   specifies what to scale
#' @return a scaled version of `x`
#' @examples
#' # see tests/testthat/test-scale_rows.R for more examples
#' m <- matrix(rnorm(50, mean = 1, sd = 2), nrow = 5,
#'             dimnames = list(LETTERS[1:5], letters[1:10]))
#' s0 <- scale_rows(m, center = TRUE, scale = FALSE)
#' all.equal(s0, t(scale(t(m), center = TRUE, scale = FALSE)))
#'
#' # mean center rows to a specific group of control samples (columns)
#' ctrl <- sample(colnames(m), 3)
#' s.ctrl <- scale_rows(m, center = ctrl, scale = FALSE)
#' ctrl.means <- Matrix::rowMeans(m[, ctrl])
#' all.equal(s.ctrl, t(scale(t(m), center = ctrl.means, scale = FALSE)))
scale_rows <- function(x, center = TRUE, scale = TRUE, ...) {
  UseMethod("scale_rows", x)
}

#' @noRd
#' @export
scale_rows.default <- function(x, center = TRUE, scale = TRUE, ...) {
  # as_matrix will produce either a `matrix`, `Matrix`,
  scale_rows(as_matrix(x, ...), center = center, scale = scale, ...)
}

#' @noRd
#' @export
#' @param base.attributes set to `TRUE` if you want to return the same
#'   attr(x, "scaled:(center|scale)") combo that `base::scale` does. Mostly
#'   for unit test purposes. When `FALSE` (default) if `center` or `scale`
#'   are not triggered, then their respeict "scaled:XXX" attribute is set
#'   to `FALSE`.
scale_rows.matrix <- function(x, center = TRUE, scale = TRUE, ...,
                              base.attributes = FALSE) {
  center. <- .target_column_idxs(x, center)
  if (test_integerish(center., lower = 1, upper = ncol(x))) {
    center. <- rowMeans(x[, center., drop = FALSE])
  }
  if (!isFALSE(center.)) {
    assert_numeric(center., len = nrow(x))
    x <- x - center.
  }

  scale. <- .target_column_idxs(x, scale)
  if (test_integerish(scale., lower = 1, upper = ncol(x))) {
    scale. <- rowSds(x[, scale., drop = FALSE], center = FALSE)
  }
  if (!isFALSE(scale.)) {
    assert_numeric(scale., len = nrow(x))
    x <- x / scale.
  }

  if (base.attributes) {
    if (!isFALSE(center.)) attr(x, "scaled:center") <- center.
    if (!isFALSE(scale.)) attr(x, "scaled:scale") <- scale.
  } else {
    attr(x, "scaled:center") <- center.
    attr(x, "scaled:scale") <- scale.
  }

  x
}

#' @noRd
#' @export
scale_rows.Matrix <- function(x, center = TRUE, scale = TRUE, ...) {
  if (!is(x, "Mnumeric")) {
    warning("Matrix object is not of type 'Mnumeric', making dense ...",
            immediate. = TRUE)
    x <- as.matrix(x)
  }
  scale_rows.matrix(x, center = center, scale = scale, ...)
}

#' @noRd
#' @export
scale_rows.sparseMatrix <- function(x, center = TRUE, scale = TRUE, ...) {
  warning("Densifying a sparseMatrix", immediate. = TRUE)
  scale_rows(as.matrix(x), center = center, scale = scale, ...)
}

#' Helper function that identifies the value that center/scale should be
#' when it's any of the logical, character, or numeric values that are
#' supported by the `center` and `scale` parameters in [scale_rows()]
#'
#' @noRd
.target_column_idxs <- function(x, target, ...) {
  # User specified row-values to use for centering or scaling
  if (test_numeric(target, min.len = nrow(x)) &&
      !test_integerish(x, lower = 1)) {
    if (test_subset(rownames(x), names(target))) {
      target <- target[rownames(x)]
    }
    if (length(target) == nrow(x)) {
      return(target)
    }
  }

  if (test_logical(target) && length(unique(target)) == 1L) {
    target <- target[1L]
  }
  if (!isFALSE(target)) {
    if (test_character(target)) {
      assert_subset(target, colnames(x))
      target <- colnames(x) %in% target
    }
    if (isTRUE(target)) {
      target <- rep(TRUE, ncol(x))
    }
    if (test_logical(target) && length(target) == ncol(x)) {
      # indicator of which columns to calculate mean from and recenter to
      target <- which(target)
    }
    assert_integerish(target, lower = 1, upper = ncol(x),
                      unique = TRUE, any.missing = FALSE)
  }

  target
}
