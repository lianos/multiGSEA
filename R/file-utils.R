##' Utility function to determine if paths listed are real directories
##'
##' @param ... character vectors containing file paths. Tilde-expansion is
##' done: see 'path.expand'.
##'
##' @return A logical vector indicating whether or not the queried directories
##' exist.
dir.exists <- function(...) {
  sapply(file.info(...)$isdir, isTRUE)
}

##' Checks that the directory provided is writable by the current user.
##'
##' This works by testing to put a temporary file into an already existing
##' directory
##'
##' @param path The path to a directory to check.
##'
##' @return \code{logical}, \code{TRUE} if \code{path} is writable by the
##' current user, otherise \code{FALSE}
dir.writable <- function(path) {
  if (!dir.exists(path)) {
    stop("The directory provided does not exist: ", path)
  }
  tmp.fn <- tempfile(tmpdir=path, fileext='.test.tmp')
  on.exit({
    if (file.exists(tmp.fn)) unlink(tmp.fn)
  })

  tryCatch({
    suppressWarnings(writeLines('test', tmp.fn))
    TRUE
  }, error=function(e) FALSE)
}
