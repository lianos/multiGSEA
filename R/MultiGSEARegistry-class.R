MultiGSEARegsitry <- function() {
  .MultiGSEARegistry()
}

registerGSEAmethod <- function(method, validate.fn, do.fn) {

}

GSEAmethods <- function() {
  .MGR@methods
}

is.registered <- function(method) {
  method %in% GSEAmethods()
}


doGSEA <- function(method, gsd, x, design, contrast, outdir, use.cache, ...) {

}

setRefClass("MultiGSEARegistry",
  fields=list(methods='character', validate='list', do='list'),
  methods=list(
    initialize=function(...) {
      ns <- asNamespace('multiGSEA')
      gsea.methods <- ls(ns, pattern='^do\\.')
      ## Which ones are implemented
      implemented <- sapply(gsea.methods, function(fn) {
        err <- tryCatch(ns[[fn]](), error=function(e) geterrmessage())
        length(grep('not.*implement', err, ignore.case=TRUE)) == 0
      })

      meths <- gsea.methods[implemented]
      .self$methods <- sapply(meths, function(method) {
        getFunction(paste0('do.', method))
      }, simplify=FALSE)
      .self$validate <- sapply(meths, function(method) {
        getFunction(paste0('validate.inputs.', method))
      }, simplify=FALSE)
    },
    methods=function() {
      .self$methods
    },
    register=function(method, validate, do, override=FALSE) {
      if (method %in% .self$methods && !override) {
        stop(sprintf("`%s` method already defined in registry", method))
      }
      if (!is.function(validate)) {
        stop("`validate` must be a function")
      }
      if (!is.function(do)) {
        stop("`do` must be a function")
      }
    },
    validateInputs=function(method, ...) {

    },
    do=function(method, ...) {

    })
)

