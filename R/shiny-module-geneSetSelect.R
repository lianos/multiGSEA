##' A module that creates a dynamic selectizeInput for a MultiGSEAResult object
##'
##' This code was inspired from the
##' \href{https://gist.github.com/MarkEdmondson1234/7e56ee7ac5caa74224327489b0849e61}{dynamicSelectShinyModule.R}
##' gist.
##'
##' @export
##' @rdname geneSetSelectModule
geneSetSelectUI <- function(id, label="Select Gene Set", mg=NULL,
                            server=TRUE) {
  requireNamespace('shiny')
  if (!is.null(mg) && !is(mg, 'MultiGSEAResultContainer')) {
    stop("`mg` must either be NULL or a MultiGSEAResultContainer")
  }
  ns <- shiny::NS(id)
  shiny::uiOutput(ns("geneset_picker"))
}

##' Returns information about the geneSetSelect object
##'
##' collection, name, stats, sep
##'
##' @export
##' @rdname geneSetSelectModule
geneSetSelect <- function(input, output, session, mgc, server=TRUE,
                          maxOptions=Inf, sep='_::_') {
  requireNamespace('shiny')
  ## Programmaticaslly create the UI from the MultiGSEAResults
  output$geneset_picker <- shiny::renderUI({
    shiny::req(mgc())
    mo <- if (is.infinite(maxOptions)) nrow(geneSets(mgc()$mg)) else maxOptions
    gs.render.select.ui(session$ns, mgc()$choices, server=server, maxOptions=mo)
  })

  if (server) {
    shiny::observeEvent(mgc(), {
      ## req(mgc()$choices)
      # sel <- if (is.null(selected)) {
      #   NULL
      # } else {
      #   paste(selected()$collection, selected()$name, sep=sep)
      # }
      updateSelectizeInput(session, "geneset", choices=mgc()$choices,
                           server=TRUE, selected=NULL)
    })
  }

  shiny::reactive({
    shiny::req(input$geneset)
    info <- input$geneset %>%
      strsplit(sep, fixed=TRUE) %>%
      unlist %>%
      sapply(as.character) %>%
      setNames(c('collection', 'name'))
    stats <- arrange_(geneSet(mgc()$mg, info[1L], info[2L]), ~ -logFC)
    list(collection=info[1L], name=info[2L], stats=stats,
         select.id=session$ns('geneset'), sep=sep)
  })
}

##' @export
##' @rdname geneSetSelectModule
updateGeneSetSelect <- function(session, id, label=NULL, choices=NULL,
                                selected=NULL, options=list(), server=FALSE) {
  stopifnot(requireNamespace('shiny'))
  childScope <- session$makeScope(id)
  shiny::withReactiveDomain(childScope, {
    mod.id <- childScope$ns('geneset')
    shiny::updateSelectizeInput(session, mod.id, label=label,
                                choices=choices, selected=selected,
                                options=options,
                                server=server)
  })
}

## Utility Functions -----------------------------------------------------------

##' Builds a selectizeInput widget that is specific to a MultiGSEAResult
##'
##' @rdname geneSetSelectModule
##'
##' @param ns the namespace function for this module
##' @param choices the output of gs.select.choices(MultiGSEAResult)
##' @param server \code{logical} to indicate whether options should be loaded
##'   on the server side (default: \code{TRUE})
##' @param maxOptions The maximum number of options to load into the dropdown
##' @return a properly wired \code{\link[shiny]{selectizeInput}}
gs.render.select.ui <- function(ns, choices, server=TRUE,
                                maxOptions=1000, sep='_::_') {
  requireNamespace('shiny')
  # predefine all options groups
  optgroups = lapply(unique(choices$collection), function(col) {
    list(value=col, label=col)
  })

  # define options to customize the selectize object
  si.opts <- list(
    placeholder='Select Gene Set',
    optgroups=optgroups,
    optgroupField='collection',
    searchField = c('label'),
    maxOptions=maxOptions,
    render=I("{
             option: function(item, escape) {
             return '<div>' + escape(item.label) + '</div>';
             }}"))

  if (server) {
    ui <- shiny::selectizeInput(ns("geneset"), label=NULL, choices=NULL, options=si.opts)
  } else {
    choices <- sapply(unique(choices$collection), function(x) {
      out <- subset(choices, collection == x)
      setNames(out$value, out$label)
    }, simplify=FALSE)
    ui <- shiny::selectizeInput(ns("geneset"), label=NULL, choices=choices)
  }

  ui
}

##' Builds a \code{data.frame} used to populate choices for selectizeInput
##'
##' Note that when returning a data.frame for the choices of a selectizeInput,
##' we need a column called "value" and a column called "label".
##' \itemize{
##'   \item{value}{
##'     the value that is sent back when an item is selected;
##'   }
##'   \item{label}{
##'     the text that appears in the selection after its triggered
##'   }
##' }
##'
##' @rdname geneSetSelectModule
##'
##' @param mg \code{MultiGSEAResult} to build options for
##' @return \code{data.frame} to populate \code{choices} of
##'   \code{selectizeInput}
gs.select.choices <- function(mg, sep='_::_') {
  geneSets(mg) %>%
    select(collection, label=name) %>%
    mutate(value=paste(collection, label, sep=sep)) %>%
    as.data.frame(stringsAsFactors=FALSE)
}
