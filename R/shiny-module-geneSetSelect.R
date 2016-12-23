##' A module that creates a dynamic selectizeInput for a MultiGSEAResult object
##'
##' This code was inspired from the
##' \href{https://gist.github.com/MarkEdmondson1234/7e56ee7ac5caa74224327489b0849e61}{dynamicSelectShinyModule.R}
##' gist.
##'
##' @export
##' @importFrom shiny NS uiOutput
##' @rdname geneSetSelectModule
geneSetSelectUI <- function(id, label="Select Gene Set") {
  ns <- NS(id)
  uiOutput(ns("geneset_picker"))
}

##' @section Module Return:
##' Returns information about the geneSetSelect object
##'
##' @export
##' @rdname geneSetSelectModule
##' @importFrom shiny renderUI req outputOptions observeEvent reactive
##' @importFrom shiny updateSelectizeInput
##'
##' @param input,output,session the shiny-required bits for the module
##' @param mgc A \code{\link{MultiGSEAResultContainer}} object
##' @param server boolean to indicate whether the genesets in the geneSetSelect
##'   widget should be rendered server side or not (Default: \code{TRUE})
##' @param maxOptions a paremeter used to customize the
##'   \code{GeneSetSelect::selectizeInput} UI element. I thought one might want
##'   to tweak this, but I just leave it as is.
##' @param sep the separater to put between the collection and name bits of a
##'   geneset. These are the values used in the gene set \code{selectizeInput}.
geneSetSelect <- function(input, output, session, mgc, server=TRUE,
                          maxOptions=Inf, sep='_::_') {
  ## Programmaticaslly create the UI from the MultiGSEAResults
  output$geneset_picker <- renderUI({
    req(mgc())
    mo <- if (is.infinite(maxOptions)) nrow(geneSets(mgc()$mg)) else maxOptions
    gs.render.select.ui(session$ns, mgc()$choices, server=server, maxOptions=mo)
  })
  outputOptions(output, "geneset_picker", suspendWhenHidden=FALSE)

  if (server) {
    observeEvent(mgc(), {
      updateSelectizeInput(session, "geneset", choices=mgc()$choices,
                           server=TRUE, selected=NULL)
    })
  }

  vals <- reactive({
    gs <- input$geneset
    if (is.null(gs) || length(gs) == 0 || nchar(gs) == 0) {
      ## HACK, just but something here if it's not selectd
      ## gs <- mgc()$choices$value[1L]
      coll <- name <- stats <- NULL
    } else {
      info <- gs %>%
        strsplit(sep, fixed=TRUE) %>%
        unlist %>%
        sapply(as.character) %>%
        setNames(c('collection', 'name'))
      coll <- info[1L]
      name <- info[2L]
      stats <- arrange_(geneSet(mgc()$mg, info[1L], info[2L]), ~ -logFC)
    }

    list(collection=coll, name=name, stats=stats,
         select.id=session$ns('geneset'), sep=sep)
  })
  return(vals)
}

##' @export
##' @rdname geneSetSelectModule
##' @importFrom shiny updateSelectizeInput
updateGeneSetSelect <- function(session, id, label=NULL, choices=NULL,
                                selected=NULL, options=list(), server=FALSE) {
  childScope <- session$makeScope(id)
  shiny::withReactiveDomain(childScope, {
    mod.id <- childScope$ns('geneset')
    updateSelectizeInput(session, mod.id, label=label,
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
##' @importFrom shiny selectizeInput
##' @param ns the namespace function for this module
##' @param choices the output of gs.select.choices(MultiGSEAResult)
##' @param server \code{logical} to indicate whether options should be loaded
##'   on the server side (default: \code{TRUE})
##' @param maxOptions The maximum number of options to load into the dropdown
##' @return a properly wired \code{\link[shiny]{selectizeInput}}
gs.render.select.ui <- function(ns, choices, server=TRUE,
                                maxOptions=1000, sep='_::_') {
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
    ui <- selectizeInput(ns("geneset"), label=NULL, choices=NULL, options=si.opts)
  } else {
    choices <- sapply(unique(choices$collection), function(x) {
      out <- subset(choices, collection == x)
      setNames(out$value, out$label)
    }, simplify=FALSE)
    ui <- selectizeInput(ns("geneset"), label=NULL, choices=choices)
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
