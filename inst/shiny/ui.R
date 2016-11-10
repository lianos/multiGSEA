library(miniUI)

dashboardPage(

  dashboardHeader(title="GSEA Explorer"),

  dashboardSidebar(
    fileInput("mgresult", 'multiGSEA Result'),

    sidebarMenu(
      menuItem("Gene Set View", tabName="genesetView"),
      menuItem("Gene View", tabName="geneView")
    )
  ),

  dashboardBody(
    tags$head(
      tags$link(rel="stylesheet", type="text/css", href="dashboard.css"),
      tags$link(rel="stylesheet", type="text/css", href="miniUI.css")
    ),

    tabItems(
      ## -----------------------------------------------------------------------
      ## tabName: genesetView
      tabItem(
        tabName="genesetView",
        tags$div(style="margin-bottom: 10px; padding: 5px; background-color: white",
                 title='GSEA Results',
                 uiOutput("gseaMethodSummary")),

        fluidRow(
          box(style="padding: 0", width=5, geneSetContrastViewUI("geneset_viewer")),
          box(width=7, mgTableBrowserUI("mg_table_browser")))
      ), ## tabItem: genesetView

      ## -----------------------------------------------------------------------
      ## tabName: geneView
      tabItem(
        tabName="geneView",
        fluidRow(
          box(width=5, mgVolcanoUI("dge_volcano")),
          box(width=7, DT::dataTableOutput("dge_stats"))
          # box(
          #   p("Gene-centric exploration across genesets coming soon"),
          #   tags$ul(
          #     tags$li("Overlap of genes across genesets"),
          #     tags$li("etc ...")),
          #   rbokehOutput("geneView_volcanoPlot")),
          # box(DT::dataTableOutput("geneView_interGenesetMembers"))
        )
      )

    )  ## dashboardBody::tabItems
  ) ## dashboardBody
)
