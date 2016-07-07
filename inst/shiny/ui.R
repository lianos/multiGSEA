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
      tags$link(rel="stylesheet", type="text/css", href="rmd_plugins.css")
    ),
    tabItems(
      ## -----------------------------------------------------------------------
      ## tabName: genesetView
      tabItem(
        tabName="genesetView",
        fluidRow(
          box(
            title="Gene Set Expression Profile",
            ## selectInput("genesetView_plot_type", "Plot Type", c('density')),
            ##h3("Gene Set Expression Profile"),
            plotlyOutput("genesetView_gseaPlot")),
          box(
            title="Gene Set Membership",
            ## h3("Gene Set Membership"),
            DT::dataTableOutput("genesetView_genesetMembers"))
          ##, box(DT::dataTableOutput("genesetView_interGenesetMembers"))
        ),

        box(
          width=12,
          title='GSEA Results',
          ## h3("GSEA Results"),
          fluidRow(
            column(3, selectInput("gseaMethod", "GSEA Methods", "")),
            column(4, sliderInput("gseaReportFDR", "FDR Cutoff", min=0, max=1, value=0.2, step=0.05))
          ),
          uiOutput("gseaResultName"),
          uiOutput("gseaMethodSummary"),
          h4("Gene Set Statistics"),
          DT::dataTableOutput("gseaResultTable")
        )
      ), ## tabItem: genesetView

      ## -----------------------------------------------------------------------
      ## tabName: genesView
      tabItem(
        tabName="geneView",
        fluidRow(
          box(plotlyOutput("geneView_volcanoPlot")),
          box(DT::dataTableOutput("geneView_interGenesetMembers"))
        )
      )

    )  ## dashboardBody::tabItems
  ) ## dashboardBody
)
