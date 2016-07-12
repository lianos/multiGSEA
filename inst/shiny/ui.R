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
      tags$link(rel="stylesheet", type="text/css", href="dashboard.css")
    ),
    tabItems(
      ## -----------------------------------------------------------------------
      ## tabName: genesetView
      tabItem(
        tabName="genesetView",
        fluidRow(
          box(
            title="Gene Set Expression Profile",
            plotlyOutput("genesetView_gseaPlot"),
            fluidRow(
              column(8,
                selectInput("genesetView_plot_type", NULL,
                            c('boxplot', 'density'), 'boxplot')),
              column(4,
                selectInput("gensetView_plot_statistic", NULL,
                            c('logFC'='logFC', 't-statistic'='t'), 'logFC'))
            )
          ),

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
          uiOutput("gseaResultLabel"),
          uiOutput("gseaMethodSummary"),
          h4("Gene Set Statistics"),
          uiOutput("resultTableMessage"),
          DT::dataTableOutput("gseaResultTable")
        )
      ), ## tabItem: genesetView

      ## -----------------------------------------------------------------------
      ## tabName: geneView
      tabItem(
        tabName="geneView",
        fluidRow(
          box(
            p("Gene-centric exploration across genesets coming soon"),
            tags$ul(
              tags$li("Overlap of genes across genesets"),
              tags$li("etc ...")),
            plotlyOutput("geneView_volcanoPlot")),
          box(DT::dataTableOutput("geneView_interGenesetMembers"))
        )
      )

    )  ## dashboardBody::tabItems
  ) ## dashboardBody
)
