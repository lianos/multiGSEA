library(miniUI)

shinyUI(fluidPage(
  title="multiGSEA Explorer",
  tags$head(
    tags$link(rel="stylesheet", type="text/css", href="dashboard.css"),
    tags$link(rel="stylesheet", type="text/css", href="miniUI.css")
  ),


  fluidRow(
    column(3, wellPanel(fileInput("mgresult", 'multiGSEA Result Upload'))),
    column(9, wellPanel(mgResultFilterUI("mg_result_filter")))),

  tabsetPanel(
    tabPanel(
      "Overview",
      fluidRow(
        column(
          12,
          tags$div(style="margin-bottom: 10px; padding: 5px; background-color: white",
                   title='GSEA Results',
                   uiOutput("gseaMethodSummary"))))),

    tabPanel(
      "GSEA Results",
      fluidRow(
        column(
          5, style="padding: 0", wellPanel(geneSetContrastViewUI("geneset_viewer"))),
        column(
          7, mgTableBrowserUI("mg_table_browser"))),
      fluidRow(
        column(
          12,
          tags$h4("Other Gene Sets with Selected Genes"),
          mgGeneSetSummaryByGeneUI('other_genesets_gsea')))
      ),

    tabPanel(
      "Differential Gene Expression",
      fluidRow(
        column(5, mgVolcanoUI("dge_volcano")),
        column(7,
               tags$div(
                 style="float:right",
                 downloadButton('download_dge_stats', 'Download Statistics')),
               DT::dataTableOutput("dge_volcano_genestats"))),
      fluidRow(
        column(
          12,
          tags$h4("Other Gene Sets with Selected Genes"),
          mgGeneSetSummaryByGeneUI('other_genesets_volcano'))))
  ) ## tabsetPanel

))
