library(shinydashboard)
library(DT)
library(plotly)
library(ggplot2)
library(scater)

##
## Define UI for application
##
ui <- dashboardPage(
  dashboardHeader(title="Gut-Stomach scRNA-seq"),
  dashboardSidebar(disable=TRUE), 
  dashboardBody(
    fluidRow(
      column(width=9,
        box(width=NULL, solidHeader=TRUE, plotlyOutput("tsnePlot")),
        box(width=NULL, solidHeader=TRUE, plotlyOutput("expressionPlot"))
      ),
      column(width=3,
        box(width=NULL, status="warning",
            fluidRow(
              column(6, selectInput("dataset", label="dataset", choices=list())),
              column(6, numericInput("perp", label="perplexity", min=1, max=500, value=5))
            ),
            radioButtons("colorBy", "Color by", choices=c("cluster", "expression"), selected="cluster"),
            hr(),
            DT::dataTableOutput("genes"),
            hr(),
            fluidRow(
              column(6, downloadLink('downloadTsnePlot', 'Download tSNE plot')),
              column(6, downloadLink('downloadExpressionPlot', 'Download expression plot'))
            ),
            fluidRow(
              column(6, numericInput("width", label="width", min=1, max=20, value=7)),
              column(6, numericInput("height", label="height", min=1, max=20, value=7))
            )
        )
      )
    )
  )
)

##
## Define server logic
##
server <- function(input, output, session) {
 
  session$onSessionEnded(stopApp)
  pal <- rainbow(5)
  pal[5] <- "grey"
  
  # show the available datasets to the interface
  observe({
    updateSelectInput(session, "dataset", label="dataset", choices=list.files(pattern="\\.RData$"))
  })
  
  # load pre-rendered data
  data <- reactiveValues()
  observe({
    if(file.exists(input$dataset)) {
      withProgress(message="Reading data...", value=0, {
        x <- local({
          load(input$dataset)
          list(sce=sce, genes=genes, top.hvg=top.hvg)
        })
        data$sce <- x$sce
        data$genes <- x$genes
        data$top.hvg <- x$top.hvg
        rm(x)
      })
    }
  })
  
  # load the table with gene names
  output$genes <- DT::renderDataTable({
    if(!is.null(data$genes))
      DT::datatable(data$genes, selection="single", options=list(pageLength=5))
  })
  
  ##
  ## TSNE logic
  ##
  # the plot
  tsnePlot <- reactiveValues(p=NULL)
  observe({
    if(length(input$genes_rows_selected) == 1) {
      withProgress(message="Calculating tSNE plot", value=0, {
        p <- plotTSNE(data$sce,
                      colour_by=if(input$colorBy == "cluster" && "clusterDTree" %in% colnames(colData(data$sce))) "clusterDTree" else data$genes$gene_id[input$genes_rows_selected],
                      shape_by="sort",
                      run_args=list(exprs_values="norm_exprs",
                                      feature_set=data$top.hvg,
                                      perplexity=input$perp,
                                      rand_seed=100)) +
          ggtitle(data$genes$gene_name[input$genes_rows_selected]) +
          theme(legend.title=element_blank())
        if(input$colorBy == "expression") {
          p <- p + scale_color_gradientn(colors=rev(colorRampPalette(pal)(100)))
        }
        p$data$label <- rownames(p$data)
        p$mapping$label <- as.name("label")
        tsnePlot$p <- p
      })
    } else {
      tsnePlot$p <- NULL
    }
  })
  
  # display the plot
  output$tsnePlot <- renderPlotly({
    if(!is.null(tsnePlot$p))
      ggplotly(tsnePlot$p, tooltip=c("colour_by", "shape_by", "label"))
    else
      ggplotly(ggplot() + geom_blank() + theme_minimal())
  })
  
  # download the plot
  output$downloadTsnePlot <- downloadHandler(
    filename=function() {
      if(length(input$genes_rows_selected) == 1)
        paste0(data$genes$gene_name[input$genes_rows_selected], ".pdf")
      else
        "empty.pdf"
    },
    content=function(con) {
      pdf(con, width=input$width, height=input$height)
      if(!is.null(tsnePlot$p))
        p <- tsnePlot$p
      else
        p <- ggplot() + geom_blank() + theme_minimal()
      print(p)
      dev.off()
    }
  ) 

  ##
  ## expression plot logic
  ##
  # the plot
  expressionPlot <- reactiveValues(p=NULL)
  observe({
    if(length(input$genes_rows_selected) == 1) {
      p <- plotExpression(data$sce,
                          features=data$genes$gene_id[input$genes_rows_selected],
                          colour_by=if(input$colorBy == "cluster" && "clusterDTree" %in% colnames(colData(data$sce))) "clusterDTree" else NULL,
                          x=if("clusterDTree" %in% colnames(colData(data$sce))) "clusterDTree" else NULL,
                          exprs_values="norm_exprs") +
             ggtitle(data$genes$gene_name[input$genes_rows_selected]) +
             theme(legend.title=element_blank())
      p$data$label <- colnames(data$sce)
      p$mapping$label <- as.name("label")
      expressionPlot$p <- p
    } else {
      expressionPlot$p <- NULL
    }
  })
  
  # display the plot
  output$expressionPlot <- renderPlotly({
    if(!is.null(expressionPlot$p))
      ggplotly(expressionPlot$p, tooltip=c("label", "Y"))
    else
      ggplotly(ggplot() + geom_blank() + theme_minimal())
  })
  
  # download the plot
  output$downloadExpressionPlot <- downloadHandler(
    filename=function() {
      if(length(input$genes_rows_selected) == 1)
        paste0(data$genes$gene_name[input$genes_rows_selected], ".pdf")
      else
        "empty.pdf"
    },
    content=function(con) {
      pdf(con, width=input$width, height=input$height)
      if(!is.null(expressionPlot$p))
        p <- expressionPlot$p
      else
        p <- ggplot() + geom_blank() + theme_minimal()
      print(p)
      dev.off()
    }
  ) 
}

# Run the application 
shinyApp(ui = ui, server = server)
