library(markdown)		# display md tables (for the mapping and fastqc stats)
source("mps.shinyrep.helpers.R")	# helper functions to display the stats
loadGlobalVars()

shinyUI(fluidPage(
	titlePanel("Multiplexed protein stability profiling"),
	sidebarLayout(
		sidebarPanel(
			# type of QC metric
			radioButtons("radio_qc","Initial QC metrics:",
						 c("SAV stats"="SAV",
						   "FastQC report"="FastQC",
						   "Duplication rate"="dupRadar",
						   "RNA types"="RNAtypes",
						   "Gene body coverage"="geneBodyCov")),
			tags$hr(),
			# type of mapping stats
			radioButtons("radio_mapping","Mapping:",
						 c("stats"="stats",
						   "STAR parms"="parms")),
			tags$hr(),
			# type of plot
			radioButtons("radio_plot","DE QC metrics:",
						c("MDS plot"="MDS",
						  "Cluster samples by top variant genes"="cluster",
						  "Sample-sample distance correlation"="corr",
						  "Variance estimation"="var",
						  "MA plot"="MA")),
			tags$hr(),
			# conditions compared
			uiOutput("radio_cont"),
			tags$hr(),
			# pval slider (FDR)
			sliderInput("slider_pval","FDR cutoff:",0,1,.05, step=.01),
			tags$hr(),
			# PDF-HTML-Word save knitr report: http://shiny.rstudio.com/gallery/download-knitr-reports.html
			radioButtons("radio_format","Document format",c("HTML","PDF","Word"),inline=TRUE),
			downloadButton("button_download")
		),
	    mainPanel(
	    	tabsetPanel(type="tabs",
	    		tabPanel("Initial QC metrics",uiOutput("qc")),
	    		tabPanel("Mapping",uiOutput("mapping")),
	    		tabPanel("DE QC metrics",plotOutput("plots")),
	    		tabPanel("DE genes",dataTableOutput("table"))

	    	),
	    	tags$hr(),
	    	textOutput("n_degenes")
		)
	)
))
