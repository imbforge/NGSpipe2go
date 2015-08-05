library(shinyIncubator)	# the progress bar
library(markdown)		# display md tables (for the mapping and fastqc stats)
source("ChIP.shinyrep.helpers.R")	# helper functions to display the stats
loadGlobalVars()

shinyUI(fluidPage(
	progressInit(),
	titlePanel("ChIP peakcalling with MACS2"),
	sidebarLayout(
		sidebarPanel(
			# type of QC metric
			radioButtons("radio_qc","Initial QC metrics:",
						 c("SAV stats"="SAV",
						   "PCR bottleneck coefficient"="PBC",
						   "FastQC report"="FastQC",
						   "IPstrength"="IPstrength",
						   "PhantomPeak"="PhantomPeak")),
			tags$hr(),
			# type of mapping stats
			radioButtons("radio_mapping","Mapping:",
						 c("stats"="stats")),
			tags$hr(),
			# conditions compared
			uiOutput("radio_cont"),
			tags$hr(),
			# PDF-HTML-Word save knitr report: http://shiny.rstudio.com/gallery/download-knitr-reports.html
			radioButtons("radio_format","Document format",c("HTML","PDF","Word"),inline=TRUE),
			downloadButton("button_download")
		),
	    mainPanel(
	    	tabsetPanel(type="tabs",
	    		tabPanel("Initial QC metrics",uiOutput("qc")),
	    		tabPanel("Mapping",uiOutput("mapping")),
	    		tabPanel("Peaks",dataTableOutput("table"))
	    	),
	    	tags$hr(),
	    	textOutput("n_peaks")
		)
	)
))
