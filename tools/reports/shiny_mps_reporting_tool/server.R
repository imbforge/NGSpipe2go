library(rmarkdown)	# for the report generator
library(knitr)		# for the report generator
source("mps.shinyrep.helpers.R")		# helper functions to generate the plots

# things that will run only once per session
# always in the global environment to ensure access (could be placed in global.R also)
loadGlobalVars()

addResourcePath("fastqc"     ,SHINYREPS_FASTQC_LOG)	     # to access from the web the files stored in the QC dir
addResourcePath("dupRadar"   ,SHINYREPS_DUPRADAR_LOG)    # to access from the web the files stored in the QC dir
#addResourcePath("RNAtypes"   ,SHINYREPS_RNATYPES_LOG)    # to access from the web the files stored in the QC dir
#addResourcePath("geneBodyCov",SHINYREPS_GENEBODYCOV_LOG) # to access from the web the files stored in the QC dir
#load(SHINYREPS_DE_EDGER,envir=.GlobalEnv)		# the outcome from the DE analysis

##
## The server function
##
shinyServer(function(input,output,session) {
	
	##
	## init the interface
	##
	# create controls
	output$radio_cont <- renderUI({
		radioButtons("radio_cont","Comparison:",colnames(conts))
	})
	radio_cont <<- 0

	# load data and format data
	withProgress(session,min=0,max=3, {
		setProgress(message="initializing the interface...")
		setProgress(value=1,detail="preparing DE genes table: formatting numbers")
		DEhelper.init("prepareDEdataTable"); Sys.sleep(.1)
		setProgress(value=2,detail="preparing DE genes table: building UCSC links")
		DEhelper.init("renderUcscGeneLinks"); Sys.sleep(.1)
		setProgress(value=3,detail="calculating sample-sample distances")
		DEhelper.init("prepareDistanceMatrix"); Sys.sleep(.1)
	})
	
	##
	## Render plots
	##
	#generate the QC HTML pages
	output$qc <- renderUI({
		switch(input$radio_qc,
			   SAV=HTML(renderMarkdown(text=DEhelper.Bustard())),
			   FastQC=HTML(renderMarkdown(text=DEhelper.Fastqc())),
			   dupRadar=HTML(renderMarkdown(text=DEhelper.dupRadar())),
			   RNAtypes=HTML(renderMarkdown(text=DEhelper.RNAtypes())),
			   geneBodyCov=HTML(renderMarkdown(text=DEhelper.geneBodyCov()))
			   )
	})

	#generate the mapping QC HTML pages
	output$mapping <- renderUI({
		switch(input$radio_mapping,
			   stats=HTML(renderMarkdown(text=DEhelper.STAR())),
			   parms=HTML(renderMarkdown(text=DEhelper.STARparms()))
			   )
	})

	#render plots using the helper functions
	output$plots <- renderPlot({
		switch(input$radio_plot,
			   MDS=DEhelper.MDS(),
			   cluster=DEhelper.cluster(),
			   corr=DEhelper.corr(),
			   var=DEhelper.var(),
			   MA=DEhelper.MAplot(i=which(colnames(conts) == input$radio_cont),fdr=input$slider_pval))
	})

	##
	## DE genes
	##
	DEgenes <- reactive({
		input$radio_cont  	# make the reactive function dependent on the contrasts radio button
		
		# in case the radio button was not created yet in the client side
		if(radio_cont != 0) {	
			DEhelper.DEgenes(i=which(colnames(conts) == input$radio_cont))
		} else {
			radio_cont <<- 1
			DEhelper.DEgenes(i=1)
		}
	})
	# render the table of DE genes
	output$table <- renderDataTable({
		DEgenes()
	})
	# and display the number of DE genes
	output$n_degenes <- renderText({
		paste("Number of DE genes:",sum(DEgenes()$FDR < input$slider_pval))
	})
	
	##
	## Downlaod report
	##
	output$button_download <- downloadHandler(
		filename <- function() {
			paste("mps.report",
				  switch(input$radio_format,HTML="html",PDF="pdf",Word="docx"),
				  sep=".")
		},
		content <- function(file) {

			# temporarily switch to the temp dir, in case you do not have write permission to the cwd
			# COMMENTED OUT: if we change the directory, we cannot access the included .md files anymore
#			src <- normalizePath("mps.report.Rmd")
#			owd <- setwd(tempdir())
#			on.exit(setwd(owd))
#			file.copy(src,"mps.report.Rmd")
			
			# render the file
			out <- render(paste0(SHINYREPS_PROJECT,"/reports/mps.report.Rmd"), {
					switch(input$radio_format,
						   HTML=html_document(),
						   PDF=pdf_document(),
						   Word=word_document())
					},
					output_file=basename(file),
					output_dir=dirname(file))
			
			file.rename(out,file)
		}
	)
})
