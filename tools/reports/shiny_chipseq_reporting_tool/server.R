library(rmarkdown)	# for the report generator
library(knitr)		# for the report generator
source("ChIP.shinyrep.helpers.R")		# helper functions to generate the plots

# things that will run only once per session
# always in the global environment to ensure access (could be placed in global.R also)
loadGlobalVars()

addResourcePath("fastqc"     ,SHINYREPS_FASTQC)	     # to access from the web the files stored in the QC dir
addResourcePath("ipstrength" ,SHINYREPS_IPSTRENGTH)	 # to access from the web the files stored in the QC dir
addResourcePath("phantompeak",SHINYREPS_PHANTOMPEAK) # to access from the web the files stored in the QC dir

##
## The server function
##
shinyServer(function(input,output,session) {
	
	##
	## init the interface
	##
	# create controls
	conts <- ChIPhelper.ComparisonsFromTargets()
	output$radio_cont <- renderUI({
 		radioButtons("radio_cont","Comparison:",conts)
	})
	radio_cont <<- 0
	
	# load data and format data
	withProgress(session,min=0,max=1, {
		setProgress(message="Initializing the interface...")
		setProgress(value=1,detail="Reading targets")
		targets <<- ChIPhelper.init("readTargets"); Sys.sleep(.1)
		setProgress(value=1,detail="Reading peaks")
		peaks   <<- ChIPhelper.init("readPeaks"); Sys.sleep(.1)
	})
	
	
	##
	## Render plots
	##
	#generate the QC HTML pages
	output$qc <- renderUI({
		switch(input$radio_qc,
			   SAV=HTML(renderMarkdown(text=ChIPhelper.Bustard())),
			   FastQC=HTML(renderMarkdown(text=ChIPhelper.Fastqc())),
			   IPstrength=HTML(renderMarkdown(text=ChIPhelper.IPstrength())),
			   PhantomPeak=HTML(renderMarkdown(text=ChIPhelper.PhantomPeak())),
			   PBC=HTML(renderMarkdown(text=ChIPhelper.PBC()))
			   )
	})

	#generate the mapping QC HTML pages
	output$mapping <- renderUI({
		switch(input$radio_mapping,
			   stats=HTML(renderMarkdown(text=ChIPhelper.Bowtie()))
			   )
	})
	
	##
	## peaks
	##
	Peaks <- reactive({
		input$radio_cont  	# make the reactive function dependent on the contrasts radio button
		
		# in case the radio button was not created yet in the client side
		if(radio_cont != 0) {	
			ChIPhelper.Peaks(i=which(conts == input$radio_cont))
		} else {
			radio_cont <<- 1
			ChIPhelper.Peaks(i=1)
		}
	})
	# render the table of DE genes
	output$table <- renderDataTable({
		Peaks()
	})
	# and display the number of DE genes
	output$n_peaks<- renderText({
		paste("Number of significant peaks:",nrow(Peaks()))
	})
	
	##
	## Downlaod report
	##
	output$button_download <- downloadHandler(
		filename <- function() {
			paste("ChIPreport",
				  switch(input$radio_format,HTML="html",PDF="pdf",Word="docx"),
				  sep=".")
		},
		content <- function(file) {

			# temporarily switch to the temp dir, in case you do not have write permission to the cwd
			# COMMENTED OUT: if we change the directory, we cannot access the included .md files anymore
#			src <- normalizePath("DEreport.Rmd")
#			owd <- setwd(tempdir())
#			on.exit(setwd(owd))
#			file.copy(src,"ChIPreport.Rmd")
			
			# render the file
			out <- render(paste0(SHINYREPS_PROJECT,"/reports/ChIPreport.Rmd"), {
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
