#############################################################################
##
## What: GREAT.R
## Who: Giuseppe Petrosino, modified by Frank Rühle
## When: 2018-02-16
##
## Script to perform Genomic Regions Enrichment analysis
##
## Args: 
## -----
## peakData=                       # path to the xls file result from MACS2
############################################################################
library("ChIPpeakAnno")	
library("rGREAT")
library("Cairo")
library("stringr")


##
# get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
 
   if(length(i <- grep(string,args,fixed=T)) == 1)
    return(do.call(convert,list(gsub(string,"",args[i]))))
  if(!is.null(default)) default else do.call(convert,list(NA))
}


args             <- commandArgs(T)
peakData         <- parseArgs(args, "peakData=") # .xls result from MACS2
ftargets         <- parseArgs(args, "targets=","targets.txt")
out              <- parseArgs(args, "out=", "GREAT_analysis")
padj             <- parseArgs(args, "padj=", 0.01, "as.numeric") # padj parameter for filtering gene list
nterms           <- parseArgs(args, "nterms=", 5, "as.numeric") # number of terms to be shown in the barplot 
db               <- parseArgs(args, "db=", "hg19") # database 
adv_upstream     <- parseArgs(args, "adv_upstream=", 5, "as.numeric") # kb upstream of the TSS 
adv_downstream   <- parseArgs(args, "adv_downstream=", 1, "as.numeric") # kb downstream of the TSS 

runstr <- paste0("Call with: Rscript GREAT.R [peakData=",peakData,"] [targets=",ftargets,"] [out=",out,"] [padj=",padj,"] [nterms=",nterms,"] [db=",db,"] [adv_upstream=",adv_upstream,"] [adv_downstream=",adv_downstream,"]")
cat(runstr)
cat("\n\n")

## check if genome assembly is supported
supportedAssemblies <- c("hg38", "hg19", "mm9", "mm10") # supported in rGREAT version >=4
if (!db %in% supportedAssemblies) {
  print(paste("Selected genome assembly", db, "is not supported by GREAT v>4.0.0. Region enrichment is supported only for genome assemblies", paste(supportedAssemblies, collapse=", ")))
  peaks <- targets <- groups <- peaks.ov <- job <- go <- pathway <- NULL
  } else {


  # load targets
  targets <- read.table(ftargets,header=T)
  
  # define peak file names depending on whether blacklist was applied and/or input control used
  if(length(list.files(peakData,pattern="_macs2_blacklist_filtered_peaks.xls")) > 0) {
    peakFilenames <- gsub("\\.vs\\.none", "", paste0(peakData, "/", targets$IPname, ".vs.", targets$INPUTname, "_macs2_blacklist_filtered_peaks.xls")) # remove '.vs.none' if no input control used
  } else {
    peakFilenames <- gsub("\\.vs\\.none", "", paste0(peakData, "/", targets$IPname, ".vs.", targets$INPUTname, "_macs2_peaks.xls"))
  }
  
  # and return the tables
  peaks <- lapply(peakFilenames, function(x) {
    x <- tryCatch(read.delim(x, comment.char="#"), error=function(e) as.data.frame(matrix(ncol=10)))
    colnames(x) <- c("chr", "start", "end", "length", "summit", "tags", "-log10 pval", "fold enrichment", "-log10 FDR", "name")
    x <- x[grepl("^chr", x$chr),] # remove non-standard chromosomes which would otherwise crash submitGreatJob
    x[order(x$chr, x$start, x$end), c(-7, -10)]
  })
  
  names(peaks) <- gsub(" vs\\. none", "", paste0(targets$IPname, " vs. ", targets$INPUTname)) # remove ' vs. none' if no input control used
  

  
  groups <- unique(targets$group)
  
  peak.ranges <- lapply(peaks, function(x){
  		  x <- GRanges(seqnames=x$chr,
  			  IRanges(x$start,
  				  end=x$end),
  			  strand="*"
  			   )
  })
  peak.groups <- targets$group
  
  for(group in groups){
  		cat(paste0("#### ", group), fill=T)
  		cat("\n", fill=T)
  		peak <- peak.ranges[peak.groups==group]
  		peaks.ov <- findOverlapsOfPeaks(peak)
  		job	<- submitGreatJob(peaks.ov$peaklist[[length(peaks.ov$peaklist)]], species = db, adv_upstream = adv_upstream, adv_downstream = adv_downstream)
  
          # create genomic region and gene associations graphs 
          CairoPDF(file=paste0(out, "/", group, "_Region-Gene_Association_Graph.pdf"))
          res = plotRegionGeneAssociationGraphs(job)
          dev.off()
          
          availableOnts <- availableOntologies(job)   
          cat("\nOntologies available: ", paste(availableOnts, collapse=", "), "\n")
          
          # GO Biological Process 
          if("GO Biological Process" %in% availableOnts) {
            go <- as.data.frame(getEnrichmentTables(job, ontology = c("GO Biological Process")))
            go_filter <- go[go$GO.Biological.Process.Binom_Adjp_BH <= padj,]
            # Filtering step used by GREAT for avoiding general terms (Binom_Fold_Enrichment) and having terms statistically significant in both tests (see Suppl. Materials in McLean C. Nature Bio, 2010)
            go_filter <- go[go$GO.Biological.Process.Binom_Fold_Enrichment >= 2 & sort(go$GO.Biological.Process.Binom_Adjp_BH) & go$GO.Biological.Process.Hyper_Adjp_BH <= 0.05,]
            if(!is.null(go_filter) && nrow(go_filter) > 0) {
       
                # create barplot showing GO category
                CairoPNG(file=paste0(out, "/", group, "_GO_Barplot.png"), width=700, height=500)
                par(mar=c(5,15,4,5))
                print(barplot(-log10(as.numeric(rev(go_filter$GO.Biological.Process.Binom_Adjp_BH[1:nterms]))), names=str_wrap(rev(go_filter$GO.Biological.Process.name[1:nterms]), 
                width = 30), horiz=T , las=1, col="blue", xlab= "-log10(Binomal FDR Q-Val)", main= "GO Biological Process"))
                dev.off()
    
                write.csv(as.data.frame(go_filter),
                         file=paste0(out, "/", group, "_GO_Enrichment.csv"))
    
    	        }
          } else {go <- NULL}
          
          # MSigDB Pathway
          if("MSigDB Pathway" %in% availableOnts) {
            pathway <- as.data.frame(getEnrichmentTables(job, ontology = c("MSigDB Pathway")))
            pathway_filter <- pathway[pathway$MSigDB.Pathway.Binom_Adjp_BH <= padj,]
            # Filtering step used by GREAT for avoiding general terms (Binom_Fold_Enrichment) and having terms statistically significant in both tests (see Suppl. Materials in McLean C. Nature Bio, 2010)
            pathway_filter <- pathway_filter[pathway_filter$MSigDB.Pathway.Binom_Fold_Enrichment >= 2 & sort(pathway_filter$MSigDB.Pathway.Binom_Adjp_BH) & pathway_filter$MSigDB.Pathway.Hyper_Adjp_BH <= 0.05,]
            if(!is.null(pathway_filter) && nrow(pathway_filter) > 0) {
    
                # create barplot showing Pathway category
                CairoPNG(file=paste0(out, "/", group, "_MSigDB_Pathway_Barplot.png"), width=700, height=500)
                par(mar=c(5,15,4,5))
                print(barplot(-log10(as.numeric(rev(pathway_filter$MSigDB.Pathway.Binom_Adjp_BH[1:nterms]))), names=str_wrap(rev(pathway_filter$MSigDB.Pathway.name[1:nterms]), 
                width = 30), horiz=T , las=1, col="blue", xlab= "-log10(Binomal FDR Q-Val)", main= "MSigDB Pathway"))
                dev.off()
    
                write.csv(as.data.frame(pathway_filter),
                         file=paste0(out, "/", group, "_MSigDB_Pathway_Enrichment.csv"))
    
    	          }
          } else {pathway <- NULL}
  }
  # save the sessionInformation
  writeLines(capture.output(sessionInfo()), paste(out, "/GREAT_session_info.txt", sep=""))
}  

  save(peaks,targets,groups,peaks.ov,job,go,pathway,file=paste0(out,"/GREAT.RData"))

  