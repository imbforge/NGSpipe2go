#####################################
##
## What: createMaserPlots.R
## Who : Sivarajan Karunanithi
## When: 27-09-2021
##
## This script takes the rMATS output and creates visualizations using the MASER package
##
## Args:
## -----
## gtf = gene_model.gtf       # gene model in gtf format
## db = ESSENTIAL_DB          # The genome DB used
## 
##
##
######################################

options(stringsAsFactors=FALSE)
library(maser)
library(rtracklayer)
library(knitr)
library(kableExtra)
library(openxlsx)
library(RColorBrewer)
##
## get arguments from the command line
##
parseArgs <- function(args,string,default=NULL,convert="as.character") {
	if(length(i <- grep(string,args,fixed=T)) == 1)
		return(do.call(convert,list(gsub(string,"",args[i]))))
	if(!is.null(default)) default else do.call(convert,list(NA))
}

args <- commandArgs(T)
gtf <- parseArgs(args,"gtf=","") # gtf gene model
db_gtf <- parseArgs(args, "db=","") # the essential_db parameter
ftype <- parseArgs(args, "ftype=", "") # which type of splicing events to consider: juncton counts (JC) or junction-exon counts (JCEC)
mincov <- parseArgs(args,"mincov=",5, convert="as.numeric") # ignore splicing events with read coverage below this count
fdr <- parseArgs(args,"fdr=","", convert="as.numeric") # FDR cut-off to select statistically significant splicing events identified by rMATS
dpsi <- parseArgs(args,"dpsi=","", convert="as.numeric") # minimum percentage spliced in (PSI) to include in plots
scripts_dir <- parseArgs(args,"scripts_dir=","") # Needed to load the modified maser scripts
rmats_dir <- parseArgs(args,"rmats_dir=","") # Location of the rMATS output folder for the current contrasts
group1 <- parseArgs(args,"group1=","") # Contrast groups
group2 <- parseArgs(args,"group2=","") # Contrast groups


##
## The maser package scripts had to be modified for minor bugs and
## to make it work for all organisms
##

source(file.path(scripts_dir,"volcanoMod.R"))

# print(c(group1,group2))
# load rMATS results - set a flag variable
error_status <- F

tryCatch( { 
	rmats <- maser(rmats_dir, cond_labels = c(group1,group2), ftype = ftype);
      }, 
      error = function(e) { 
	      print("Error occured while reading rMATS results (possibly no splicing events were found)");
	      error_status <<- T;
      }
)

if(!error_status) {

  # process gtf file
  ens_gtf <- rtracklayer::import.gff(gtf)
  GenomeInfoDb::genome(ens_gtf) <- db_gtf

  for(e in c("A3SS_events", "A5SS_events", "SE_events", "RI_events", "MXE_events")) { 
	# use gene IDs if no symbols available
	slot(rmats, e)$geneSymbol <- ifelse(is.na(slot(rmats, e)$geneSymbol), slot(rmats, e)$GeneID, slot(rmats, e)$geneSymbol)
  }

  # modAB: rMATS attaches "chr" to all chromosome names, which do not start with "chr" anyway.
  #        we could either do the same in the gtf annotation or remove the added "chr" from the
  #        rMATS results. let's try the second option here:
  rmats_mod <- rmats
  for(e in c("A3SS_gr", "A5SS_gr", "SE_gr", "RI_gr", "MXE_gr")) {

      # modify chromosome names in each GRanges object in e's slot
      slot(rmats_mod,e) <- GRangesList(lapply(slot(rmats,e),
                                              function(gr) {
                                                  seqlevels(gr)[!seqlevels(gr) %in% seqlevels(ens_gtf)] <- 
                                                          seqlevels(ens_gtf)[match(seqlevels(gr),paste0("chr",seqlevels(ens_gtf)),nomatch=0)]
                                                  return(gr)
                                              }))
  }
  # modAB: as a consequence, options(ucscChromosomeNames=FALSE) has to be set to allow for arbitrary chromosome identifiers
  #        when plotting (in plotTranscriptsMod.R). otherwise, functions createAnnotationTrack_event() -> 
  #        createAnnotationTrackA3SS_event() -> Gviz::AnnotationTrack() in plotTranscriptsMod.R will through an error

  # Filtering events: Low coverage splicing junctions are commonly found in RNA-seq data and lead to low confidence PSI levels.
  # We can remove low coverage events using filterByCoverage(), which may significantly reduced the number of splicing events.
  # modAB: use modified rmats object
  #rmats_filt <- maser::filterByCoverage(rmats, avg_reads = mincov)
  rmats_filt <- maser::filterByCoverage(rmats_mod, avg_reads = mincov)

  # The function topEvents() allows to select statistically significant events given a FDR cutoff and minimum PSI change.
  rmats_top <- maser::topEvents(rmats_filt, fdr = fdr, deltaPSI = dpsi)
  #print(rmats_top)

  pdf(paste0(rmats_dir,"/Global_SplicingEvents.pdf"))
  n_splice_events=0
  for(e in c("A3SS", "A5SS", "SE", "RI", "MXE")) {
	  n_splice_events = n_splice_events + nrow(slot(rmats_top, paste0(e,"_events")))
  }
  if(n_splice_events > 0) {
	  plotspldist <- maser::splicingDistribution(rmats_filt, fdr = fdr, deltaPSI = dpsi)
          print(plotspldist)
  }

  for(e in c("A3SS", "A5SS", "SE", "RI", "MXE")) {
	if(nrow(slot(rmats_top, paste0(e,"_events"))) > 0) {
		#plotPCA <- maser::pca(rmats_filt, type = e)
		#print(plotPCA)
		plotVolcano <- volcanoMod(rmats_filt, fdr = fdr, deltaPSI = dpsi, type = e, title = e)
		print(plotVolcano)
		top <- summary(rmats_top, type=e)
		#write the results for each event and contrast
		write.csv(top, file=paste0(rmats_dir,"/","Significant_",e,"_events.csv"), row.names = F)
		write.xlsx(top, file=paste0(rmats_dir,"/","Significant_",e,"_events.xlsx"), row.names = F, overwrite = T)	
	}
  }
  dev.off()

  save(rmats_filt, ens_gtf, file=paste0(rmats_dir,"/","maser_processed.RData"))
}
