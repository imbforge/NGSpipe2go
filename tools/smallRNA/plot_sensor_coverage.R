library("data.table")
library("ggplot2")
library("ggbio")
library("GenomicRanges")
library("scales")

# ==========================================================================
# Arguments
# ==========================================================================
# sensor_bed <- "/fsimb/groups/imb-kettinggr/genomes/Caenorhabditis_elegans/Ensembl/WBcel235/Annotation/Sensor/21U_sensor_to_plot.bed"
# normalization_factors_path <- "/fsimb/groups/imb-kettinggr/adomingues/projects/imb_ketting_2017_10_sensor/results/mapped/tracks/normalization_factors.txt"
# cov_files <- list.files(pattern='*.22G.minus.cov')
# base_name <- "all_sensor_strains"

args <- commandArgs(trailingOnly = TRUE)
sensor_bed <- args[1]
normalization_factors_path <- args[2]
cov_files <- args[3:length(args)]
print(args)

strand <- "minus"


# source("https://raw.githubusercontent.com/koundy/ggplot_theme_Publication/master/R/ggplot_theme_Publication.R")
theme_Publication <- function(base_size=14, base_family="Helvetica") {
      library(grid)
      library(ggthemes)
      (theme_foundation(base_size=base_size, base_family=base_family)
       + theme(plot.title = element_text(face = "bold",
                                         size = rel(1.2), hjust = 0.5),
               text = element_text(),
               panel.background = element_rect(colour = NA),
               plot.background = element_rect(colour = NA),
               panel.border = element_rect(colour = NA),
               axis.title = element_text(face = "bold",size = rel(1)),
               axis.title.y = element_text(angle=90,vjust =2),
               axis.title.x = element_text(vjust = -0.2),
               axis.text = element_text(), 
               axis.line.x = element_line(colour="black"),
               axis.line.y = element_line(colour="black"),
               axis.ticks = element_line(),
               panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               # legend.position = "bottom",
               # legend.direction = "vertical",
               legend.key.size= unit(0.2, "cm"),
               legend.spacing = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold")
       ))
      
}

scale_fill_Publication <- function(...){
      library(scales)
      discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
      
}

scale_colour_Publication <- function(...){
      library(scales)
      discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
      
}

publication_colors <- c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")

# dir.create(file.path("figure"), showWarnings = FALSE)

# ==========================================================================
# Read coverage data and parse file name
# ==========================================================================
read_and_summarize <- function(cov_file){
	print(cov_file)
	cov_tmp <- fread(cov_file)
  exp <- basename(cov_file)
	setnames(cov_tmp, c('Chromosome', 'Pos', 'Count'))
	cov_tmp[, Strain:=gsub('(\\w+_\\w+)_(\\w+)_(\\w+)\\.(\\w+)\\.(\\w+)\\.cov$', '\\1', exp), ]
	cov_tmp[, RNAclass:=gsub('(\\w+_\\w+)_(\\w+)_(\\w+)\\.(\\w+)\\.(\\w+)\\.cov$', '\\4', exp), ]
	cov_tmp[, Mapping:=gsub('(\\w+_\\w+)_(\\w+)_(\\w+)\\.(\\w+)\\.(\\w+)\\.cov$', '\\5', exp), ]
	cov_tmp[, Replicate:=gsub('(\\w+_\\w+)_(\\w+)_(\\w+)\\.(\\w+)\\.(\\w+)\\.cov$', '\\2', exp), ]
	cov_tmp[, Exp:=gsub('(\\w+_\\w+)_(\\w+)_(\\w+)\\.(\\w+)\\.(\\w+)\\.cov$', '\\1_\\2_\\3', exp), ]
	return(cov_tmp)
}

covs <- lapply(cov_files, read_and_summarize)
counts <- rbindlist(covs, fill=TRUE)
n_conds <- length(unique(counts$Strain))

# ==========================================================================
# Sensor coordinates
# ==========================================================================
sensor_coor <- fread(
      sensor_bed,
      col.names=c('chromosome', 'start', 'end', 'id', 'score', 'strand')
      )
sensor_coor$text_position <- (sensor_coor$start + sensor_coor$end)/2


genes <- makeGRangesFromDataFrame(
   sensor_coor,
   keep.extra.columns=TRUE,
   ignore.strand=FALSE,
   )


highl_lines <- c(
	sensor_coor[id=="mCherry"]$start,
	sensor_coor[id=="mCherry"]$end,
	sensor_coor[id=="21U"]$start,
	sensor_coor[id=="21U"]$end
	)


# ==========================================================================
# normalization
# ==========================================================================
fac <- fread(normalization_factors_path)
fac[, V1:=gsub(".All", "", V1)]


counts_rpm <- merge(counts, fac, by.x='Exp', by.y='V1')
setnames(counts_rpm, "V2", "libSize")
counts_rpm[, RPM:=Count/libSize*10^6, ]
counts_rpm[, Exp:=NULL,]

# ==========================================================================
# plot strand
# ==========================================================================
cov_minus  <- counts_rpm[Mapping == strand] 

p_cov <- ggplot(cov_minus, aes(x=Pos, y=RPM, colour=Replicate)) +
	geom_line(size = 1) +
	geom_vline(xintercept=highl_lines, linetype="dotted", colour="gray") +
	facet_wrap( ~ Strain, ncol=1) +
	theme_Publication() +
	theme(legend.position=c(0,1),legend.justification=c(0,1)) +
	scale_colour_manual(values = publication_colors[1:3]) +
	ylab('Reads per million \n(non-structural reads)') +
	xlab('Position along the 21ur loci')
# p_cov

# ==========================================================================
# plot sensor
# ==========================================================================

sensor_cols <-c( 
	"#ef3b2c",
	'#a6cee3',
	"#ef3b2c",
	'#a6cee3',
	'#a6cee3'
)

p_sensor <- ggplot(genes) +
	geom_rect(aes(fill = id)) +
	geom_text(aes(x=text_position, y=1, label=id)) + 
  	labs(x=NULL, y=NULL) +
	theme_Publication()  + 
	theme(legend.position="none") +
	scale_fill_manual(values=sensor_cols)

# p_sensor


# ==========================================================================
# assemble plots
# ==========================================================================

final <- tracks(
  p_cov, p_sensor,
  heights = c(n_conds*2, 0.5),
  # title=base_name,
  xlim=genes
  ) 
# print(final)
ggsave(paste0("coverage.replicates.pdf"), final, width = 10, height=n_conds*2.5)
ggsave(paste0("coverage.replicates.png"), final, width = 10, height=n_conds*2.5)


# ==========================================================================
# Standard deviation
# ==========================================================================

counts_cond <- counts_rpm[,list(Mean_reps=mean(RPM), SD=sd(RPM)), by=c('Pos', 'Strain', 'Mapping')][Mapping == strand]

p_cov_sd <- ggplot(counts_cond, aes(x=Pos, y=Mean_reps)) +
	geom_ribbon(aes(ymin=Mean_reps-SD, ymax=Mean_reps+SD), fill="black", alpha=0.2) +
	geom_line(colour="black", size = 1) +
	geom_vline(xintercept=highl_lines, linetype="dotted", colour="gray", size=1) +
	facet_wrap( ~ Strain, ncol=1) +
	theme_Publication() +
	theme(legend.position=c(0,1),legend.justification=c(0,1)) +
	scale_colour_Publication() +
	scale_fill_Publication() +
	ylab('Reads per million \n(non-structural reads)') +
	xlab('Position along the 21ur loci')


# ==========================================================================
# assemble plots
# ==========================================================================

final <- tracks(
  p_cov_sd, p_sensor,
  heights = c(n_conds*2, 0.5),
  # title=base_name,
  xlim=genes
  )
# print(final)
ggsave(paste0("coverage.SD.pdf"), final, width = 10, height=n_conds*2.5)
ggsave(paste0("coverage.SD.png"), final, width = 10, height=n_conds*2.5)
