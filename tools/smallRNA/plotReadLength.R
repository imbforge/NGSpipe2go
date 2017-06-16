library('data.table')
library('ggplot2')
library('dplyr')
library('scales')
library('RColorBrewer')
library("scales")


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
               panel.grid.major = element_line(colour="#f0f0f0"),
               panel.grid.minor = element_blank(),
               legend.key = element_rect(colour = NA),
               legend.position = "bottom",
               legend.direction = "horizontal",
               legend.key.size= unit(0.2, "cm"),
               legend.spacing = unit(0, "cm"),
               legend.title = element_text(face="italic"),
               plot.margin=unit(c(10,5,5,5),"mm"),
               strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
               strip.text = element_text(face="bold"),
               plot.caption=element_text(size=12)
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


dir.create(file.path("figure"), showWarnings = FALSE)


readWithFileName <- function(file_path) {
  dt <- fread(file_path)
  dt$file <- gsub(
    "^([^.]*).*", "\\1",
    basename(file_path)
    )
  return(dt)
}

directory <- "./"
files <- list.files(directory, '.readlength.txt', recursive=TRUE)

counts <- lapply(files[!grepl('family|class', files)], readWithFileName)
counts <- rbindlist(counts, use.names=TRUE, fill=FALSE, idcol=NULL)

len <- counts %>%
   splitstackshape:::cSplit('file', '_', drop=TRUE, direction='wide')  %>%
   setnames(c('Count', 'Length', 'Genotype', 'DevStage', 'Replicate', 'Treatment')) %>%
   mutate(Sample=paste(Genotype, Replicate, sep="_")) %>%
   setDT
len

## generate colors for each replicate sample:
## http://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# cols3 <- c("#c75f65", "#949b48", "#9475c5")
# cols3 <- c("#aaa04a", "#8a5ab6", "#8d7779")
cols3 <- c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")[1:2]
cols <- c()
for (col in cols3){
colfunc <- colorRampPalette(c(col, "white"))
   cols <- c(
            cols,
            colfunc(15)[c(1,5,8)]
            )
}


ggplot(len, aes(x=factor(Length), y=Count, color=Sample)) +
   geom_line(aes(group=Sample)) + 
   geom_point() +
   scale_y_continuous(labels = comma) +
   scale_color_manual(values=cols) +
   facet_grid(Treatment ~ ., scales="free_y") +
   labs(
      x='Read length (bp)',
      y='Number of Reads',
      caption=paste('Source:', basename(getwd()))) +
   # xlab('Read length (bp)') +
   # ylab('Number of Reads') +
   theme_Publication() 
ggsave('figure/AllReadsLengthDistribution.pdf')
ggsave('figure/AllReadsLengthDistribution.png')

## calculate percentage

len <- len %>%
   group_by(Sample, Treatment) %>%
   mutate(Perc=Count/sum(Count)*100)

# len %>%
#    group_by(Sample, Treatment) %>%
#    summarize(Total=sum(Perc))

ggplot(len, aes(x=factor(Length), y=Perc, color=Sample)) +
   geom_line(aes(group=Sample)) + 
   geom_point() +
   scale_y_continuous(labels = comma) +
   scale_color_manual(values=cols) +
   facet_grid(Treatment ~ ., scales="free_y") +
   labs(
      x='Read length (bp)',
      y='% of Reads',
      caption=paste('Source:', basename(getwd()))) +
   theme_Publication()
ggsave('figure/PercentageReadsLengthDistribution.pdf')
ggsave('figure/PercentageReadsLengthDistribution.png')
