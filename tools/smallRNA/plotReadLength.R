library('data.table')
library('ggplot2')
library('dplyr')
library("makeitprettier")
library('scales')
library('RColorBrewer')
library("scales")


readWithFileName <- function(file_path) {
  dt <- fread(file_path)
  dt$file <- gsub(
    "^([^.]*).*", "\\1",
    basename(file_path)
    )
  return(dt)
}

args <- commandArgs(trailingOnly = TRUE)
input_dir <- args[1]
output_dir <- args[1]

files <- list.files(
  input_dir,
  '.readlength.txt',
  recursive=TRUE,
  full.names=TRUE
  )
out_fig_path <- paste0(output_dir, "/plots")
dir.create(file.path(out_fig_path), showWarnings = FALSE)

counts <- lapply(files[!grepl('family|class', files)], readWithFileName)
counts <- rbindlist(counts, use.names=TRUE, fill=FALSE, idcol=NULL)

len <- counts %>%
   setDT() %>%
   setnames(c("V1", "V2"), c("Count", "Length")) %>%
   .[, Sample := gsub("\\..*", "", basename(as.character(file)))] %>%
   .[]
len


ggplot(len, aes(x=factor(Length), y=Count, color=Sample)) +
   geom_line(aes(group=Sample)) + 
   geom_point(size = 2) +
   scale_y_continuous(labels = comma) +
   scale_color_prettier() +
   labs(
      x='Read length (bp)',
      y='Number of Reads'
      ) +
   theme_redl(base_size = 18) 
  
save_plot(
  paste0(out_fig_path, '/AllReadsLengthDistribution'),
  save_data = TRUE,
  width = 10,
  height = 7
)

## calculate percentage

len <- len %>%
   group_by(Sample) %>%
   mutate(Perc=Count/sum(Count)*100)

ggplot(len, aes(x=factor(Length), y=Perc, color=Sample)) +
   geom_line(aes(group=Sample)) + 
   geom_point(size = 2) +
   scale_y_continuous(labels = comma) +
   scale_color_prettier() +
   labs(
      x='Read length (bp)',
      y='% of Reads in library'
      ) +
   theme_redl(base_size = 18) 
  
save_plot(
  paste0(out_fig_path, '/PercentageReadsLengthDistribution'),
  save_data = TRUE,
  width = 10,
  height = 7
)
