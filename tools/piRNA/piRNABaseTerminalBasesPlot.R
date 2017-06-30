#!/usr/bin/env Rscript

##############################################################################
### This script will plot the number of nucleotides around a piRNA read
###############################################################################


libs=c(
   "ggplot2",
   "RColorBrewer",
   "scales",
   "reshape2",
   "plyr",
   'grid',
   "gridExtra")

lapply(libs, require, character.only=T)
theme_set(theme_bw(16))
## create plot directory:
dir.create(file.path("figure"), showWarnings = FALSE)

readCleanNAs <- function(filename){
   counts <- read.csv(filename)
   counts[is.na(counts)] <- 0
   cols <- c('A', 'C', 'G', 'T', 'Position')
   counts <- counts[ ,colnames(counts) %in% cols]
   names(counts) <- gsub('T', 'U', names(counts))
   return(counts)
}


plotter <- function(df, title){
   p <- ggplot(df, aes(x=variable, y=value, fill=nucleotides)) +
      geom_bar(stat='identity', position = "fill") +
      scale_y_continuous(labels = percent_format()) +
      ggtitle(title) +
      ylab("Nucleotide count") +
      xlab("position relative to the piRNA sequence\n(+: inside piRNA)") +
      guides(fill=guide_legend(title=NULL)) +
      scale_fill_brewer(palette="Paired") +
      theme(axis.text.x=element_text(size=rel(0.8)), plot.title=element_text(vjust=1), axis.title.x=element_text(size=rel(0.8), vjust=-0.5), axis.title.y=element_text(size=rel(0.8), vjust=1.5), axis.text.y=element_text(size=rel(0.8)), legend.text=element_text(size=rel(0.8)))
   return(p)
}


# source('multiplot.R')

#'Plot multiple plots in a single pane
#'
#'ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
#' @import grid ggplot2
#' @export
#'
#' @param ... Two or more ggplot2 objects
#' @param  plotlist (optional) a list of ggplot2 objects
#' @param  cols Number of columns in layout
#' @param  layout A matrix specifying the layout. If present, 'cols' is ignored. See Details
#' @param  title Optional title as a character string
#' @param  widths a vector of relative column widths eg. c(3,2)
#' @param  heights a vector of relative column heights eg. c(3,2)
#' @param  titlefont The font of the title
#' @param  titleface The font face (1 = normal, 2 = bold, 3 = italic, 4 = bold italic)
#' @param  titlesize The size of the title font
#'
#' @details If plotting three plots and the layout is something like
#'   matrix(c(1,2,3,3), nrow=2, byrow=TRUE), then plot 1 will go in the upper
#'   left, 2 will go in the upper right, and 3 will go all the way across the
#'   bottom.  To save, you must use the desired device (eg \code{png()}), or
#'   save from the RStudio Viewer.
#'
#' Borrowed and modified from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#'Antonio: copied from https://github.com/ateucher/useful_code/blob/master/R/multiplot.r
#'
#' @return NULL (invisibly)
#' @examples \dontrun{
#' library("ggplot2")
#' plot1 <- ggplot(iris, aes(x = Species, y = Sepal.Length)) +
#'    geom_bar(stat = "identity")
#' plot2 <- ggplot(mtcars, aes(x = mpg, y = disp)) +
#'    geom_smooth()
#' multiplot(plot1, plot2, cols = 2, widths = c(3,2), title = "My two unrelated plots")
#' multiplot(plot1, plot2, cols = 1, heights = c(10,2), title = "My two unrelated plots")
#' myplots <- list(plot1, plot2, plot1)
#' multiplot(plotlist = myplots, layout =matrix(c(1,2,3,3), nrow=2),
#'      heights = c(1,3), widths = c(3,4), title = "My three unrelated plots")
#' ## Adjusting fonts
#' library(extrafont)
#' loadfonts()
#' multiplot(plotlist = myplots, layout =matrix(c(1,2,3,3), nrow=2),
#'           heights = c(1,3), widths = c(3,4), title = "My three unrelated plots",
#'           titlefont = "Wingdings", titleface = 4, titlesize = 20)
#'}
multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL, widths=NULL, heights=NULL,
                      title=NULL, titlefont = "", titleface = 1, titlesize = 16) {
  ggplot2:::theme_set(theme_bw(16))
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }

  if (!is.null(title)) { # Add a narrow row at the top for the title
    layout <- rbind(rep(0,ncol(layout)),layout)
    if (is.null(heights)) {
      plotrows <- nrow(layout)-1
      rowheights <- c(0.1, rep(1,plotrows)/plotrows)
    } else {
      rowheights <- c(0.1, heights/sum(heights))
    }
  } else {
    if (is.null(heights)) {
      rowheights <- rep(1,nrow(layout))
    } else {
      rowheights <- heights
    }
  }

  if (is.null(widths)) {
    colwidths <- rep(1, cols)
  } else {
    colwidths <- widths
  }

  if (numPlots==1) {

    return(plots[[1]] + labs(title=title))

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout),
                                               widths=colwidths,
                                               heights=rowheights)))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }

    if (!is.null(title)) {
      grid.text(title, vp = viewport(layout.pos.row = 1
                                     , layout.pos.col = 1:ncol(layout)),
                gp = gpar(fontfamily = titlefont, fontface = titleface,
                          fontsize = titlesize))
    }

  }
return(invisible(NULL))
}

# plotter(counts_m, exp)

counts_files <- list.files(pattern='.count.csv')
exps <- gsub('.count.csv', '', counts_files)

counts_l <- lapply(counts_files, readCleanNAs)
names(counts_l) <- exps

resize_fac <- nrow(counts_l[[1]])/10


counts_all <- list()
counts_all_p <- list()
for(exp in rev(names(counts_l))){
   print(exp)
   counts <- counts_l[[exp]]
   positions <- counts$Position
   counts$Position <- NULL
   counts <- data.frame(t(as.matrix(counts)))
   # len <- length(colnames(counts))/2
   nucleotides <- rownames(counts)
   counts <- data.frame(cbind(counts), nucleotides)

   colnames(counts) <- c(positions, 'nucleotides')
   counts_m <- melt(counts, id.vars='nucleotides')
   counts_m$nucleotides <- factor(counts_m$nucleotides,
   levels=c("A", "U", "C", "G"))

   ## size of the line to represent piRNA
   line_len <- length(levels(counts_m$variable))/2

   if(grepl('3prime', exp)){
      # positions <- rev(positions)
      colnames(counts) <- c(positions, 'nucleotides')
      counts_m <- melt(counts, id.vars='nucleotides')
      counts_m$nucleotides <- factor(counts_m$nucleotides,
   levels=c("A", "U", "C", "G"))
      counts_all[[exp]] <- counts_m
      counts_all_p[[exp]] <- plotter(counts_m, exp) +
      annotate("segment", x=0.5, xend=line_len + 0.5, y=1.01, yend=1.01, color="red", size = 1.5) +
      annotate("text", x=line_len - 1.5, y=1.035, label = "piRNA", color="red", size=3*resize_fac) +
      theme(text = element_text(size=9 * resize_fac))
   } else {
      colnames(counts) <- c(positions, 'nucleotides')
      counts_m <- melt(counts, id.vars='nucleotides')
      counts_m$nucleotides <- factor(counts_m$nucleotides,
   levels=c("A", "U", "C", "G"))
      counts_all[[exp]] <- counts_m
      counts_all_p[[exp]] <- plotter(counts_m, exp) +
      annotate("segment", x=line_len + 0.5, xend=2*line_len  + 0.5, y=1.01, yend=1.01, color="red", size=1.5) +
      annotate("text", x=line_len + 2.5, y=1.035, label = "piRNA", color="red", size=3*resize_fac) +
      theme(text = element_text(size=9 * resize_fac))
   }
}


pdf_file <- paste(
   'figure/',
   basename(getwd()),
   '.NucleotideDistributionOnPiRNA.pdf',
   sep='')

png_file <- paste(
   'figure/',
   basename(getwd()),
   '.NucleotideDistributionOnPiRNA.png',
   sep='')

pdf(pdf_file,
   width=12*resize_fac,
   height=8*resize_fac)
multiplot(plotlist=counts_all_p,layout=matrix(c(1,2,3,4), nrow=2, byrow=TRUE), title=basename(getwd()), titlesize = 12 * resize_fac)
# do.call("grid.arrange", c(counts_all_p, ncol=2, sub = textGrob("TITLE BELOW", gp=gpar(cex=2))))
dev.off()

png(png_file,
   width=720*resize_fac,
   height=640*resize_fac)
multiplot(plotlist=counts_all_p,layout=matrix(c(1,2,3,4), nrow=2, byrow=TRUE), title=basename(getwd()), titlesize = 12 * resize_fac)
# do.call("grid.arrange", c(counts_all_p, ncol=2))
dev.off()


# g1 <- do.call("arrangeGrob", c(counts_all_p, ncol=2)) +
#    ggtitle("testing 1.. 2")
# ggsave("figure/tes.pdf", g1)


# highlights <- data.frame(cyl=c(8))

# ggplot() +
#   geom_rect(aes(xmin=-Inf, xmax=3, ymin=-Inf, ymax=Inf), fill='red', alpha=0.2) +
#   geom_bar(data = mtcars, aes(x=factor(gear)), position="dodge", fill = 'black') +
#   facet_wrap(~cyl)

counts_all_p[[1]] +
  annotate("segment", x=0.5, xend=5.5, y=-1, yend=-1, color="red") +
  annotate("text", x = 4.5, y = max(counts_m$value)*0.1, label = "piRNA", color="red")


counts_all_p[[2]] +
  annotate("segment", x=5.5, xend=10.5, y=-1, yend=-1, color="red") +
  annotate("text", x = 6.5, y = max(counts_m$value)*0.1, label = "piRNA", color="red")
