# The scripts were taken from the maser package and modified by Frank Ruehle (see # mod) and Anke Busch (see # modAB)
volcanoMod <- function (events, type = c("A3SS", "A5SS", "SE", "RI", "MXE"), 
          fdr = 0.05, deltaPSI = 0.1, title = "", top=25) 
{
  if (!is(events, "Maser")) {
    stop("Parameter events has to be a maser object.")
  }
  type <- match.arg(type)
  events <- as(events, "list")
  IncLevelDifference <- NULL
  Status <- NULL
  stats <- events[[paste0(type, "_", "stats")]]
  #we add the gene id and gene name to the stats to refer to it later
  event_info <- as.data.frame(events[[paste0(type, "_events")]]) 
  stats$geneSymbol <- event_info$geneSymbol[match(stats$ID, event_info$ID)] 

  cond1 <- dplyr::filter(stats, FDR < fdr & IncLevelDifference > deltaPSI)
  cond2 <- dplyr::filter(stats, FDR < fdr & IncLevelDifference < (-1 * deltaPSI))
  status <- rep("Not significant", times = nrow(stats))
  #Ignored the up/down change made by FR, as we want to show the sample names directly
  #status[stats$ID %in% cond1$ID] <- "up"    # events$conditions[1] # mod
  #status[stats$ID %in% cond2$ID] <- "down"  # events$conditions[2] # mod
  status[stats$ID %in% cond1$ID] <- events$conditions[1]
  status[stats$ID %in% cond2$ID] <- events$conditions[2]
  # modAB: stats is NOT sorted based on FDR column, i.e. extracting the row number
  #        after the last occurrence of FDR==0 and setting all FDR of 0 to the FDR seen
  #        in this row, will set them to a random value and NOT the next largest one after 0
  #        --> keep FDRs of 0, they will result in Inf values when -log10(FDR) and be plotted
  #            at the very top of the plot
  #FDR <- stats$FDR
  #idx_zero <- which(stats$FDR == 0)
  #idx_min_nonzero <- max(which(stats$FDR == 0)) + 1
  #FDR[idx_zero] <- FDR[idx_min_nonzero]
  log10pval <- -log10(stats$FDR)
  #Ignored the up/down change made by FR, as we want to show the sample names directly
  #plot.df <- data.frame(ID = stats$ID, deltaPSI = stats$IncLevelDifference, 
  #                      log10pval = log10pval, Status = factor(status, levels = c("Not significant", "up", "down" # mod: up and down instead groups
  #                                                                                #events$conditions[1], events$conditions[2]
  #                                                                               )))
  plot.df <- data.frame(ID = stats$ID, deltaPSI = stats$IncLevelDifference, 
                        log10pval = log10pval, Status = factor(status, 
							       levels = c("Not significant", events$conditions[1], events$conditions[2])),
                        geneSymbol=stats$geneSymbol)
  #order the results according to log10pval 
  plot.df <- plot.df[order(plot.df$log10pval, decreasing=T),]
  #amount of significant results
  sig_num <- sum(plot.df$Status!="Not significant")
  # select top sign. events
  sig.df <- plot.df[plot.df$Status!="Not significant",][1:min(top, sig_num),]

  # modAB: change maser color scheme:
  # modAB:   * always plot "Not significant" in grey
  # modAB:   * use the same group colors as for PCAs and heatmaps in DESeq2 section
  # modAB:     (the reference group (i.e. A in B.vs.A) is colored in brewer.pal(9,"Set1")[1])
  #if (length(unique(status)) < 3) {    # modAB
  #  colors <- c("blue", "red")         # modAB
  #}  else {                            # modAB
  #  colors <- c("grey", "blue", "red") # modAB
  #}                                    # modAB
  ggplot(plot.df, aes(x = deltaPSI, y = log10pval, colour = Status)) + 
    geom_point(aes(colour = Status),alpha=0.7) + # modAB: added transparency
    #scale_colour_manual(values = colors) + # modAB: use colors as given below
    scale_colour_manual(values = c("grey",brewer.pal(9,"Set1")[1:2]),
                        breaks = c("Not significant",events$conditions[2],events$conditions[1]),
                        labels = c("not significant",events$conditions[2],events$conditions[1])) +   
    geom_text_repel(data=sig.df,
                mapping=aes(deltaPSI, log10pval, label=geneSymbol),
                max.overlaps=50,
                show.legend=F) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(face = "plain", colour = "black", size = 12), 
          axis.title.y = element_text(face = "plain", colour = "black", size = 12), 
          plot.caption = element_text(hjust = 0.5),
          panel.grid.minor = element_blank(), 
          plot.background = element_blank()) + 
    labs(title = title, y = "-log10 adj. p-value",
         x = "delta PSI",
         caption=paste("The gene symbol is displayed for the top", min(top,sig_num), "events")
    ) # mod: switched axis labels
}
