# The scripts were taken from the maser package and modified by Frank Ruehle (see # mod) and Anke Busch (see # modAB)
volcanoMod <- function (events, type = c("A3SS", "A5SS", "SE", "RI", "MXE"), 
          fdr = 0.05, deltaPSI = 0.1, title = "") 
{
  if (!is(events, "Maser")) {
    stop("Parameter events has to be a maser object.")
  }
  type <- match.arg(type)
  events <- as(events, "list")
  IncLevelDifference <- NULL
  Status <- NULL
  stats <- events[[paste0(type, "_", "stats")]]
  cond1 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference > deltaPSI)
  cond2 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference < (-1 * deltaPSI))
  status <- rep("Not significant", times = nrow(stats))
  #Ignored the up/down change made by FR, as we want to show the sample names directly
  #status[stats$ID %in% cond1$ID] <- "up"    # events$conditions[1] # mod
  #status[stats$ID %in% cond2$ID] <- "down"  # events$conditions[2] # mod
  status[stats$ID %in% cond1$ID] <- events$conditions[1]
  status[stats$ID %in% cond2$ID] <- events$conditions[2]
  FDR <- stats$FDR
  idx_zero <- which(stats$FDR == 0)
  idx_min_nonzero <- max(which(stats$FDR == 0)) + 1
  FDR[idx_zero] <- FDR[idx_min_nonzero]
  log10pval <- -1 * log10(FDR)
  #Ignored the up/down change made by FR, as we want to show the sample names directly
  #plot.df <- data.frame(ID = stats$ID, deltaPSI = stats$IncLevelDifference, 
  #                      log10pval = log10pval, Status = factor(status, levels = c("Not significant", "up", "down" # mod: up and down instead groups
  #                                                                                #events$conditions[1], events$conditions[2]
  #                                                                               )))
  plot.df <- data.frame(ID = stats$ID, deltaPSI = stats$IncLevelDifference, 
                        log10pval = log10pval, Status = factor(status, 
							       levels = c("Not significant", events$conditions[1], events$conditions[2])))
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
                        labels = c("Not significant",events$conditions[2],events$conditions[1])) +
    theme_bw() + 
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(face = "plain", colour = "black", size = 12), 
          axis.title.y = element_text(face = "plain", colour = "black", size = 12), 
          panel.grid.minor = element_blank(), 
          plot.background = element_blank()) + 
    labs(title = title, y = "Log10 Adj. Pvalue", x = "Delta PSI") # mod: switched axis labels
}
