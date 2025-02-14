
#
# rMATS plotting functions, adapted/extended from maser
# 

# splicingDistribution_mod: modified version of function maser::splicingDistribution 
# - maser::splicingDistribution() create one barplot summarizing the events with
#   higher inclusion in per group, events are given as percentage of all events in
#   all event types per group
# - splicingDistribution_mod() create the barplots are in maser::splicingDistribution()
#   and additionally create the same barplot, but instead of plotting percentages per
#   group, absolute number of events are plotted. in addition to these two barplots,
#   two additional ones are created. they show the number or percentage of events 
#   with higher in inclusion in one or the other group summarized per event type
splicingDistribution_mod <- function (events, fdr = 0.05, deltaPSI = 0.1) 
{
    if (!is(events, "Maser")) {
        stop("Parameter events has to be a maser object.")
    }
    events <- as(events, "list")
    as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
    nevents_cond1 <- rep(0, length(as_types))
    nevents_cond2 <- rep(0, length(as_types))

    for (i in 1:length(as_types)) {
        stats <- events[[paste0(as_types[i], "_", "stats")]]
        cond1 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference > deltaPSI)
        cond2 <- dplyr::filter(stats, FDR < fdr, IncLevelDifference < (-1 * deltaPSI))
        nevents_cond1[i] <- length(cond1$ID)
        nevents_cond2[i] <- length(cond2$ID)
    }
    nevents_prop1 <- if(sum(nevents_cond1,na.rm=T)==0) { nevents_cond1 } else { nevents_cond1/sum(nevents_cond1,na.rm=T) }
    nevents_prop2 <- if(sum(nevents_cond2,na.rm=T)==0) { nevents_cond2 } else { nevents_cond2/sum(nevents_cond2,na.rm=T) }

    nevents_type_prop1 <- ifelse(nevents_cond1+nevents_cond2==0,0,nevents_cond1/(nevents_cond1 + nevents_cond2))
    nevents_type_prop2 <- ifelse(nevents_cond1+nevents_cond2==0,0,nevents_cond2/(nevents_cond1 + nevents_cond2))

    condition <- c(rep(events$conditions[1], length(as_types)), 
                   rep(events$conditions[2], length(as_types)))
    condition <- factor(condition, 
                        levels = c(events$conditions[1], events$conditions[2]))

    df.plot <- data.frame(Condition  = condition, 
                          Type       = c(as_types, as_types), 
                          Count      = c(nevents_cond1, nevents_cond2),
                          Proportion = c(nevents_prop1, nevents_prop2),
                          Prop_type  = c(nevents_type_prop1, nevents_type_prop2))

    p.perc.cond <- ggplot(df.plot, aes(x = Condition, y = Proportion, colour = Type, fill = Type)) + 
        geom_bar(stat = "identity", alpha = 0.6) + 
        theme_bw() + 
        theme(legend.title = element_blank(), 
              panel.grid = element_blank()) + 
        scale_y_continuous(labels = scales::percent_format()) +
        labs(x = "",
             y = "% of splicing events") + 
        scale_fill_brewer(palette = "Set2") + 
        scale_color_brewer(palette = "Set2") + 
        coord_flip()

    p.count.cond <- ggplot(df.plot, aes(x = Condition, y = Count, colour = Type, fill = Type)) +
        geom_bar(stat = "identity", alpha = 0.6) +
        theme_bw() +
        theme(legend.title = element_blank(),
              panel.grid = element_blank()) + 
        labs(x = "",
             y = "# of splicing events") +
        scale_fill_brewer(palette = "Set2") + 
        scale_color_brewer(palette = "Set2") + 
        coord_flip()

    p.perc.type <- ggplot(df.plot, aes(x = Type, y = Prop_type, alpha = Condition, fill = Type, linetype = Condition)) +
        geom_bar(stat = "identity", color="gray50", width=0.8) +
        theme_bw() +
        theme(legend.position = "top",
              panel.grid = element_blank()) + 
        scale_y_continuous(labels = scales::percent_format()) +
        labs(x = "",
             y = "% of splicing events") +
        scale_fill_brewer(palette = "Set2") + 
        scale_alpha_manual(values=c(0.6,1)) + 
	scale_linetype_manual(values=c("dashed","solid")) +
        guides(alpha=guide_legend(nrow=1,reverse=T,title="higher inclusion in"),
	       linetype=guide_legend(nrow=1,reverse=T,title="higher inclusion in"),
               fill="none") + 
        coord_flip()
 
    p.count.type <- ggplot(df.plot, aes(x = Type, y = Count, alpha = Condition, fill = Type, linetype = Condition)) +
        geom_bar(stat = "identity", color="gray50", width=0.8) +
        theme_bw() +
        theme(legend.position = "top",
              panel.grid = element_blank()) + 
        labs(x = "",
             y = "# of splicing events") +
        scale_fill_brewer(palette = "Set2") + 
        scale_alpha_manual(values=c(0.6,1)) + 
	scale_linetype_manual(values=c("dashed","solid")) + 
        guides(alpha=guide_legend(nrow=1,reverse=T,title="higher inclusion in"),
	       linetype=guide_legend(nrow=1,reverse=T,title="higher inclusion in"),
               fill="none") +
        coord_flip()

    return(list(p.perc.cond  = p.perc.cond,
                p.count.cond = p.count.cond,
                p.perc.type  = p.perc.type,
                p.count.type = p.count.type))

}





