
#
# rMATS filter functions, adapted from maser
# 

# maser_mod: modified version of function maser::maser 
# - maser::maser() only reads inclusion counts
# - maser_mod() reads both inclusion and exclusion counts
maser_mod <- function (path, cond_labels, ftype = c("JCEC", "JC")) 
{
    ftype <- match.arg(ftype)
    rmats_out <- list.files(path, pattern = paste0(ftype, ".txt"), 
                            full.names = FALSE)
    mats <- new("Maser")
    counts.col <- c("IJC_SAMPLE_1", "SJC_SAMPLE_1","IJC_SAMPLE_2", "SJC_SAMPLE_2")
    
    if (!grepl("/$", path)) {
        path <- paste0(path, "/")
    }
    for (f in rmats_out) {
        events <- read.table(paste0(path, f), sep = "\t", stringsAsFactors = FALSE, header = TRUE)
        type <- unlist(strsplit(f, ".", fixed = TRUE))[1]
        inc1 <- strsplit(events[, counts.col[1]], ",")
        exc1 <- strsplit(events[, counts.col[2]], ",")
        inc2 <- strsplit(events[, counts.col[3]], ",")
        exc2 <- strsplit(events[, counts.col[4]], ",")
        reads.inc1 <- suppressWarnings(matrix(as.numeric(unlist(inc1)), 
                                              nrow = length(inc1), ncol = length(inc1[[1]]), byrow = TRUE))
        reads.inc2 <- suppressWarnings(matrix(as.numeric(unlist(inc2)), 
                                              nrow = length(inc2), ncol = length(inc2[[1]]), byrow = TRUE))
        reads.exc1 <- suppressWarnings(matrix(as.numeric(unlist(exc1)), 
                                              nrow = length(exc1), ncol = length(exc1[[1]]), byrow = TRUE))
        reads.exc2 <- suppressWarnings(matrix(as.numeric(unlist(exc2)), 
                                              nrow = length(exc2), ncol = length(exc2[[1]]), byrow = TRUE))
        reads.mat <- cbind(reads.inc1+reads.exc1, reads.inc2+reads.exc2)
        rownames(reads.mat) <- events$ID
        col_names <- c(paste0(cond_labels[1], "_", seq(1, length(inc1[[1]]), 
                                                       1)), paste0(cond_labels[2], "_", seq(1, length(inc2[[1]]), 
                                                                                            1)))
        colnames(reads.mat) <- col_names
        slot(mats, paste0(type, "_", "counts")) <- reads.mat
        inc1 <- strsplit(events[, "IncLevel1"], ",")
        inc2 <- strsplit(events[, "IncLevel2"], ",")
        reads.inc1 <- suppressWarnings(matrix(as.numeric(unlist(inc1)), 
                                              nrow = length(inc1), ncol = length(inc1[[1]]), byrow = TRUE))
        reads.inc2 <- suppressWarnings(matrix(as.numeric(unlist(inc2)), 
                                              nrow = length(inc2), ncol = length(inc2[[1]]), byrow = TRUE))
        reads.mat <- cbind(reads.inc1, reads.inc2)
        rownames(reads.mat) <- events$ID
        colnames(reads.mat) <- col_names
        slot(mats, paste0(type, "_", "PSI")) <- reads.mat
        mats@n_cond1 <- length(inc1[[1]])
        mats@n_cond2 <- length(inc2[[1]])
        mats@conditions <- cond_labels
        slot(mats, paste0(type, "_", "stats")) <- events[, c("ID", 
                                                             "PValue", "FDR", "IncLevelDifference")]
        grl <- maser:::create_GRanges(events, type)
        slot(mats, paste0(type, "_", "gr")) <- grl
        slot(mats, paste0(type, "_", "events")) <- events[, c("ID", 
                                                              "GeneID", "geneSymbol")]
    }
    return(mats)
}

# filterByCoverage_mod: modified version of function maser::filterByCoverage 
# - maser::filterByCoverage() keeps events with an average inclusion count
#                             (over all considered samples) > avg_reads
# - filterByCoverage_mod() keeps events with an average inclusion and average
#                          exclusion count >= avg_reads for both groups
filterByCoverage_mod <- function (events, cond_labels, avg_reads = 10) 
{
    ID <- NULL
    if (!is(events, "Maser")) {
        stop("Parameter events has to be a maser object.")
    }
    as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
    events <- as(events, "list")
    events_new <- list()
    lapply(as_types, function(type) {
        counts <- events[[paste0(type, "_", "counts"),drop=F]] # drop=F needed in case events[[paste0(type, "_", "counts")]] has only 1 row
        counts_1 <- counts[,unlist(lapply(colnames(counts),function(n) {s <- unlist(strsplit(n,"_")); paste(s[1:(length(s)-1)],collapse="_") }))==cond_labels[1]]
        counts_2 <- counts[,unlist(lapply(colnames(counts),function(n) {s <- unlist(strsplit(n,"_")); paste(s[1:(length(s)-1)],collapse="_") }))==cond_labels[2]]
        res_id_1 <- rownames(counts_1)[rowMeans(counts_1) >= avg_reads]
        res_id_2 <- rownames(counts_2)[rowMeans(counts_2) >= avg_reads]
        res_id <- intersect(res_id_1,res_id_2)
        slots <- maser:::filterByIds(type, events, res_id)
        lapply(names(slots), function(attrib) {
            events_new[[paste0(attrib)]] <<- slots[[paste0(attrib)]]
        })
    })
    events_new[["n_cond1"]] <- events$n_cond1
    events_new[["n_cond2"]] <- events$n_cond2
    events_new[["conditions"]] <- events$conditions
    return(as(events_new, "Maser"))
}




