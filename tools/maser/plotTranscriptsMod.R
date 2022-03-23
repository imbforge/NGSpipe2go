# The scripts were taken from the maser package and modified by Frank Ruehle (see # mod)
plotTranscriptsMod <- function (events, type = c("A3SS", "A5SS", "SE", "RI", "MXE"), 
          event_id, gtf, zoom = FALSE, show_PSI = TRUE, title="")  # modFR add title
{
  is_strict = TRUE
  if (!is(events, "Maser")) {
    stop("Parameter events has to be a maser object.")
  }
  if (!class(gtf) == "GRanges") {
    stop(cat("\"gtf\" should be a GRanges class."))
  }
  if (any(!grepl("chr", GenomeInfoDb::seqlevels(gtf)))) {
    GenomeInfoDb::seqlevels(gtf) <- paste0("chr", GenomeInfoDb::seqlevels(gtf))
  }
  # std_chr <- c(paste0("chr", seq(1:22)), "chrX", "chrY") # modFR outcomment chromosome filtering
  # if (any(!seqlevels(gtf) %in% std_chr)) {
  #   GenomeInfoDb::seqlevels(gtf, pruning.mode = "coarse") <- std_chr
  # }
  type <- match.arg(type)
  events <- as(events, "list")
  annot <- events[[paste0(type, "_", "events")]]
  if (length(unique(annot$geneSymbol)) > 1) {
    stop(cat("Multiple genes found. Use geneEvents() to select \n             gene-specific AS events."))
  }
  grl <- events[[paste0(type, "_", "gr")]]
  idx.event <- grep(as.numeric(event_id), grl[[1]]$ID)
  eventGr <- lapply(names(grl), function(exon) {
    return(grl[[exon]][idx.event])
  })
  eventGr <- GRangesList(eventGr)
  names(eventGr) <- names(grl)
  eventTrack <- createAnnotationTrack_event(eventGr, type)
  gtf_exons <- gtf[gtf$type == "exon", ]
  txnTracks <- createAnnotationTrack_transcripts(eventGr, gtf_exons, 
                                                 type, is_strict)
  if (show_PSI) {
    PSI <- events[[paste0(type, "_", "PSI")]]
    PSI_event <- PSI[idx.event, , drop = FALSE]
    groups <- factor(c(rep(events$conditions[1], events$n_cond1), 
                       rep(events$conditions[2], events$n_cond2)), levels = events$conditions)
    psiTrack <- createPSITrack_event(eventGr, PSI_event, 
                                     groups, type, zoom)
    trackList <- list(psiTrack, eventTrack, txnTracks$inclusionTrack, 
                      txnTracks$skippingTrack)
  } else {
    trackList <- list(eventTrack, txnTracks$inclusionTrack, 
                      txnTracks$skippingTrack)
  }
  if (zoom) {
    Gviz::plotTracks(trackList, col.line = NULL, col = NULL, main = title, # modFR add title
                     Inclusion = "orange", Skipping = "purple", Retention = "orange", 
                     Non_Retention = "purple", MXE_Exon1 = "orange", MXE_Exon2 = "purple", 
                     A5SS_Short = "orange", A5SS_Long = "purple", A3SS_Short = "orange", 
                     A3SS_Long = "purple", from = start(range(unlist(eventGr))) - 
                       500, to = end(range(unlist(eventGr))) + 500)
  }  else {
    Gviz::plotTracks(trackList, col.line = NULL, col = NULL,  main = title, # modFR add title
                     Inclusion = "orange", Skipping = "purple", Retention = "orange", 
                     Non_Retention = "purple", MXE_Exon1 = "orange", MXE_Exon2 = "purple", 
                     A5SS_Short = "orange", A5SS_Long = "purple", A3SS_Short = "orange", 
                     A3SS_Long = "purple")
  }
}



createAnnotationTrack_event <- function(eventGr, type){
  
  if(type == "A3SS") {
    event_track <- createAnnotationTrackA3SS_event(eventGr)
  }
  
  if(type == "A5SS") {
    event_track <- createAnnotationTrackA5SS_event(eventGr)
  }
  
  if (type == "SE"){
    event_track <- createAnnotationTrackSE_event(eventGr)
  }
  
  if (type == "RI"){
    event_track <- createAnnotationTrackRI_event(eventGr)
  }
  
  if (type == "MXE"){
    event_track <- createAnnotationTrackMXE_event(eventGr)
  }
  
  return(event_track)
  
}

createAnnotationTrack_transcripts <- function(eventGr, gtf_exons, 
                                              type, is_strict){
  
  if(type == "A3SS") {
    txn_tracks <- createAnnotationTrackA3SS_transcripts(eventGr, gtf_exons)
  }
  
  if(type == "A5SS") {
    txn_tracks <- createAnnotationTrackA5SS_transcripts(eventGr, gtf_exons)
  }
  
  if (type == "SE"){
    txn_tracks <- createAnnotationTrackSE_transcripts(eventGr, gtf_exons,
                                                      is_strict)
  }
  
  if (type == "RI"){
    txn_tracks <- createAnnotationTrackRI_transcripts(eventGr, gtf_exons,
                                                      is_strict)
  }
  
  if (type == "MXE"){
    txn_tracks <- createAnnotationTrackMXE_transcripts(eventGr, gtf_exons,
                                                       is_strict)
  }
  
  return(txn_tracks)
  
}
#' @importFrom dplyr filter
createExonTable <- function(gtf_exons, ids){
  
  transcript_id <- NULL
  # Recover exons of transcripts for the inclusion track using transcript
  #IDs
  res <- dplyr::filter(as.data.frame(gtf_exons), 
                       transcript_id %in% ids)
  
  # Create data frame for transcript track - follow the model 
  #from data(geneModels)
  res.df <- res[, c("seqnames", "start", "end", "strand", "exon_id",
                    "transcript_id")]  # modFR
                #  "transcript_name")] # modFR
  colnames(res.df) <- c("chromosome","start","end","strand","exon", 
                        "transcript")
  
  return(res.df)
}

#' @importFrom Gviz GeneRegionTrack
#' @importFrom GenomicRanges GRanges
createTxnTrack <- function(res.df, trackLabel, featureName){
  
  if (nrow(res.df) > 0){ 
    res.df$feature <- featureName
    txnTrack <- Gviz::GeneRegionTrack(range = res.df, 
                                      name = trackLabel, 
                                      transcriptAnnotation = "transcript",
                                      col = NULL,
                                      col.line = NULL)  
  }else {
    txnTrack <- Gviz::GeneRegionTrack(range = GRanges(),
                                      name = trackLabel, 
                                      transcriptAnnotation = "transcript",
                                      col = NULL,
                                      col.line = NULL)  
  }
  return(txnTrack)
  
}

#' @importFrom Gviz GeneRegionTrack
createAnnotationTrackSE_transcripts <- function(eventGr, gtf_exons,
                                                is_strict){
  
  tx_ids <- mapTranscriptsSEevent(eventGr, gtf_exons, is_strict)
  
  # Inclusion track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_3exons)
  inclusionTrack <- createTxnTrack(res.df, "Inclusion", "Inclusion")
  
  # Skipping track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_2exons)
  skippingTrack <- createTxnTrack(res.df, "Skipping", "Skipping")
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackRI_transcripts <- function(eventGr, gtf_exons,
                                                is_strict){
  
  tx_ids <- mapTranscriptsRIevent(eventGr, gtf_exons, is_strict)
  
  # Retention track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_retention)
  retention_Track <- createTxnTrack(res.df, "Retention", "Retention")
  
  # Non-retention track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_nonRetention)
  nonRetention_Track <- createTxnTrack(res.df, "Non-retention", "Non_Retention")
  
  txn_tracks <- list("inclusionTrack" = retention_Track,
                     "skippingTrack" = nonRetention_Track)
  return(txn_tracks)
  
}

createAnnotationTrackMXE_transcripts <- function(eventGr, gtf_exons,
                                                 is_strict){
  
  tx_ids <- mapTranscriptsMXEevent(eventGr, gtf_exons, is_strict)
  
  # MXE Exon 1 track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_mxe_exon1)
  inclusionTrack <- createTxnTrack(res.df, "MXE Exon 1", "MXE_Exon1")
  
  # MXE Exon 2 track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_mxe_exon2)
  skippingTrack <- createTxnTrack(res.df, "MXE Exon 2", "MXE_Exon2")
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackA5SS_transcripts <- function(eventGr, gtf_exons){
  
  #reverse strand becomes A3SS
  if (as.character(strand(eventGr[1])) == "+"){
    tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
  }else{
    tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
  }
  
  # Short exon track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_short)
  inclusionTrack <- createTxnTrack(res.df, "A5SS Short", "A5SS_Short")
  
  # Long track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_long)
  skippingTrack <- createTxnTrack(res.df, "A5SS Long", "A5SS_Long")
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

createAnnotationTrackA3SS_transcripts <- function(eventGr, gtf_exons){
  
  #reverse strand becomes A5SS
  if (as.character(strand(eventGr[1])) == "+"){
    tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
  }else{
    tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
  }
  
  # Short exon track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_short)
  inclusionTrack <- createTxnTrack(res.df, "A3SS Short", "A3SS_Short")
  
  # Long exon track
  res.df <- createExonTable(gtf_exons, tx_ids$txn_long)
  skippingTrack <- createTxnTrack(res.df, "A3SS Long", "A3SS_Long")
  
  txn_tracks <- list("inclusionTrack" = inclusionTrack,
                     "skippingTrack" = skippingTrack)
  return(txn_tracks)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
createAnnotationTrackSE_event <- function(eventGr){
  
  transcript_id <- NULL
  trackGr <- c(unlist(eventGr), unlist(eventGr[2:3]))
  trackGr$group <- rep(c("Inclusion", "Skipping"), c(3, 2))
  trackGr$type <- rep("Exon skipping", 5)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                       groupAnnotation = "group", 
                                       shape = "box",
                                       stacking = "squish", 
                                       id = "Exon skipping",
                                       col = NULL,
                                       col.line = NULL)
  Gviz::feature(event_track) <- rep(c("Inclusion", "Skipping"), c(3, 2))
  
  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
createAnnotationTrackRI_event <- function(eventGr){
  
  transcript_id <- NULL
  trackGr <- c(eventGr$exon_ir, eventGr$exon_upstream, 
               eventGr$exon_downstream)
  trackGr$group <- rep(c("Retention", "Non-retention"), c(1, 2))
  trackGr$type <- rep("Intron retention", 3)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                       groupAnnotation = "group", 
                                       shape = "box",
                                       stacking = "squish", 
                                       id = "Intron retention",
                                       col = NULL,
                                       col.line = NULL)
  
  Gviz::feature(event_track) <- rep(c("Retention", "Non_Retention"),
                                    c(1, 2))
  
  
  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
#' 
createAnnotationTrackMXE_event <- function(eventGr){
  
  trackGr <- c(eventGr$exon_upstream, eventGr$exon_1, 
               eventGr$exon_downstream, eventGr$exon_upstream,
               eventGr$exon_2, eventGr$exon_downstream)
  
  trackGr$group <- rep(c("MXE_Exon1", "MXE_Exon2"), c(3, 3))
  trackGr$type <- rep("Mutually Exclusive Exons", 3)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                       groupAnnotation = "group", shape = "box",
                                       stacking = "squish", 
                                       id = "Mutually Exclusive Exons",
                                       col = NULL,
                                       col.line = NULL)
  
  Gviz::feature(event_track) <- rep(c("MXE_Exon1", "MXE_Exon2"), c(3, 3))
  
  
  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
createAnnotationTrackA5SS_event <- function(eventGr){
  
  trackGr <- c(eventGr$exon_short, eventGr$exon_flanking,
               eventGr$exon_long, eventGr$exon_flanking)
  trackGr$group <- rep(c("A5SS Short", "A5SS Long"), c(2, 2))
  trackGr$type <- rep("A5SS", 4)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                       groupAnnotation = "group", 
                                       shape = "box",
                                       stacking = "squish", id = "A5SS",
                                       col = NULL,
                                       col.line = NULL)
  Gviz::feature(event_track) <- rep(c("A5SS_Short", "A5SS_Long"), c(2, 2))
  
  
  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom Gviz feature
createAnnotationTrackA3SS_event <- function(eventGr){
  
  trackGr <- c(eventGr$exon_flanking, eventGr$exon_short,
               eventGr$exon_flanking, eventGr$exon_long)
  trackGr$group <- rep(c("A3SS Short", "A3SS Long"), c(2, 2))
  trackGr$type <- rep("A3SS", 4)
  
  event_track <- Gviz::AnnotationTrack(trackGr, name = "Event", 
                                       groupAnnotation = "group", 
                                       shape = "box",
                                       stacking = "squish", 
                                       id = "A3SS",
                                       col = NULL,
                                       col.line = NULL)
  Gviz::feature(event_track) <- rep(c("A3SS_Short", "A3SS_Long"), c(2, 2))
  
  
  return(event_track)
  
}

#' @importFrom Gviz AnnotationTrack
#' @importFrom GenomicRanges GRanges
#' @importFrom parallel mclapply
createUniprotKBtracks <- function(eventGr, features, protein_ids, ncores = 1){
  
  options(ucscChromosomeNames=FALSE)
  
  if(.Platform$OS.type == "windows"){
    ncores = 1
  }
  
  uniprotTracks <- mclapply(seq_along(features), function(i){
    
    feature_gr <- createGRangesUniprotKBtrack(features[i])
    ovl_gr <- overlappingFeatures(feature_gr, eventGr)
    #ovl_gr_filt <- ovl_gr[ovl_gr$Uniprot_ID %in% protein_ids, ] 
    ovl_gr_filt <- ovl_gr
    
    uniq_features <- match(unique(as.vector(ovl_gr_filt$Name)), 
                           as.vector(ovl_gr_filt$Name))
    
    ovl_gr_filt_uniq <- ovl_gr_filt[uniq_features, ]
    
    if (length(ovl_gr_filt_uniq) > 0 ){
      
      track <- Gviz::AnnotationTrack(range = ovl_gr_filt_uniq, 
                                     name = features[i], 
                                     id = ovl_gr_filt_uniq$Name, 
                                     showFeatureId = TRUE,
                                     fontcolor.feature = "darkblue",
                                     #fill = "aliceblue", 
                                     shape = "box",
                                     col = NULL,
                                     col.line = NULL)  
      
    }else {
      track <- Gviz::AnnotationTrack(range = GRanges(), name = features[i])
    }
    
    return(track)
    
  }, mc.cores = ncores)
  
  return(uniprotTracks)  
  
}



createPSITrack_event <- function(eventGr, PSI_event, groups, type, zoom){
  
  if(type == "A3SS" || type == "A5SS" ) {
    psi_track <- createPSIDataTrack(eventGr, PSI_event, groups, zoom, 
                                    eventGr$exon_long)
  }
  
  if (type == "SE"){
    psi_track <- createPSIDataTrack(eventGr, PSI_event, groups, zoom, 
                                    eventGr$exon_target)
  }
  
  if (type == "RI"){
    psi_track <- createPSIDataTrack(eventGr, PSI_event, groups, zoom, 
                                    eventGr$exon_ir)
  }
  
  if (type == "MXE"){
    psi_track <- createPSITrackMXE_event(eventGr, PSI_event, groups, zoom)
  }
  
  return(psi_track)
  
}

#' @importFrom Gviz DataTrack
#' @importFrom GenomicRanges ranges
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges values
createPSIDataTrack <- function(eventGr, PSI_event, groups, zoom, exonGr){
  
  if(zoom){
    trackGr <- exonGr  
  }else{
    #create space for boxplot plotting
    trackGr <- range(unlist(eventGr))
    start(trackGr) <- start(trackGr) - 200
    end(trackGr) <- end(trackGr) + 200
  }
  
  values(trackGr) <- PSI_event 
  psi_track <- Gviz::DataTrack(trackGr, 
                               name = "PSI", 
                               groups = groups,
                               feature = groups,
                               legend = TRUE,
                               type = c("boxplot"),
                               fill = c("#0000ff", "#ff0000"),
                               col = c("#0000ff", "#ff0000"),
                               col.line = c("#0000ff", "#ff0000"))
  
  return(psi_track)
  
}


#' @importFrom Gviz DataTrack
#' @importFrom GenomicRanges ranges
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges values
createPSITrackMXE_event <- function(eventGr, PSI_event, groups, zoom){
  
  if(zoom){
    trackGr <- c(eventGr$exon_1, eventGr$exon_2)
    start(trackGr) <- start(trackGr) - 50
    end(trackGr) <- end(trackGr) + 50
  }else{
    #create space for boxplot plotting
    trackGr <- c( range(c(range(eventGr$exon_upstream), 
                          range(eventGr$exon_1))),
                  
                  range(c(range(eventGr$exon_2), 
                          range(eventGr$exon_downstream)))
    )
    start(trackGr) <- start(trackGr) - 200
    end(trackGr) <- end(trackGr) + 200
  }
  #rMATS PSI levels are for Exon 1
  data <- rbind(PSI_event, 1-PSI_event )
  values(trackGr) <- data
  psi_track <- Gviz::DataTrack(trackGr, 
                               name = "PSI", 
                               groups = groups,
                               legend = TRUE,
                               type = c("boxplot"),
                               fill = c("blue", "red"),
                               col = c("blue", "red"))
  return(psi_track)
  
  
}




#' Mapping of splice events to UniprotKB protein features.
#' 
#' @param events a maser object with transcript and protein identifiers.
#' @param tracks a character vector indicating valid UniprotKB features or 
#' categories.
#' @param by a character vector, possible values 
#' are \code{c("feature", "category")}.
#' @param ncores number of cores for multithreading (available only in OSX and Linux 
#' machines). If Windows, \code{ncores} will be set to 1 automatically.
#' @return a maser object with protein feature annotation.
#' @details This function performs mapping of splicing events to protein
#'  features available in the UniprotKB database. Annotation tracks of protein
#'  features mapped to the hg38 build of the human genome are retrieved from the
#'  public UniprotKB FTP. The function will overlap exons involved in the splice
#'  event with the feature genomic coordinates retrieved from UniprotKB.
#' 
#' Annotation can be executed either by feature or category. If categories are 
#' provided, all features within the category group will be included for 
#' annotation.
#' 
#' Thus, batch annotation is enabled either by using \code{by = category} or 
#' by providing mutilple features in the \code{tracks} argument.
#' 
#' Visualization of protein features can be done 
#' using \code{\link{plotUniprotKBFeatures}}.
#'  
#' @examples
#' ## Create the maser object
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' 
#' ## Ensembl GTF annotation for SRSF6
#' gtf_path <- system.file("extdata", file.path("GTF", "SRSF6_Ensembl85.gtf"),
#'  package = "maser")
#' ens_gtf <- rtracklayer::import.gff(gtf_path)
#' 
#' ## Retrieve gene specific splice events
#' srsf6_events <- geneEvents(hypoxia_filt, geneS = "SRSF6")
#' 
#' ## Map splicing events to transcripts
#' srsf6_mapped <- mapTranscriptsToEvents(srsf6_events, ens_gtf)
#' 
#' ## Annotate splice events with protein domains
#' srsf6_annot <- mapProteinFeaturesToEvents(srsf6_mapped, tracks = "domain")
#' head(annotation(srsf6_annot, "SE"))
#' 
#' @seealso \code{\link{plotUniprotKBFeatures}}
#' @export
#' @import GenomicRanges
#' @importFrom dplyr filter
#' @importFrom parallel mclapply
#' @import methods

mapProteinFeaturesToEvents <- function(events, tracks, by = c("feature", 
                                                              "category"), ncores = 1){
  
  by <- match.arg(by)
  
  if(!is(events, "Maser")){
    stop("Parameter events has to be a maser object.")
  }
  
  if(.Platform$OS.type == "windows"){
    ncores = 1
  }
  
  df <- availableFeaturesUniprotKB()
  
  if(by == "feature"){
    
    if (!any(tracks %in% as.vector(df$Name))){
      stop(cat("\"tracks\" arg is invalid."))
    }
    features <- tracks
  }
  
  Category <- NULL
  
  if(by == "category"){
    
    if (!any(tracks %in% as.vector(df$Category))){
      stop(cat("\"tracks\" arg is invalid."))
    }
    
    df_filt <- dplyr::filter(df, Category %in% tracks)
    features <- as.vector(df_filt$Name)
  }
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  events <- as(events, "list")
  
  # Add UniprotKB features annotation
  events_with_features <- events
  
  # Create GRanges for each UniprotKB feature in features argument
  features_Gr <- mclapply(features, createGRangesUniprotKBtrack, 
                          mc.cores = ncores)
  features_Gr <- GRangesList(features_Gr)
  names(features_Gr) <- features
  
  for (type in as_types){
    
    # Retrieve Events annotation
    annot <- events[[paste0(type,"_","events")]]
    grl <- events[[paste0(type,"_","gr")]]
    idx.cols <- grep("^list_ptn_", colnames(annot))
    
    if (nrow(annot) == 0){
      next
    }
    
    #Create matrix for storing results
    annot_uniprotkb <- matrix("NA", nrow = nrow(annot), ncol = length(features))
    colnames(annot_uniprotkb) <- features
    
    lapply(seq(1, nrow(annot)), function(i){ #for each event
      
      eventGr <- GRangesList()
      for (exon in names(grl)){
        eventGr[[paste0(exon)]] <- grl[[paste0(exon)]][i]
      }
      
      protein_ids <- unique(c(annot[i, idx.cols[1]], annot[i, idx.cols[2]]))
      
      lapply(seq_along(features), function(j){ #for each feature
        
        ovl_gr <- overlappingFeatures(features_Gr[[j]], eventGr)
        #ovl_gr_filt <- ovl_gr[ovl_gr$Uniprot_ID %in% protein_ids, ]
        ovl_gr_filt <- ovl_gr
        
        aux <- unique(as.character(ovl_gr_filt$Name))
        aux <- gsub(" ", "", aux)
        res <- paste0(aux, collapse = ",")
        if (!res == ""){
          annot_uniprotkb[i,j] <<- res  
        }
        
      })
      
    })
    
    # for (i in 1:nrow(annot)) {
    #   
    #   # Genomic ranges of alternative splicing event
    #   eventGr <- GRangesList()
    #   for (exon in names(grl)){
    #     eventGr[[paste0(exon)]] <- grl[[paste0(exon)]][i]
    #   }
    #   
    #   protein_ids <- unique(c(annot[i, idx.cols[1]], annot[i, idx.cols[2]]))
    #   
    #   
    #   for (j in 1:length(features)) {
    #     ovl_gr <- overlappingFeatures(features_Gr[[j]], eventGr)
    #     ovl_gr_filt <- ovl_gr[ovl_gr$Uniprot_ID %in% protein_ids, ] 
    #     
    #     aux <- unique(as.character(ovl_gr_filt$Name))
    #     aux <- gsub(" ", "", aux)
    #     res <- paste0(aux, collapse = ",")
    #     if (!res == ""){
    #       annot_uniprotkb[i,j] <- res  
    #     }
    #     
    #   } #all features
    #   
    #   
    # }#all events
    
    #write annotation to maser object
    events_with_features[[paste0(type,"_","events")]] <- 
      cbind(annot, annot_uniprotkb)
    
  }#all types
  
  return(as(events_with_features, "Maser"))
}


mapENSTtoUniprotKB <- function(enst_ids){
  
  if (enst_ids == ""){
    return(c(""))
  }
  
  aux <- strsplit(enst_ids,",")[[1]]
  tokens <- strsplit(aux, "\\.")
  
  enst_trans <- vapply(seq_along(aux),
                       function(i) tokens[[i]][[1]], character(1))
  
  idx.map <- match(enst_trans, UKB_ENST_map$ENST_ID)
  
  return(as.vector(UKB_ENST_map$UniprotKB_ID[idx.map]))
  
}

mapProteinsToEvents <- function(events){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  
  # Add UniprotKB ID
  events_with_ptn <- events
  
  for (type in as_types){
    
    # Retrieve Events annotation
    annot <- events[[paste0(type,"_","events")]]
    idx.cols <- grep("^txn_", colnames(annot))
    
    if (nrow(annot) == 0){
      next
    }
    
    list_ptn_a <- rep("", nrow(annot))
    list_ptn_b <- rep("", nrow(annot))
    
    for (i in 1:nrow(annot)) {
      
      
      list_ptn_a[i] <- paste(mapENSTtoUniprotKB(annot[i,idx.cols[1]]), 
                             collapse = ",")
      list_ptn_b[i] <- paste(mapENSTtoUniprotKB(annot[i,idx.cols[2]]), 
                             collapse = ",")
      
    } #for all events in annot
    
    annot[["list_ptn_a"]] <- list_ptn_a
    annot[["list_ptn_b"]] <- list_ptn_b
    
    events_with_ptn[[paste0(type,"_","events")]] <- annot
    
    
  } #for each AS type
  
  return(events_with_ptn)
  
}


#' Mapping of splice events to Ensembl transcripts.
#' 
#' @param events a maser object.
#' @param gtf a \code{GRanges} object obtained from an Ensembl or Gencode GTF
#'  file using the hg38 build of the human genome.
#' @param ncores number of cores for multithreading (available only in OSX and Linux 
#' machines). If Windows, \code{ncores} will be set to 1 automatically.
#' @return a maser object with transcript and protein identifiers.
#' @details This function performs mapping of splice events in the maser object
#'  to Ensembl transcripts by overlapping exons involved in the splice event to
#'  the transcript models provided in the GTF. 
#'  
#'  Each type of splice event requires a specific mapping procedure
#'   (described below).
#'  
#'  The mapping will also add Uniprot identifiers when the ENST transcript 
#'  encodes for a protein. 
#'    
#'  Visualization of affected transcripts can be done 
#'  using \code{\link{plotTranscripts}}.
#'  
#'   \describe{
#'     \item{\strong{Exon skipping}}{}
#'     \item{Inclusion transcript(s)}{Transcript(s) overlapping the cassette
#'      exon, as well both flanking exons (i.e upstream and downstream exons).}
#'     \item{Skipping transcript(s)}{Transcript(s) overlapping both flanking 
#'     exons but not the cassettte exon.}
#'   }
#'   
#'   \describe{
#'     \item{\strong{Intron retention}}{}
#'     \item{Retention transcript(s)}{Transcript(s) overlapping exactly the
#'      retained intron.}
#'     \item{Skipping transcript(s)}{Transcript(s) where intron is spliced out
#'      and overlapping both flanking exons.}
#'   }
#'   
#'   \describe{
#'   \item{\strong{Mutually exclusive exons}}{}
#'     \item{Exon1 transcript(s)}{Transcript(s) overlapping the first exon 
#'     and both flanking exons.}
#'     \item{Exon2 transcript(s)}{Transcript(s) overlapping the second exon and
#'      both flanking exons.}
#'   }
#'   
#'   \describe{
#'     \item{\strong{Alternative 3' and 5' splice sites}}{}
#'     \item{Short exon transcript(s)}{Transcript(s) overlapping both short and 
#'                        downstream exons.}
#'     \item{Long exon transcript(s)}{Transcript(s) overlapping both long and
#'      downstream exons.}
#'   }
#'   
#' @examples
#' ## Create the maser object
#' path <- system.file("extdata", file.path("MATS_output"), package = "maser")
#' hypoxia <- maser(path, c("Hypoxia 0h", "Hypoxia 24h"))
#' hypoxia_filt <- filterByCoverage(hypoxia, avg_reads = 5)
#' 
#' ## Ensembl GTF annotation for SRSF6
#' gtf_path <- system.file("extdata", file.path("GTF", 
#'  "Ensembl85_examples.gtf.gz"), package = "maser")
#' ens_gtf <- rtracklayer::import.gff(gtf_path)
#'  
#' ## Retrieve gene specific splice events
#' srsf6_events <- geneEvents(hypoxia_filt, geneS = "SRSF6")
#' 
#' ## Map splicing events to transcripts
#' srsf6_mapped <- mapTranscriptsToEvents(srsf6_events, ens_gtf)
#' head(annotation(srsf6_mapped, "SE"))
#' 
#' @seealso \code{\link{plotTranscripts}}
#' @export
#' @import GenomicRanges
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom parallel mclapply
#' @importFrom dplyr inner_join
#' @import methods
#' 
mapTranscriptsToEvents <- function(events, gtf, ncores = 1){
  
  is_strict = TRUE
  
  if(!is(events, "Maser")){
    stop("Parameter events has to be a maser object.")
  }
  
  if (!class(gtf) == "GRanges"){
    stop(cat("\"gtf\" should be a GRanges object."))
  }
  
  if(.Platform$OS.type == "windows"){
    ncores = 1
  }
  
  #Add chr to seqnames - necessary for Gviz plots and compatible with maser()
  if(any(!grepl("chr", GenomeInfoDb::seqlevels(gtf)))){
    GenomeInfoDb::seqlevels(gtf) <- paste0("chr", 
                                           GenomeInfoDb::seqlevels(gtf)) 
  }
  
  gtf_exons <- gtf[gtf$type=="exon",]
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  events <- as(events, "list")
  
  # Add transcripts to events using mapping functions based on sequence overlap 
  events_with_txn <- events
  
  for (type in as_types){
    
    # Retrieve Events annotation
    annot <- events[[paste0(type,"_","events")]]
    
    if (nrow(annot) == 0){
      next
    }
    
    # Retrieve Events gene ranges
    grl <- events[[paste0(type,"_","gr")]]
    
    list_txn <- mclapply(seq_along(annot$ID), function(i) {
      
      # Genomic ranges of alternative splicing events
      eventGr <- lapply(names(grl), function(exon){
        grl[[exon]][i]
      })
      eventGr <- GRangesList(eventGr)
      names(eventGr) <- names(grl)
      
      if(type == "SE") {
        tx_ids <- mapTranscriptsSEevent(eventGr, gtf_exons, is_strict)
        return(data.frame(ID = annot$ID[i],
                          txn_3exons = paste(tx_ids$txn_3exons, collapse = ","),
                          txn_2exons = paste(tx_ids$txn_2exons, collapse = ","),
                          stringsAsFactors = FALSE)
        )
        
      }
      
      if(type == "MXE") {
        tx_ids <- mapTranscriptsMXEevent(eventGr, gtf_exons, is_strict)
        return(data.frame(ID = annot$ID[i],
                          txn_mxe_exon1 = paste(tx_ids$txn_mxe_exon1, collapse = ","),
                          txn_mxe_exon2 = paste(tx_ids$txn_mxe_exon2, collapse = ","),
                          stringsAsFactors = FALSE)
        )
      }
      
      if(type == "RI") {
        tx_ids <- mapTranscriptsRIevent(eventGr, gtf_exons, is_strict)
        return(data.frame(ID = annot$ID[i],
                          txn_nonRetention = paste(tx_ids$txn_nonRetention, collapse = ","),
                          txn_retention = paste(tx_ids$txn_retention, collapse = ","),
                          stringsAsFactors = FALSE)
        )
      }
      
      if(type == "A5SS") {
        
        #reverse strand becomes A3SS
        if (as.character(strand(eventGr[1])) == "+"){
          tx_ids <<- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
        }else{
          tx_ids <<- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
        }
        
        return(data.frame(ID = annot$ID[i],
                          txn_short = paste(tx_ids$txn_short, collapse = ","),
                          txn_long = paste(tx_ids$txn_long, collapse = ","),
                          stringsAsFactors = FALSE)
        )
        
      }
      
      if(type == "A3SS") {
        
        #reverse strand becomes A5SS
        if (as.character(strand(eventGr[1])) == "+"){
          tx_ids <<- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
        }else{
          tx_ids <<- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
        }
        
        return(data.frame(ID = annot$ID[i],
                          txn_short = paste(tx_ids$txn_short, collapse = ","),
                          txn_long = paste(tx_ids$txn_long, collapse = ","),
                          stringsAsFactors = FALSE)
        )
        
      }
      
    }, mc.cores = ncores) #for all events in annot
    
    df_txn <- do.call(rbind, list_txn)
    
    # for (i in 1:nrow(annot)) {
    #   
    #   # Genomic ranges of alternative splicing events
    #   eventGr <- lapply(names(grl), function(exon){
    #     grl[[exon]][i]
    #   })
    #   eventGr <- GRangesList(eventGr)
    #   names(eventGr) <- names(grl)
    #   
    #   if(type == "SE") {
    #     tx_ids <- mapTranscriptsSEevent(eventGr, gtf_exons, is_strict)
    #     list_txn_a[i] <- paste(tx_ids$txn_3exons, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_2exons, collapse = ",")
    #   }
    #   
    #   if(type == "MXE") {
    #     tx_ids <- mapTranscriptsMXEevent(eventGr, gtf_exons, is_strict)
    #     list_txn_a[i] <- paste(tx_ids$txn_mxe_exon1, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_mxe_exon2, collapse = ",")
    #   }
    #   
    #   if(type == "RI") {
    #     tx_ids <- mapTranscriptsRIevent(eventGr, gtf_exons, is_strict)
    #     list_txn_a[i] <- paste(tx_ids$txn_nonRetention, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_retention, collapse = ",")
    #   }
    #   
    #   if(type == "A5SS") {
    #     
    #     #reverse strand becomes A3SS
    #     if (as.character(strand(eventGr[1])) == "+"){
    #       tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
    #     }else{
    #       tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
    #     }
    #     
    #     list_txn_a[i] <- paste(tx_ids$txn_short, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_long, collapse = ",")
    #   }
    #   
    #   if(type == "A3SS") {
    #     
    #     #reverse strand becomes A5SS
    #     if (as.character(strand(eventGr[1])) == "+"){
    #       tx_ids <- mapTranscriptsA3SSevent(eventGr, gtf_exons)  
    #     }else{
    #       tx_ids <- mapTranscriptsA5SSevent(eventGr, gtf_exons)  
    #     }
    # 
    #     list_txn_a[i] <- paste(tx_ids$txn_short, collapse = ",")
    #     list_txn_b[i] <- paste(tx_ids$txn_long, collapse = ",")
    #   }
    #   
    # } #for all events in annot
    
    
    # Update Events annotation
    annot_new <- dplyr::inner_join(annot, df_txn, by = "ID")
    events_with_txn[[paste0(type,"_","events")]] <- annot_new
    
    
  } #for each event type
  
  events_with_ptn <- mapProteinsToEvents(events_with_txn)
  
  return(as(events_with_ptn, "Maser"))
  
}

#' @import GenomicRanges
mapTranscriptsSEevent <- function(eventGr, gtf_exons, is_strict = TRUE){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_upstream, 
                                        gtf_exons, type = "any")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_target, 
                                        gtf_exons, type = "any")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_downstream, 
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  
  #obtain intron range for inclusion event and skipping event
  intron.skipping <- GenomicRanges::GRanges(
    seqnames = seqnames(eventGr$exon_target),
    ranges = IRanges::IRanges(
      start = end(eventGr$exon_upstream) + 1,  
      end = start(eventGr$exon_downstream) - 1),
    strand = strand(eventGr$exon_target)
  )
  
  intron.inclusion <- GenomicRanges::GRanges(
    seqnames = seqnames(eventGr$exon_target),
    ranges = IRanges::IRanges(
      start = c(end(eventGr$exon_upstream) + 1,
                end(eventGr$exon_target) + 1),  
      end = c(start(eventGr$exon_target) - 1,
              start(eventGr$exon_downstream) -1)
    ),
    strand = strand(eventGr$exon_target)
  )
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron.inclusion <- GenomicRanges::findOverlaps(intron.inclusion, 
                                                      gtf_exons, type = "any")
  mytx.ids.intron.inclusion <- 
    gtf_exons$transcript_id[subjectHits(ovl.intron.inclusion)]
  
  ovl.intron.skipping <- GenomicRanges::findOverlaps(intron.skipping, 
                                                     gtf_exons, type = "any")
  mytx.ids.intron.skipping <- 
    gtf_exons$transcript_id[subjectHits(ovl.intron.skipping)]
  
  #decide wich transcripts to plot in inclusion and skipping tracks
  if (is_strict){
    mytx.ids.3exons <- intersect(mytx.ids.e1, mytx.ids.e3)#has both flank exons
    mytx.ids.3exons <- intersect(mytx.ids.3exons, mytx.ids.e2) #and target exon
    mytx.ids.2exons <- intersect(mytx.ids.e1, mytx.ids.e3)
    
  }else {
    mytx.ids.3exons <- union(mytx.ids.e1, mytx.ids.e3)#has either flaking exons
    mytx.ids.3exons <- intersect(mytx.ids.3exons, mytx.ids.e2) #and target exon  
    mytx.ids.2exons <- union(mytx.ids.e1, mytx.ids.e3)
  }
  
  mytx.ids.3exons <- setdiff(mytx.ids.3exons, mytx.ids.intron.inclusion)
  #mytx.ids.2exons <- setdiff(mytx.ids.2exons, mytx.ids.3exons)
  mytx.ids.2exons <- setdiff(mytx.ids.2exons, mytx.ids.intron.skipping)
  
  tx_ids <- list()
  tx_ids[["txn_3exons"]] <- mytx.ids.3exons
  tx_ids[["txn_2exons"]] <- mytx.ids.2exons
  
  return(tx_ids)  
  
}

#' @import GenomicRanges
mapTranscriptsRIevent <- function(eventGr, gtf_exons, is_strict = TRUE){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_upstream, 
                                        gtf_exons, type = "any")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_ir, 
                                        gtf_exons, type = "equal")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_downstream, 
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  
  #obtain intron range from the retention event
  intron <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_ir),
                                   ranges = IRanges::IRanges(
                                     start = end(eventGr$exon_upstream) + 1,  
                                     end = start(eventGr$exon_downstream) - 1),
                                   strand = strand(eventGr$exon_ir)
  )
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron <- GenomicRanges::findOverlaps(intron, gtf_exons, type = "any")
  mytx.ids.intron <- gtf_exons$transcript_id[subjectHits(ovl.intron)]
  
  #decide wich transcripts to plot in retention and non-retention tracks
  if (is_strict){
    #has both upstream and downstream exons
    tx.ids.nonRetention <- intersect(mytx.ids.e1, mytx.ids.e3) 
    
  }else {
    #has either upstream and downstream exons
    tx.ids.nonRetention <- union(mytx.ids.e1, mytx.ids.e3) 
  }
  
  tx.ids.nonRetention <- setdiff(tx.ids.nonRetention, mytx.ids.intron)
  tx.ids.Retention <- mytx.ids.e2
  
  tx_ids <- list()
  tx_ids[["txn_nonRetention"]] <- tx.ids.nonRetention
  tx_ids[["txn_retention"]] <- tx.ids.Retention
  return(tx_ids)
  
}

#' @import GenomicRanges
mapTranscriptsMXEevent <- function(eventGr, gtf_exons, is_strict = TRUE){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_1, 
                                        gtf_exons, type = "any")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_2, 
                                        gtf_exons, type = "any")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_upstream,
                                        gtf_exons, type = "any")
  ovl.e4 <- GenomicRanges::findOverlaps(eventGr$exon_downstream,
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  mytx.ids.e4 <- gtf_exons$transcript_id[subjectHits(ovl.e4)]
  
  #obtain intron range for inclusion event and skipping event
  intron.mxe.exon1 <- GenomicRanges::GRanges(
    seqnames = seqnames(eventGr$exon_1),
    ranges = IRanges::IRanges(
      start = c(end(eventGr$exon_upstream) + 1,
                end(eventGr$exon_1) + 1),  
      end = c(start(eventGr$exon_1) - 1,
              start(eventGr$exon_downstream) -1)
    ),
    strand = strand(eventGr$exon_1)
  )
  
  intron.mxe.exon2 <- GenomicRanges::GRanges(
    seqnames = seqnames(eventGr$exon_2),
    ranges = IRanges::IRanges(
      start = c(end(eventGr$exon_upstream) + 1,
                end(eventGr$exon_2) + 1),  
      end = c(start(eventGr$exon_2) - 1,
              start(eventGr$exon_downstream) -1)
    ),
    strand = strand(eventGr$exon_2)
  )
  
  #find transcripts with exons overlapping intronic regions
  ovl.mxe.exon1 <- GenomicRanges::findOverlaps(intron.mxe.exon1, 
                                               gtf_exons, type = "any")
  mytx.ids.intron1 <- gtf_exons$transcript_id[subjectHits(ovl.mxe.exon1)]
  
  ovl.mxe.exon2 <- GenomicRanges::findOverlaps(intron.mxe.exon2, 
                                               gtf_exons, type = "any")
  mytx.ids.intron2 <- gtf_exons$transcript_id[subjectHits(ovl.mxe.exon2)]
  
  
  #decide wich transcripts to plot in inclusion and skipping tracks
  if (is_strict){
    #has both flanking exons
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.e3, mytx.ids.e4) 
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.mxe.exon1, mytx.ids.e1) #and exon1
    
    #has both flanking exons
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.e3, mytx.ids.e4) 
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.mxe.exon2, mytx.ids.e2) #and exon2
    
  }else {
    #has either flanking exons
    mytx.ids.mxe.exon1 <- union(mytx.ids.e3, mytx.ids.e4) 
    mytx.ids.mxe.exon1 <- intersect(mytx.ids.mxe.exon1, mytx.ids.e1) #and exon 1
    
    #has both flanking exons
    mytx.ids.mxe.exon2 <- union(mytx.ids.e3, mytx.ids.e4) 
    mytx.ids.mxe.exon2 <- intersect(mytx.ids.mxe.exon2, mytx.ids.e2) #and exon2
  }
  
  #remove transcripts with exons in intronic regions
  mytx.ids.mxe.exon1 <- setdiff(mytx.ids.mxe.exon1, mytx.ids.intron1)
  mytx.ids.mxe.exon2 <- setdiff(mytx.ids.mxe.exon2, mytx.ids.intron2)
  
  tx_ids <- list()
  tx_ids[["txn_mxe_exon1"]] <- mytx.ids.mxe.exon1
  tx_ids[["txn_mxe_exon2"]] <- mytx.ids.mxe.exon2
  return(tx_ids)
  
}

#' @import GenomicRanges
mapTranscriptsA5SSevent <- function(eventGr, gtf_exons){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_short, 
                                        gtf_exons, type = "equal")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_long, 
                                        gtf_exons, type = "equal")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_flanking, 
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  
  #obtain intron range for short event and long event
  intron.short <- GenomicRanges::GRanges(
    seqnames = seqnames(eventGr$exon_short),
    ranges = IRanges::IRanges(
      start = end(eventGr$exon_short) + 1,  
      end = start(eventGr$exon_flanking) - 1),
    strand = strand(eventGr$exon_short)
  )
  
  intron.long <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_long),
                                        ranges = IRanges::IRanges(
                                          start = end(eventGr$exon_long) + 1,  
                                          end = start(eventGr$exon_flanking) - 1),
                                        strand = strand(eventGr$exon_long)
  )
  
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron.short <- GenomicRanges::findOverlaps(intron.short, 
                                                  gtf_exons, type = "any")
  mytx.ids.intron.short <- 
    gtf_exons$transcript_id[subjectHits(ovl.intron.short)]
  
  ovl.intron.long <- GenomicRanges::findOverlaps(intron.long, 
                                                 gtf_exons, type = "any")
  mytx.ids.intron.long <- gtf_exons$transcript_id[subjectHits(ovl.intron.long)]
  
  #decide wich transcripts to plot in short and long tracks
  mytx.ids.short <- intersect(mytx.ids.e1, mytx.ids.e3)
  mytx.ids.long <- intersect(mytx.ids.e2, mytx.ids.e3)
  
  #remove transcripts with exons overlapping intronic regions
  mytx.ids.short <- setdiff(mytx.ids.short, mytx.ids.intron.short)
  mytx.ids.long <- setdiff(mytx.ids.long, mytx.ids.intron.long)
  
  tx_ids <- list()
  tx_ids[["txn_short"]] <- mytx.ids.short
  tx_ids[["txn_long"]] <- mytx.ids.long
  return(tx_ids)
  
}

#' @import GenomicRanges
mapTranscriptsA3SSevent <- function(eventGr, gtf_exons){
  
  # Transcripts overlapping with splicing event 
  ovl.e1 <- GenomicRanges::findOverlaps(eventGr$exon_short, 
                                        gtf_exons, type = "equal")
  ovl.e2 <- GenomicRanges::findOverlaps(eventGr$exon_long, 
                                        gtf_exons, type = "equal")
  ovl.e3 <- GenomicRanges::findOverlaps(eventGr$exon_flanking, 
                                        gtf_exons, type = "any")
  
  mytx.ids.e1 <- gtf_exons$transcript_id[subjectHits(ovl.e1)]
  mytx.ids.e2 <- gtf_exons$transcript_id[subjectHits(ovl.e2)]
  mytx.ids.e3 <- gtf_exons$transcript_id[subjectHits(ovl.e3)]
  
  #obtain intron range for short event and long event
  intron.short <- GenomicRanges::GRanges(
    seqnames = seqnames(eventGr$exon_short),
    ranges = IRanges::IRanges(
      start = end(eventGr$exon_flanking) + 1,  
      end = start(eventGr$exon_short) - 1),
    strand = strand(eventGr$exon_short)
  )
  
  intron.long <- GenomicRanges::GRanges(seqnames = seqnames(eventGr$exon_long),
                                        ranges = IRanges::IRanges(
                                          start = end(eventGr$exon_flanking) + 1,  
                                          end = start(eventGr$exon_long) - 1),
                                        strand = strand(eventGr$exon_long)
  )  
  
  #find transcripts with exons overlapping intronic regions
  ovl.intron.short <- GenomicRanges::findOverlaps(intron.short, 
                                                  gtf_exons, type = "any")
  mytx.ids.intron.short <- 
    gtf_exons$transcript_id[subjectHits(ovl.intron.short)]
  
  ovl.intron.long <- GenomicRanges::findOverlaps(intron.long, 
                                                 gtf_exons, type = "any")
  mytx.ids.intron.long <- gtf_exons$transcript_id[subjectHits(ovl.intron.long)]
  
  #decide wich transcripts to plot in short and long tracks
  mytx.ids.short <- intersect(mytx.ids.e1, mytx.ids.e3)
  mytx.ids.long <- intersect(mytx.ids.e2, mytx.ids.e3)
  
  #remove transcripts with exons overlapping intronic regions
  mytx.ids.short <- setdiff(mytx.ids.short, mytx.ids.intron.short)
  mytx.ids.long <- setdiff(mytx.ids.long, mytx.ids.intron.long)
  
  tx_ids <- list()
  tx_ids[["txn_short"]] <- mytx.ids.short
  tx_ids[["txn_long"]] <- mytx.ids.long
  return(tx_ids)
}
