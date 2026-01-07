// Job resource limits. Adjust to your needs, though defaults are usually good enough.
// NOTES:
//    * the commands are *sorted*. Please, keep the order!!!
//    * check the documentation here: http://docs.bpipe.org/Language/Config/
//
config {
  executor="slurm"
  short_queue="groups".execute().text =~ /imb-bioinfocf/ ? "bcfshort" : "short"
  long_queue="groups".execute().text =~ /imb-bioinfocf/ ? "bcflong" : "long"
  queue=short_queue   // default queue
  commands {
    AddR {
      walltime="04:00:00"
      procs="1"
      memory="50"
    }
    AddUMIBarcodeToFastq {
      walltime="02:00:00"
      procs="1"
      memory="2"
    }
    assignSouporcellCluster { 
      walltime="03:00:00" 
      procs="1" 
      memory="128"
    }
    Bam2FastQ {
      walltime="1:00:00"
      procs="1"
      memory="2"
    }
    Bam2bw {
      walltime="01:00:00"
      procs="1"
      memory="2"
    }
    Bam2bwStrand {
      walltime="1:00:00"
      procs="1"
      memory="1"
    }
    Bam2bwStrandPE {
      walltime="1:00:00"
      procs="1"
      memory="1"
    }
    bam2bw {
      walltime="01:00:00"
      procs="1"
      memory="2"
    }
    bamCoverage {
      walltime="20:00:00"
      queue=long_queue
      procs="4"
      memory="32"
    }
    BAMindexer {
      walltime="01:00:00"
      procs="1"
      memory="1"
    }
    BamQC {
      walltime="02:00:00"
      procs="1"
      memory="2"
    }
    BaseRecalibration {
      walltime="48:00:00"
      queue=long_queue
      procs="8"
      memory="50"
    }
    excludedRegions_filter {
      walltime="01:00:00"
      procs="1"
      memory="32"
    }
    bowtie1 {
      queue=long_queue
      walltime="24:00:00"
      procs="4"
      memory="24"
    }
    bowtie2 {
      queue=long_queue
      walltime="24:00:00"
      procs="4"
      memory="24"
    }
    bowtie1_sRNA {
      walltime="04:00:00"
      procs="4"
      memory="24"
    }
    BWA_pe {
      queue=long_queue
      walltime="24:00:00"
      procs="16"
      memory="32"
    }
    CatFastQ {
      walltime="1:00:00"
      procs="1"
      memory="2"
    }
    cellranger_count {
      queue=long_queue
      walltime="20:00:00"
      procs="8"
      memory="256"
    }
    cellranger_aggr {
      queue=long_queue
      walltime="20:00:00"
      procs="8"
      memory="256"
    }
    cellrangeratac_count {
      queue=long_queue
      walltime="20:00:00"
      procs="8"
      memory="256"
    }
    cellrangeratac_aggr {
      queue=long_queue
      walltime="20:00:00"
      procs="8"
      memory="256"
    }
    cellrangerarc_count {
      queue=long_queue
      walltime="20:00:00"
      procs="8"
      memory="256"
    }
    cellrangerarc_aggr {
      queue=long_queue
      walltime="20:00:00"
      procs="8"
      memory="256"
    }
    CollectPlots {
      walltime="01:00:00"
      procs="1"
      memory="2"
    }
    CollectVariantCallingMetrics {
      walltime="03:00:00"
      procs="1"
      memory="32"
    }
    CombinedStats {
      walltime="01:00:00"
      procs="1"
      memory="4"
    }
    collectBpipeLogs {
      walltime="00:45:00"
      procs="1"
      memory="64"
    }
    count_breaks {
      walltime="04:00:00"
      procs="4"
      memory="16"
    }
    count_breaks_strandless {
      walltime="01:00:00"
      procs="1"
      memory="1"
    }
    CountNonStrutReads {
      walltime="03:00:00"
      procs="1"
      memory="1"
    }
    CountReadLengths {
      walltime="03:00:00"
      procs="1"
      memory="1"
    }
    CountReads {
      walltime="02:00:00"
      procs="1"
      memory="4"
    }
    CountReadsSummary {
      walltime="04:00:00"
      procs="1"
      memory="50"
    }
    CRmotifCounts {
      walltime="04:00:00"
      procs="1"
      memory="128"
    }
    CTannoSeurat {
      queue=long_queue
      walltime="16:00:00"
      procs="1"
      memory="256"
    }
    CTannoMarker {
      queue=long_queue
      walltime="16:00:00"
      procs="1"
      memory="256"
    }
    Cutadapt {
      walltime="04:00:00"
      procs="1"
      memory="4"
    }
    CutadaptStats {
      walltime="01:00:00"
      procs="1"
      memory="1"
    }
    de_bioc {
      queue=long_queue
      walltime="6:00:00"
      procs="4"
      memory="128"
    }
    DE_DESeq2 {
      walltime="01:00:00"
      procs="1"
      memory="4"
    }
    DE_DESeq2_miRNAmature {
      walltime="01:00:00"
      procs="1"
      memory="4"
    }
    DE_DESeq2_MM {
      walltime="01:00:00"
      procs="1"
      memory="4"
    }
    DE_edgeR {
      walltime="01:00:00"
      procs="1"
      memory="4"
    }
    DedupStats {
      walltime="01:00:00"
      procs="1"
      memory="1"
    }
    demux_gt {
      queue=long_queue
      walltime="12:00:00"
      procs="30"
      memory="64"
    }
    demux_hto { 
      walltime="03:00:00" 
      procs="1" 
      memory="32"
    }
    diffbind2 {
      walltime="01:00:00"
      procs="8"
      memory="64"
    }
    diffbind3 {
      walltime="01:00:00"
      procs="8"
      memory="64"
    }
    diffExprSeurat {
      queue=long_queue
      walltime="12:00:00"
      procs="8"
      memory="128"
    }
    diffPeaks {
      queue=long_queue
      walltime="15:00:00"
      procs="8"
      memory="128"
    }
    dupRadar {
      walltime="05:00:00"
      procs="4"
      memory="10"
    }
    DNAaccess {
      walltime="04:00:00"
      procs="4"
      memory="128"
    }
    extend {
      walltime="01:00:00"
      procs="4"
      memory="8"
    }
    FastQC {
      walltime="02:00:00"
      procs="1"
      memory="2"
    }
    FastQQualityFilter {
      walltime="02:00:00"
      procs="1"
      memory="4"
    }
    FastQQualityFilterStats {
      walltime="01:00:00"
      procs="1"
      memory="1"
    }
    FastqScreen {
      walltime="01:00:00"
      procs="4"
      memory="8"
    }
    FastxTrimmer {
      walltime="02:00:00"
      procs="1"
      memory="4"
    }
    filbowtie2unique {
      walltime="05:00:00"
      procs="8"
      memory="12"
    }
    filtCells_bioc {
      walltime="03:00:00"
      procs="4"
      memory="128"
    }
    Filter2HTSeq {
      walltime="01:00:00"
      procs="1"
      memory="1"
    }
    filter2htseq {
      walltime="00:30:00"
      procs="1"
      memory="1"
    }
    FilterAndMergeSJtab {
      walltime="01:00:00"
      procs="1"
      memory="1"
    }
    FilterChr {
      walltime="02:00:00"
      procs="4"
      memory="1"
    }
    FilterDuplicates {
      walltime="02:00:00"
      procs="1"
      memory="4"
    }
    FilterRNAClasses {
      walltime="1:00:00"
      procs="1"
      memory="1"
    }
    filter_smallRNA_counts {
      walltime="00:30:00"
      procs="1"
      memory="1"
    }
    findmarkers_bioc {
      walltime="4:00:00"
      procs="4"
      memory="128"
    }
    GatherBQSRReports {
      walltime="03:00:00"
      procs="4"
      memory="16"
    }
    geneBodyCov2 {
      walltime="03:00:00"
      procs="4"
      memory="32"
    }
    GenerateStarIndexFromSJ {
      walltime="04:00:00"
      procs="8"
      memory="32"
    }
    GenomicsDBImport {
      walltime="04:00:00"
      procs="8"
      memory="32"
    }
    GenotypeGVCFs {
      walltime="04:00:00"
      procs="8"
      memory="32"
    }
    GO_Enrichment {
      walltime="01:00:00"
      procs="4"
      memory="16"
    }
    GREAT {
      walltime="01:00:00"
      procs="1"
      memory="1"
    }
    grn {
      queue=long_queue
      walltime="24:00:00"
      procs="4"
      memory="254"  
    } 
    hclust_bioc {
      queue=long_queue
      walltime="24:00:00"
      procs="4"
      memory="256"
    }
    HTseqCount {
      walltime="04:00:00"
      procs="1"
      memory="4"
    }
    igraph_bioc {
      queue=long_queue
      walltime="24:00:00"
      procs="4"
      memory="256"
    }
    IndelRealignment {
      walltime="24:00:00"
      queue=long_queue
      procs="8"
      memory="50"
    }
    inferexperiment {
      walltime="02:00:00"
      procs="1"
      memory="4"
    }
    InsertSize {
      walltime="04:00:00"
      procs="1"
      memory="10"
    }
    ipstrength {
      walltime="02:00:00"
      procs="2"
      memory="16"
    }
    kmeans_bioc {
      queue=long_queue
      walltime="24:00:00"
      procs="4"
      memory="256"
    }
    macs2 {
      walltime="04:00:00"
      procs="1"
      memory="8"
    }
    MappingStats {
      walltime="02:00:00"
      procs="1"
      memory="4"
    }
    MappingStatsPlot {
      walltime="02:00:00"
      procs="1"
      memory="4"
    }
    MarkDups {
      walltime="04:00:00"
      procs="1"
      memory="5"
    }
    MarkDups2 {
      walltime="04:00:00"
      procs="1"
      memory="10"
    }
    miRDeep2 {
      walltime="02:00:00"
      procs="2"
      memory="2"
    }
    miRDeep2Mapper {
      walltime="04:00:00"
      procs="8"
      memory="2"
    }
    MirDeep2 {
      walltime="04:00:00"
      procs="2"
      memory="4"
    }
    MirDeep2Mapper {
      walltime="04:00:00"
      procs="4"
      memory="8"
    }
    motifActivity {
      queue=long_queue
      walltime="08:00:00"
      procs="1"
      memory="128"  
    } 
    motifEnrich {
      queue=long_queue
      walltime="08:00:00"
      procs="1"
      memory="128"  
    } 
    motifFootprinting {
      queue=long_queue
      walltime="24:00:00"
      procs="1"
      memory="128"  
    } 
    MULTIQC {
      queue=long_queue
      walltime="12:00:00"
      procs="4"
      memory="128"
    }
    norm_bioc {
      queue=long_queue
      walltime="06:00:00"
      procs="4"
      memory="128"
    }
    NucleotideSignature {
      walltime="1:00:00"
      procs="1"
      memory="25"
    }
    pattern_filtering {
      walltime="04:00:00"
      procs="2"
      memory="1"
    }
    pbc {
      walltime="02:00:00"
      procs="1"
      memory="24"
    }
    peaks2genes {
      queue=long_queue
      walltime="08:00:00"
      procs="1"
      memory="128"  
    } 
    peak_annotation {
      walltime="02:00:00"
      procs="1"
      memory="8"  
    } 
    pear {
      walltime="2:00:00"
      procs="1"
      memory="8"
    }
    phantompeak {
      walltime="05:00:00"
      procs="8"
      memory="16"
    }
    PingPongPro {
      walltime="1:00:00"
      procs="1"
      memory="2"
    }
    PingPongSignal {
      walltime="1:00:00"
      procs="1"
      memory="25"
    }
    PlotReadLengths {
      walltime="1:00:00"
      procs="1"
      memory="1"
    }
    PlotSmallRNAclasses {
      walltime="02:00:00"
      procs="1"
      memory="4"
    }
    qc_bioc {
      walltime="03:00:00"
      procs="4"
      memory="128"
    }
    qualimap {
      walltime="04:00:00"
      procs="1"
      memory="10"
    }
    readAggrData_bioc {
      walltime="04:00:00"
      procs="4"
      memory="128"
    }
    RepEnrich {
      queue=long_queue
      walltime="24:00:00"
      procs="16"
      memory="8"
    }
    RepEnrichPE {
      queue=long_queue
      walltime="24:00:00"
      procs="16"
      memory="8"
    }
    RmDups {
      walltime="04:00:00"
      procs="1"
      memory="8"
    }
    rMATS {
      walltime="03:00:00"
      procs="4"
      memory="8"
    }
    rnaseqc {
      walltime="1:00:00"
      procs="4"
      memory="4"
    }
    rnatypes {
      walltime="0:10:00"
      procs="4"
      memory="4"
    }
    sc_filter {
      walltime="03:00:00"
      procs="4"
      memory="128"
    }
    sc_integrateRNA {
      walltime="04:00:00"
      procs="4"
      memory="128"
    }
    sc_integrateATAC {
      walltime="04:00:00"
      procs="4"
      memory="128"
    }
    sc_qc {
      walltime="03:00:00"
      procs="4"
      memory="128"
    }
    sc_readAggrData {
      walltime="04:00:00"
      procs="4"
      memory="128"
    }
    sc_readIndivSamplesAndMerge{
      walltime="04:00:00"
      procs="4"
      memory="128"
    }
    scType_bioc {
      queue=long_queue
      walltime="06:00:00"
      procs="4"
      memory="128"
    }
    SCTransform {
      walltime="04:00:00"
      procs="4"
      memory="128"
    }
    SelectUnMapped {
      walltime="04:00:00"
      procs="4"
      memory="2"
    }
    SelectUniqMappers {
      walltime="04:00:00"
      procs="4"
      memory="8"
    }
    SequenceBias {
      walltime="0:10:00"
      procs="1"
      memory="4"
    }
    shinyReports {
      walltime="00:05:00"
      procs="1"
      memory="1"
    }
    snpEff {
      walltime="3:00:00"
      procs="1"
      memory="16"
    }
    splitpipe_all {
      queue=long_queue
      walltime="20:00:00"
      procs="8"
      memory="256"
    }
    splitpipe_comb {
      queue=long_queue
      walltime="20:00:00"
      procs="8"
      memory="256"
    }
    SplitNCigarReads {
      queue=long_queue
      procs="8"
      memory="24"
    }
    SplitReadStrands {
      walltime="02:00:00"
      procs="4"
      memory="12"
    }
    STAR {
      walltime="04:00:00"
      procs="4"
      memory="40"
    }
    STAR_pe {
      walltime="04:00:00"
      procs="8"
      memory="6"
    }
    STAR_pe_2nd {
      walltime="04:00:00"
      procs="8"
      memory="32"
    }
    SubreadCount {
      walltime="01:00:00"
      procs="4"
      memory="2"
    }
    subread2rnatypes {
      walltime="02:00:00"
      procs="4"
      memory="10"
    }
    subread_count {
      walltime="00:30:00"
      procs="4"
      memory="4"
    }
    subread_miRNAmature_count {
      walltime="00:30:00"
      procs="4"
      memory="4"
    }
    tpm {
      walltime="00:30:00"
      procs="1"
      memory="4"
    }
    trackhub {
      walltime="00:05:00"
      procs="1"
      memory="1"
    }
    trackhub_config {
      walltime="00:05:00"
      procs="1"
      memory="1"
    }
    TrimUMIs {
      walltime="01:00:00"
      procs="1"
      memory="2"
    }
    upsetPlot {
      queue=long_queue
      walltime="16:00:00"
      procs="1"
      memory="32"
    }
    umi_filtering {
      walltime="01:00:00"
      procs="2"
      memory="4"
    }
    umicount {
      walltime="05:00:00"
      procs="1"
      memory="20"
    }
    umicount_tab {
      walltime="05:00:00"
      procs="1"
      memory="20"
    }
    umidedup {
      walltime="05:00:00"
      procs="1"
      memory="5"
    }
    ValidateVariants {
      walltime="2:00:00"
      procs="8"
      memory="20"
    }
    VariantCallHC {
      walltime="48:00:00"
      queue=long_queue
      procs="2"
      memory="20"
    }
    VariantCallUG {
      walltime="24:00:00"
      queue=long_queue
      procs="8"
      memory="20"
    }
    VariantEval {
      walltime="24:00:00"
      queue=long_queue
      procs="8"
      memory="20"
    }
    VariantFiltration {
      walltime="02:00:00"
      procs="8"
      memory="20"
    }
    VariantFuseHC {
      walltime="24:00:00"
      queue=long_queue
      procs="2"
      memory="20"
    }
    VariantScoreRecalibration {
      walltime="04:00:00"
      procs="1"
      memory="25"
    }
    wnn {
      queue=long_queue
      walltime="08:00:00"
      procs="1"
      memory="128"
    }
  }
}
