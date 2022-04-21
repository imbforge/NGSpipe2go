// Job resource limits. Adjust to your needs, though defaults are good enough.
// Contents of this file may be overriden through the pipeline's bpipe.config.groovy
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
    BAMindexer { 
      walltime="01:00:00" 
      procs="1" 
      memory="1"
    }
    BWA_pe { 
      queue="bcflong" 
      walltime="24:00:00" 
      procs="16" 
      memory="32"
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
    BamQC { 
      walltime="02:00:00" 
      procs="1" 
      memory="2"
    }
    VariantScoreRecalibration { 
      walltime="04:00:00" 
      procs="1" 
      memory="25"
    }
    BaseRecalibration { 
      walltime="48:00:00" 
      queue="bcflong" 
      procs="8" 
      memory="50"
    }
    Bowtie_se { 
      walltime="02:00:00" 
      procs="8" 
      memory="16"
    }
    CatFastQ { 
      walltime="1:00:00" 
      procs="1" 
      memory="2"
    }
    CollectPlots { 
      walltime="01:00:00" 
      procs="1" 
      memory="2"
    }
    CombinedStats { 
      walltime="01:00:00" 
      procs="1" 
      memory="4"
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
    Cutadapt { 
      walltime="02:00:00" 
      procs="1" 
      memory="4"
    }
    CutadaptStats { 
      walltime="01:00:00" 
      procs="1" 
      memory="1" 
    }
    DE_DESeq2 { 
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
    FastQScreen { 
      walltime="01:00:00" 
      procs="4" 
      memory="8"
    }
    FastxTrimmer { 
      walltime="02:00:00" 
      procs="1" 
      memory="4"
    }
    Filter2HTSeq { 
      walltime="01:00:00" 
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
    GO_Enrichment { 
      walltime="01:00:00" 
      procs="4" 
      memory="8"
    }
    GREAT { 
      walltime="01:00:00" 
      procs="1" 
      memory="1"
    }
    GenerateStarIndexFromSJ { 
      walltime="04:00:00" 
      procs="8" 
      memory="32"
    }
    HTseqCount { 
      walltime="04:00:00" 
      procs="1" 
      memory="4"
    }
    IndelRealignment { 
      walltime="24:00:00" 
      queue="bcflong" 
      procs="8" 
      memory="50"
    }
    InsertSize { 
      walltime="04:00:00" 
      procs="1" 
      memory="10"
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
      memory="50"
    }
    MarkDups2 { 
      walltime="04:00:00" 
      procs="1" 
      memory="10"
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
    MPSprofiling { 
      queue=long_queue 
      walltime="8:00:00" 
      procs="1" 
      memory="128"
    }
    MULTIQC {
      queue="bcflong"
      walltime="12:00:00"
      procs="4"
      memory="50"
    }
    NucleotideSignature { 
      walltime="1:00:00" 
      procs="1" 
      memory="25"
    }
    pear {
      walltime="2:00:00"
      procs="1"
      memory="8"
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
    RepEnrich { 
      queue="long" 
      walltime="24:00:00" 
      procs="16" 
      memory="8"
    }
    RepEnrichPE { 
      queue="long" 
      walltime="24:00:00" 
      procs="16" 
      memory="8"
    }
    RmDups { 
      walltime="02:00:00" 
      procs="1" 
      memory="8"
    }
    STAR { 
      walltime="04:00:00" 
      procs="4" 
      memory="48"
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
    SplitNCigarReads { 
      queue="long" 
      procs="8" 
      memory="24"
    }
    SplitReadStrands { 
      walltime="02:00:00" 
      procs="4" 
      memory="12"
    }
    SubreadCount { 
      walltime="01:00:00" 
      procs="4" 
      memory="2"
    }
    TrimUMIs { 
      walltime="01:00:00" 
      procs="1" 
      memory="2"
    }
    VariantCallHC { 
      walltime="48:00:00" 
      queue="bcflong" 
      procs="2" 
      memory="20"
    }
    VariantCallUG { 
      walltime="24:00:00" 
      queue="bcflong" 
      procs="8" 
      memory="20"
    }
    VariantEval { 
      walltime="24:00:00" 
      queue="bcflong" 
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
      queue="bcflong" 
      procs="2" 
      memory="20"
    }
    bam2bw { 
      walltime="01:00:00" 
      procs="1" 
      memory="2"
    }
    bamCoverage { 
      walltime="03:00:00" 
      procs="8" 
      memory="16"
    }
    blacklist_filter { 
      walltime="01:00:00" 
      procs="1" 
      memory="1"
    }
    bowtie2_pe { 
      queue="long" 
      walltime="24:00:00" 
      procs="4" 
      memory="24"
    }
    bowtie_se { 
      walltime="04:00:00" 
      procs="8" 
      memory="16"
    }
    collectBpipeLogs { 
      walltime="00:05:00" 
      procs="1" 
      memory="1"
    }
    diffbind { 
      walltime="01:00:00" 
      procs="4" 
      memory="8"
    }
    dupRadar { 
      walltime="02:00:00" 
      procs="4" 
      memory="1"
    }
    extend { 
      walltime="01:00:00" 
      procs="4" 
      memory="8"
    }
    filbowtie2unique { 
      walltime="05:00:00" 
      procs="8" 
      memory="12"
    }
    filter2htseq { 
      walltime="00:30:00" 
      procs="1" 
      memory="1"
    }
    geneBodyCov2 { 
      walltime="01:00:00" 
      procs="4" 
      memory="16"
    }
    inferexperiment { 
      walltime="02:00:00" 
      procs="1" 
      memory="4"
    }
    ipstrength { 
      walltime="02:00:00" 
      procs="2" 
      memory="8"
    }
    macs2 { 
      walltime="04:00:00" 
      procs="1" 
      memory="8"
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
    pbc { 
      walltime="02:00:00" 
      procs="1" 
      memory="24"
    }
    peak_annotation { 
      walltime="02:00:00" 
      procs="1" 
      memory="8"   
    }   
    phantompeak { 
      walltime="05:00:00" 
      procs="8" 
      memory="16"
    }
    qualimap { 
      walltime="04:00:00" 
      procs="1" 
      memory="10"
    }
    rnatypes { 
      walltime="0:10:00" 
      procs="4" 
      memory="4"
    }
    shinyReports { 
      walltime="00:05:00" 
      procs="1" 
      memory="1"
    }
    subread2rnatypes { 
      walltime="01:00:00" 
      procs="4" 
      memory="1"
    }
    subread_count { 
      walltime="00:30:00" 
      procs="4" 
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
 
    }
}
