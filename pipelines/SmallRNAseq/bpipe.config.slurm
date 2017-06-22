executor="slurm"
queue="short"
commands {
   FastQC {
      walltime="02:00:00"
      procs="1"
      memory="2"
      options="--job-name=FastQC"
   }
   Cutadapt {
      walltime="02:00:00"
      procs="1"
      memory="4"
      options="--job-name=Cutadapt"
   }
   FastQQualityFilter {
      walltime="02:00:00"
      procs="1"
      memory="4"
      options="--job-name=FastQQualityFilter"
   }
   FastxTrimmer {
      walltime="02:00:00"
      procs="1"
      memory="4"
      options="--job-name=FastXtrimmer"
   }
   FilterDuplicates {
      walltime="02:00:00"
      procs="1"
      memory="4"
      options="--job-name=FilterDuplicates"
   }
   Bowtie_se {
      walltime="04:00:00"
      procs="8"
      memory="8"
      options="--job-name=bowtie"
   }
   RepEnrich {
      queue="long"
      procs="8"
      memory="8"
      options="--job-name=RepEnrich"
   }
   RepEnrichPE {
      queue="long"
      procs="8"
      memory="8"
      options="--job-name=RepEnrichPE"
   }
   BAMindexer {
      walltime="01:00:00"
      procs="1"
      memory="1"
      options="--job-name=BAMindexer"
   }
   DedupStats {
      walltime="01:00:00"
      procs="1"
      memory="1"
      options="--job-name=DedupStats"
   }
   MappingStatsPlot {
      walltime="02:00:00"
      procs="1"
      memory="4"
      options="--job-name=MappingStatsPlot"
   }
   CountReads {
      walltime="02:00:00"
      procs="1"
      memory="4"
      options="--job-name=CountReads"
   }
   HTseqCount {
      walltime="04:00:00"
      procs="1"
      memory="4"
      options="--job-name=htseq_count"
     }
   CountReadsSummary {
      queue="long"
      procs="1"
      memory="50"
      options="--job-name=CountReadsSummary"
   }
   SplitReadStrands {
      walltime="02:00:00"
      procs="8"
      memory="3"
      options="--job-name=SplitReadStrands"
   }
   SelectUniqMappers {
      walltime="04:00:00"
      procs="8"
      memory="2"
      options="--job-name=SelectUniqMappers"
   }
   SelectUnMapped {
      walltime="04:00:00"
      procs="8"
      memory="2"
      options="--job-name=SelectUnMapped"
   }
   Bam2bw {
      walltime="01:00:00"
      procs="1"
      memory="2"
      options="--job-name=bam2bw"
   }
   CatFastQ {
      queue="long"
      procs="1"
      memory="2"
      options="--job-name=CatFastQ"
   }
   Bam2FastQ {
      queue="long"
      procs="1"
      memory="2"
      options="--job-name=Bam2FastQ"
   }
   NucleotideSignature {
      procs="1"
      memory="25"
      options="--job-name=NucleotideSignature"
   }
   PingPongSignal {
      procs="1"
      memory="25"
      options="--job-name=PingPongSignal"
   }
   PingPongPro {
      procs="1"
      memory="2"
      options="--job-name=PingPongPro"
   }
   TrimUMIs {
      procs="1"
      memory="2"
      options="--job-name=TrimUMIs"
   }
   MirDeep2Mapper {
      walltime="04:00:00"
      procs="8"
      memory="8"
      options="--job-name=MirDeep2Mapper"
   }
   MirDeep2 {
      walltime="04:00:00"
      procs="2"
      memory="4"
      options="--job-name=mirDeep2"
   }
   CollectPlots {
      walltime="01:00:00"
      procs="1"
      memory="2"
      options="--job-name=CollectPlots"
   }
   FilterRNAClasses {
      procs="1"
      memory="1"
      options="--job-name=FilterRNAClasses"
   }
   SubReadCount {
      walltime="0:30"
      procs="4"
      memory="2"
      options="--job-name=SubReadCount"
   }
      CountReadLengths {
      walltime="03:00:00"
      procs="1"
      memory="1"
      options="--job-name=CountReadLengths"
   }
      PlotReadLengths {
      procs="1"
      memory="1"
      options="--job-name=PlotReadLengths"
   }
      CountNonStrutReads {
      walltime="03:00:00"
      procs="1"
      memory="1"
      options="--job-name=CountNonStrutReads"
   }
}
