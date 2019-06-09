//rule for task bowtie_se from catalog ChIPseq, version 1
//desc: Align single end reads
Bowtie_se = {
   doc title: "Bowtie SE alignment",
      desc:  "Align single end reads",
      constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
      bpipe_version: "tested with bpipe 0.9.8.7",
      author: "Sergi Sayols"

   output.dir = MULTIMAP_OUT_DIR

   def BOWTIE_FLAGS = " -q --sam"  +
                       " "   + BOWTIE_QUALS    +
                       " "   + BOWTIE_BEST     +
                       " -p " + Integer.toString(BOWTIE_THREADS) +
                       " -v " + Integer.toString(BOWTIE_MM)       +
                       " -M " + Integer.toString(BOWTIE_MULTIREPORT) 

   transform(".cutadapt.highQ.deduped.trimmed.fastq.gz") to (".bam") {

      def SAMPLE_NAME = output.prefix

      exec """
            module load bowtie/${BOWTIE_VERSION} &&
            module load samtools/${SAMTOOLS_VERSION} &&
            
      if [ ! -e $TMP ]; then
        mkdir -p $TMP;
      fi &&

      SAMPLENAME_BASE=\$(basename ${SAMPLE_NAME}) &&

      echo 'BOWTIE_FLAGS' $BOWTIE_FLAGS > $output.dir/\${SAMPLENAME_BASE}.bowtie.log &&
      echo 'BOWTIE_REF' $BOWTIE_REF >> $output.dir/\${SAMPLENAME_BASE}.bowtie.log && 

      zcat $input | bowtie $BOWTIE_FLAGS $BOWTIE_REF - 2>> $output.dir/\${SAMPLENAME_BASE}.bowtie.log | awk '{if (\$1~/^@/) print; else {if(\$5 == 255) print \$0"\tNH:i:1"; else print \$0"\tNH:i:2";}}' | samtools view -bhSu - | samtools sort -@ $BOWTIE_THREADS -T $TMP/\$(basename $output.prefix)_bowtie1_sort - -o $output

      ""","Bowtie_se"
   }
}
