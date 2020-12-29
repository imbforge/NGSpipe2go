//rule for task bowtie_se from catalog ChIPseq, version 1
//desc: Align single end reads
Bowtie2_se = {
   doc title: "Bowtie2 SE alignment",
      desc:  "Align single end reads",
      constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
      bpipe_version: "tested with bpipe 0.9.8.7",
      author: "Antonio Domingues"

   output.dir = MULTIMAP_OUT_DIR

   def BOWTIE2_FLAGS = " --end-to-end " +
                       " --no-1mm-upfront "  +
                       " -a " +
                       " --no-unal " +
                       " -N "   + BOWTIE2_SEED_MM    +
                       " -L "   + BOWTIE2_SEED_LENGTH     +
                       " -D "   + BOWTIE2_SEED_EXT     +
                       " -R "   + BOWTIE2_RESEED     +
                       " --gbar "   + BOWTIE2_GAPS     +
                       " --rdg "   + BOWTIE2_RGAP_PEN     +
                       " --rfg "   + BOWTIE2_FGAP_PEN     +
                       " --mp "   + BOWTIE2_MM_PEN     +
                       " --score-min "   + BOWTIE2_MIN_SCORE     +
                       " -p " + Integer.toString(BOWTIE2_THREADS) +
                       " -v " + Integer.toString(BOWTIE2_MM)       +
                       " -M " + Integer.toString(BOWTIE2_MULTIREPORT) +
                       " --trim5 " + Integer.toString(BOWTIE2_TRIMM5)  +
                       " --trim3 " + Integer.toString(BOWTIE2_TRIMM3)

   transform(".deduped_barcoded.fastq.gz") to (".bam") {

      def SAMPLE_NAME = input.prefix.prefix

      exec """
         module load bowtie2/${BOWTIE2_VERSION} &&
         module load samtools/${SAMTOOLS_VERSION} &&

         zcat $input | bowtie2 $BOWTIE2_FLAGS -x $BOWTIE2_REF -U - 2> ${SAMPLE_NAME}.bt.log | samtools view -bhSu - | samtools sort -@ $BOWTIE2_THREADS - -o $output
      ""","Bowtie2_se"
   }
}
