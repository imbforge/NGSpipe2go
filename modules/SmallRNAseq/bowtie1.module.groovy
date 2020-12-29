//rule for task bowtie_se from catalog ChIPseq, version 1
//desc: Align single end reads
Bowtie_se = {
   doc title: "Bowtie SE alignment",
      desc:  "Align single end reads",
      constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
      bpipe_version: "tested with bpipe 0.9.8.7",
      author: "Sergi Sayols and modified by Antonio Domingues"

   output.dir = MULTIMAP_OUT_DIR

   def BOWTIE_FLAGS = " -q --sam"  +
                       " "   + BOWTIE_QUALS    +
                       " "   + BOWTIE_BEST     +
                       " -p " + Integer.toString(BOWTIE_THREADS) +
                       " -v " + Integer.toString(BOWTIE_MM)       +
                       " -M " + Integer.toString(BOWTIE_MULTIREPORT) +
                       " --trim5 " + Integer.toString(BOWTIE_TRIMM5)  +
                       " --trim3 " + Integer.toString(BOWTIE_TRIMM3)

  // keep only the file name without any extensions.
   def SAMPLE_NAME = input.split("/")[-1].split("\\.")[0]
   
   produce(SAMPLE_NAME + ".bam") {
   branch.totalBams = output

      exec """
         module load bowtie/${BOWTIE_VERSION}     &&
         module load samtools/${SAMTOOLS_VERSION} &&

         zcat $input | bowtie $BOWTIE_FLAGS $BOWTIE_REF - 2> $output.dir/${SAMPLE_NAME}.bt.log | samtools view -bhSu - | samtools sort -@ $BOWTIE_THREADS - -o $output
      ""","Bowtie_se"
   }
}
