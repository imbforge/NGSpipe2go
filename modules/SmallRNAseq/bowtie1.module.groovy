//rule for task bowtie_se from catalog ChIPseq, version 1
//desc: Align single end reads
Bowtie_se = {
   doc title: "Bowtie SE alignment",
      desc:  "Align single end reads",
      constraints: "Only works with compressed input. Samtools multithreaded version expected (>=0.1.19).",
      bpipe_version: "tested with bpipe 0.9.8.7",
      author: "Sergi Sayols"

   output.dir = MAPPED

   def BOWTIE_FLAGS = " -q --sam"  +
                       " "   + BOWTIE_QUALS    +
                       " "   + BOWTIE_BEST     +
                       " -p " + Integer.toString(BOWTIE_THREADS) +
                       " -v " + Integer.toString(BOWTIE_MM)       +
                       " -m " + Integer.toString(BOWTIE_MULTIMAP) +
                       " --trim5 " + Integer.toString(BOWTIE_TRIMM5)  +
                       " --trim3 " + Integer.toString(BOWTIE_TRIMM3)

   transform(".fastq.gz") to (".bam") {

      def SAMPLE_NAME = input.prefix.prefix

      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi                                          &&

         echo 'VERSION INFO'  1>&2 &&
         echo \$(bowtie --version) 1>&2 &&
         echo '/VERSION INFO' 1>&2 &&

         echo "Sample is: ${SAMPLE_NAME}" &&

         zcat $input | $BOWTIE_PATH $BOWTIE_FLAGS $BOWTIE_REF - 2> ${SAMPLE_NAME}.bt.log | samtools view -bhSu - | samtools sort -@ $BOWTIE_THREADS - $output.prefix
      ""","Bowtie_se"
   }
}
