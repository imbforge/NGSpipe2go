//rule for task STAR_pe from catalog RNAseq, version 1
//desc: Align single end reads
STAR_pe = {
   doc title: "STAR SE/PE alignment",
      desc:  "Align single end and pair end reads",
      constraints: "Only works with compressed input. Set all global vars.",
      author: "Sergi Sayols modified by N.Kreim and Antonio Domingues"

   output.dir = OUTDIR

   // create the TMP folder if it doesn't exists
   F_TMP = new File(TMP)
   if(! F_TMP.exists()) {
      F_TMP.mkdirs()
   }
   // create the LOGS/STAR folder if it doesn't exists
   F_LOG = new File(LOGS + "/STAR_1stPass")
   if(! F_LOG.exists()) {
      F_LOG.mkdirs()
   }

   // code chunk
   //before we start we have to define the output filename
   // def OUTPUTFILE = input1
   // OUTPUTFILE = (OUTPUTFILE =~ /_R1.fastq.gz/).replaceFirst("")
   def EXP = input1.split("/")[-1].replaceAll("_R1.fastq.gz", "")

   produce(EXP + ".SJ.out.tab") {
      // flags
      // def SAMPLE = new File(OUTPUTFILE)
      def STAR_FLAGS = " --runMode alignReads" +
                " --limitGenomeGenerateRAM " + STAR_MAXRAM +
                " --limitIObufferSize " + STAR_BUFSIZE +
                " --genomeDir " + ESSENTIAL_STAR_REF +
                " --runThreadN " + STAR_THREADS +
                " --outFilterMismatchNmax " + STAR_MM +
                " --outFilterMultimapNmax " + STAR_MULTIMAP +
                " --genomeLoad NoSharedMemory" +
                " --alignIntronMin " + STAR_MININTRO +
                " --outStd SAM" +
                " --outSAMattributes Standard" +
                " --outSJfilterReads Unique" +
                " --outFileNamePrefix " + output.dir + "/" + EXP + "." +
                " --outTmpDir " + TMP + "/" + EXP +
                // " --sjdbGTFfile " + ESSENTIAL_GENESGTF +
                " --readFilesCommand zcat"

      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi                                          &&

         STAR $STAR_FLAGS --readFilesIn $inputs | ${TOOL_SAMTOOLS}/samtools view -bhSu -F 256
- | ${TOOL_SAMTOOLS}/samtools sort -@ $STAR_THREADS - $EXP &&
         rm ${EXP}.bam &&

         rm -rf ${TMP}/${EXP}

      ""","STAR_pe"
   }
}
