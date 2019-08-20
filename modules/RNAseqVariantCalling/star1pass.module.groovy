// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/star1pass.vars.groovy"

STAR_pe = {
   doc title: "STAR SE/PE alignment",
      desc:  "Align single end and pair end reads",
      constraints: "Only works with compressed input. Set all global vars.",
      author: "Sergi Sayols modified by N.Kreim and Antonio Domingues"

   output.dir = OUTDIR

   // create the LOGS/STAR folder if it doesn't exists
   def F_LOG = new File(LOGS + "/STAR_1stPass")
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
                // " --sjdbGTFfile " + ESSENTIAL_GENESGTF +
                " --readFilesCommand zcat"

      def TOOL_ENV = prepare_tool_env("star", tools["star"]["version"], tools["star"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
      def PREAMBLE = get_preamble("STAR_pe")

      exec """
         ${TOOL_ENV} &&
         ${PREAMBLE} &&

         STAR $STAR_FLAGS --outTmpDir \${TMP}/$EXP --readFilesIn $inputs | samtools view -bhSu -F 256 - | samtools sort -@ $STAR_THREADS - $EXP &&
         rm ${EXP}.bam &&
         rm -rf \${TMP}/${EXP}
      ""","STAR_pe"
   }
}
