//  rule to generate start index using SJ.out.tab, version 1
GenerateStarIndexFromSJ = {
   doc title: "GenerateStarIndexFromSJ",
   desc: "Creates a new STAR genome using SJ identified in a previous mapping step. Part of the 2-step mapping. Based on https://code.google.com/p/rna-star/issues/detail?id=7",
   constraints: "STAR STAR_2.4.2a",
   author: "Antonio Domingues"

   output.dir = OUTDIR_2ND_INDEX

   produce("sjdbInfo.txt"){
      def int OVERHANG
      OVERHANG = ESSENTIAL_READLENGTH.toInteger() - 1

      def STAR_FLAGS = " --runMode genomeGenerate" +
          " --genomeDir " + OUTDIR_2ND_INDEX +
          " --genomeFastaFiles " + ESSENTIAL_GENOME_REF +
          " --sjdbFileChrStartEnd " + SJDBFILE +
          " --sjdbOverhang " + OVERHANG.toString() +
          " --runThreadN " + STAR_THREADS +
          " --limitGenomeGenerateRAM " + STAR_MAXRAM +
          " --limitIObufferSize " + STAR_BUFSIZE +
          // " --outFilterMismatchNmax " + STAR_MM +
          // " --outFilterMultimapNmax " + STAR_MULTIMAP +
          // " --genomeLoad NoSharedMemory" +
          // " --alignIntronMin " + STAR_MININTRO +
          // " --outStd SAM" +
          // " --outSAMattributes Standard" +
          // " --outSJfilterReads Unique" +
          // " --outFileNamePrefix " + output.dir + "/" + EXP + "." +
          // " --outTmpDir " + TMP + "/" + EXP +
          // " --sjdbGTFfile " + ESSENTIAL_GENESGTF +
          " --readFilesCommand zcat"

      exec """
            STAR $STAR_FLAGS
      ""","GenerateStarIndexFromSJ"
   }
}
