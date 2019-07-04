FastxTrimmer = {
	doc title: "Read trimming",
		desc:  "Trims reads up to a given position. Useful when input files with the same read length are needed.",
		constraints: "Only supports compressed FASTQ files",
		author: "Antonio Domingues"

	output.dir = FASTQ_TRIMMER_OUTDIR
   def FASTQ_TRIMMER_FLAGS= " -f "   + FIRST_BASE  +
                            " -l "   + LAST_BASE  +
                            FASTQ_TRIMMER_OTHER

   transform(".fastq.gz") to (".23bp.fastq.gz") {
      exec """
         module load fastx_toolkit/${FASTX_VERSION} &&

         zcat $input | fastx_trimmer $FASTQ_TRIMMER_FLAGS -o $output
      ""","FastxTrimmer"
   }
}
