NucleotideSignature = {
   doc title: "Ping-Pong Nucleotide Signature",
         desc:  "Determines the nucleotide bias of small RNA reads. It is expected that libraries enriched in piRNAs will have a typical 5U (Ziwi bound) and/or 10A bias (Zili bound). Generates tables of nucleotide frequency and plots",
         constraints: "",
         author: "Antonio Domingues"

   output.dir = PIRNA_SIGNATURE_OUTDIR
   def SAMPLE_NAME = input.split("/")[-1].replaceAll(".bam", "")
   def FEATURES_NAME = FEATURES_PATH.split("/")[-1].replaceAll(".bed", "")
   def OUTNAME = SAMPLE_NAME + "." + FEATURES_NAME
   def OUT_FOLDER = output.dir + "/" + OUTNAME

   produce(
      OUT_FOLDER + "/figure/" + OUTNAME + ".NucleotideDistributionOnPiRNA.pdf"
      ) {

      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi &&
         source ${TOOL_PYTHONENV} &&
         cd $PIRNA_SIGNATURE_OUTDIR &&

         python $PIRNA_SIGNATURE_TOOL_PATH -f $FASTA_PATH -b $input -u 20 -d 20 -g danRer7 -i $FEATURES_PATH -o $OUT_FOLDER

      ""","NucleotideSignature"
   }
}
