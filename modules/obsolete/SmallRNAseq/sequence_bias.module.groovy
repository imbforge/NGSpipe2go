load MODULE_FOLDER + "SmallRNAseq/sequence_bias.vars.groovy"

SequenceBias = {
   doc title: "Length and first nucleotide frequency ",
         desc:  "Summarizes the length and first nucleotide frequence for mapped reads (small RNAs).",
         constraints: "",
         author: "Antonio Domingues"

   output.dir = BIAS_OUTDIR

   transform(".bam") to (".nuc_bias.txt"){

      exec """
         module load bedtools/${BEDTOOLS_VERSION} &&
         module load pysam/${PYSAM_VERSION} &&

         python ${BIAS_TOOL_PATH} -i $input -o $output

      ""","SequenceBias"
   }
}
