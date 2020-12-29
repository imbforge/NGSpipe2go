FiltermiRNA = {
   doc title: "Filter miRNAs",
         desc:  "Filter alignments to select miRNAs only",
         constraints: "",
         author: "Antonio Domingues"

   output.dir = FILTER_MIRNA_OUTDIR

   transform(".bam") to (".miRNA.bam"){

      exec """
         module load bedtools/${BEDTOOLS_VERSION} &&
         module load htseq/${HTSEQ_VERSION} &&

         python $FILTER_CLASSES_TOOL_PATH -i $input -m 20 -M 24 -o stdout | intersectBed -a stdin -b $FILTER_CLASSES_MIRNA_REF -s -f 1.0 > $output1 &&
         samtools index $output1

      ""","FiltermiRNA"
   }
   forward input
}
