FilterNonStructuralReads = {
   doc title: "Remove non-structural Reads",
         desc:  "Remove non-structural reads. Created for C. elegans projects, but useful for others",
         constraints: "",
         author: "Antonio Domingues"

   output.dir = FILTER_NONSTUC_OUTDIR

   transform(".bam") to (".non_structural.bam"){

      exec """
         module load bedtools/${BEDTOOLS_VERSION} &&

         bedtools intersect -a $input -b $FILTER_STUCTURAL_REF -v -s -f 0.9 > $output &&

         samtools index $output

      ""","FilterNonStructuralReads"
   }
}
