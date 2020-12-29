FilterSensorClasses = {
   doc title: "Filter small RNA classes",
         desc:  "Filter alignments to select C. elegans small RNA classes: 21U, 22G, 26G",
         constraints: "",
         author: "Antonio Domingues"

   output.dir = FILTER_CLASSES_OUTDIR

   transform(".bam") to (".sensor.22G.bam"){

      exec """
         module load bedtools/${BEDTOOLS_VERSION} &&
         module load htseq/${HTSEQ_VERSION} &&

         python $FILTER_CLASSES_TOOL_PATH -i $input -m 22 -M 22 -o stdout | intersectBed -a stdin -b $FILTER_CLASSES_SENSOR_REF -S > $output &&

         samtools index $output1 

      ""","FilterSensorClasses"
   }
   forward input, output
}
