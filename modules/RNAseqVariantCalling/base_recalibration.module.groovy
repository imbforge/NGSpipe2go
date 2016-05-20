// rule for task GATK BaseRecalibrator from the RNAseq variant calling pipeline
BaseRecalibration = {
   doc title: "GATK BaseRecalibrator",
       desc: "Recalibrate Base Qualities in BAM files, using GATK.",
       constraints: "GATK version >= 3.5",
       author: "Antonio Domingues"

   output.dir = OUTDIR_STAR2ND

   def JAVA_FLAGS = "-Xmx" + RECAL_MAXMEM
   def GATK_FLAGS  = " -R " + GATK_REF +
                     " -nct " + GATK_THREADS

   transform (".bam") to (".recalibration.table", ".recalibrated.bam"){

      exec """
            if [ -n "\$LSB_JOBID" ]; then
               export TMPDIR=/jobdir/\${LSB_JOBID};
            fi &&

            echo 'VERSION INFO'  1>&2 &&
            echo \$(${TOOL_JAVA}/java -jar ${TOOL_GATK}/GenomeAnalysisTK.jar --version) 1>&2 &&
            echo '/VERSION INFO' 1>&2 &&

            ${TOOL_JAVA}/java -jar ${TOOL_GATK}/GenomeAnalysisTK.jar -T BaseRecalibrator -I $input -o $output1 -knownSites ${VCF_REF}  ${GATK_FLAGS} &&

            ${TOOL_JAVA}/java ${JAVA_FLAGS} -jar ${TOOL_GATK}/GenomeAnalysisTK.jar -T PrintReads -I $input -BQSR $output1 -o $output2 ${GATK_FLAGS}

      ""","BaseRecalibration"
   }
}
