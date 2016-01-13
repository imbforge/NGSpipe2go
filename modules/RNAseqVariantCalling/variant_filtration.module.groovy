// rule for task GATK BaseRecalibrator from the RNAseq variant calling pipeline
VariantFiltration = {
   doc title: "GATK HaplotypeCaller",
       desc: "Filter variants following bast practices:http://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq. Note that values are hardcoded.",
       constraints: "GATK version >= 3.5",
       author: "Antonio Domingues"

   output.dir = OUTPUT_HC
   def GATK_FLAGS  = " -R " + GATK_REF

   transform (".vcf.gz") to (".filtered.vcf.gz") {

      exec """
            if [ -n "\$LSB_JOBID" ]; then
               export TMPDIR=/jobdir/\${LSB_JOBID};
            fi &&

            echo 'VERSION INFO'  1>&2 &&
            echo \$(${TOOL_JAVA}/java -jar ${TOOL_GATK}/GenomeAnalysisTK.jar --version) 1>&2 &&
            echo '/VERSION INFO' 1>&2 &&

            ${TOOL_JAVA}/java -jar ${TOOL_GATK}/GenomeAnalysisTK.jar -T VariantFiltration -V $input -o $output -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" ${GATK_FLAGS}

      ""","VariantFiltration"
   }
}
