// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/variantCall_HC.vars.groovy"

VariantCallHC = {
   doc title: "GATK HaplotypeCaller",
       desc: "Call variants, using GATK HaplotypeCaller.",
       constraints: "GATK version >= 3.5",
       author: "Antonio Domingues"

   output.dir = OUTPUT_HC

   def JAVA_FLAGS = "-Xmx" + HC_MAXMEM
   def GATK_FLAGS  = " -R " + GATK_REF +
                     " -nct " + GATK_THREADS +
                     " -stand_call_conf " + MIN_SCORE_CALL +
                     " -stand_emit_conf " + MIN_SCORE_EMIT +
                     " -dontUseSoftClippedBases"

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble("VariantCallHC")

    transform (".rg.duprm.split.recalibrated.bam") to (".UG.vcf.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java $JAVA_FLAGS -jar \${gatk} -T HaplotypeCaller -I $input -o $output --dbsnp $VCF_REF  $GATK_FLAGS
        ""","VariantCallHC"
   }
}
