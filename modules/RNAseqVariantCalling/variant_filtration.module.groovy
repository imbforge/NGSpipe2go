// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/variant_filtration.vars.groovy"

VariantFiltration = {
   doc title: "GATK HaplotypeCaller",
       desc: "Filter variants following bast practices:http://gatkforums.broadinstitute.org/gatk/discussion/3891/calling-variants-in-rnaseq. Note that values are hardcoded.",
       constraints: "GATK version >= 3.5",
       author: "Antonio Domingues"

   output.dir = OUTPUT_HC

   def JAVA_FLAGS = "-Xmx" + VARFILT_MAXMEM
   def GATK_FLAGS  = " -R " + GATK_REF

   def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                  prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
   def PREAMBLE = get_preamble("VariantFiltration")

   transform (".vcf.gz") to (".filtered.vcf.gz") {
      exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java $JAVA_FLAGS -jar \${gatk} -T VariantFiltration -V $input -o $output -window 35 -cluster 3 -filterName FS -filter "FS > 30.0" -filterName QD -filter "QD < 2.0" $GATK_FLAGS
      ""","VariantFiltration"
   }
}
