VariantCallHC = {
   doc title: "GATK HaplotypeCaller",
       desc: "Call variants, using GATK HaplotypeCaller.",
       constraints: "GATK version >= 3.5",
       author: "Antonio Domingues"

    output.dir = VariantCallHC_vars.outdir

    def HaplotypeCaller_FLAGS =
        " -dontUseSoftClippedBases" +
        (VariantCallHC_vars.threads        ? " -nct "    + VariantCallHC_vars.threads  : "" ) +
        (VariantCallHC_vars.gatk_ref       ? " -R "      + VariantCallHC_vars.gatk_ref : "" ) +
        (VariantCallHC_vars.vcf_ref        ? " --dbsnp " + VariantCallHC_vars.vcf_ref  : "" ) +
        (VariantCallHC_vars.min_score_call ? " -stand_call_conf " + VariantCallHC_vars.min_score_call : "") +
        (VariantCallHC_vars.min_score_emit ? " -stand_emit_conf " + VariantCallHC_vars.min_score_emit : "")

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble("VariantCallHC")

    transform (".rg.duprm.split.recalibrated.bam") to (".HC.vcf.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${VariantCallHC_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T HaplotypeCaller $HaplotypeCaller_FLAGS -I $input -o $output
        ""","VariantCallHC"
   }
}
