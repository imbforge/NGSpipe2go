VariantFuseHC = {
    doc title: "GATK Fuse variants spread over multiple bam files",
        desc:  "Create single-sample gVCFs from multiple bam files (see: https://gatkforums.broadinstitute.org/gatk/discussion/3893/calling-variants-on-cohorts-of-samples-using-the-haplotypecaller-in-gvcf-mode)",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel"

    output.dir = VariantFuseHC_vars.outdir
    def VariantFuseHC_FLAGS =
        (VariantFuseHC_vars.threads        ? " -nct "    + VariantFuseHC_vars.threads : "" ) +
        (VariantFuseHC_vars.bwa_ref        ? " -R "      + VariantFuseHC_vars.bwa_ref : "" ) +
        (VariantFuseHC_vars.known_variants ? " --dbsnp " + VariantFuseHC_vars.known_variants : "" ) +
        (VariantFuseHC_vars.refconf        ? " --emitRefConfidence "       + VariantFuseHC_vars.refconf   : "" ) +
        (VariantFuseHC_vars.indextype      ? " --variant_index_type "      + VariantFuseHC_vars.indextype : "" ) +
        (VariantFuseHC_vars.indexparm      ? " --variant_index_parameter " + VariantFuseHC_vars.indexparm : "" ) +
        (VariantFuseHC_vars.extra          ? " "         + VariantFuseHC_vars.extra : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, subdir:"", input:new File(input1.prefix).getName())

    transform (".dupmarked.realigned.recalibrated.bam") to (".HC.vcf.gz") {
        // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${VariantFuseHC_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T HaplotypeCaller $VariantFuseHC_FLAGS -I $input -o $output
        ""","VariantFuseHC"
    }
}

