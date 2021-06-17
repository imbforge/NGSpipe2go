VariantCallUG = {
    doc title: "GATK Variant Calling UG",
        desc:  "Call variants in BAM files using GATK UnifiedGenotyper",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "Oliver Drechsel"

    output.dir = VariantCallUG_vars.outdir

    def UnifiedGenotyper_FLAGS =
        " -glm BOTH " +
        (VariantCallUG_vars.call_region    ? " -L "      + VariantCallUG_vars.call_region    : "" ) +
        (VariantCallUG_vars.threads        ? " -nt "     + VariantCallUG_vars.threads        : "" ) + // this is not a typo!
        (VariantCallUG_vars.threads        ? " -nct "    + VariantCallUG_vars.threads        : "" ) + // check tool multithread options
        (VariantCallUG_vars.bwa_ref        ? " -R "      + VariantCallUG_vars.bwa_ref        : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(module:"VariantCallUG", branch:branch, branch_outdir:"")

    transform (".dupmarked.realigned.recalibrated.bam") to (".UG.vcf.gz") {
    // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_genotyper_UnifiedGenotyper.php
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${VariantCallUG_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T UnifiedGenotyper $UnifiedGenotyper_FLAGS -I $input -o $output
        ""","VariantCallUG"
    }
}

