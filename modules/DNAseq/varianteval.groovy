VariantEval = {
    doc title: "GATK Base Quality Recalibration",
        desc:  "Recalibrate Base Qualities in BAM files, using GATK",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Oliver Drechsel"

    output.dir = VariantEval_vars.outdir

    def VariantEval_FLAGS = 
            (VariantEval_vars.threads        ? " -nt "     + VariantEval_vars.threads        : "" ) +
            (VariantEval_vars.bwa_ref        ? " -R "      + VariantEval_vars.bwa_ref        : "" ) +
            (VariantEval_vars.known_variants ? " --dbsnp " + VariantEval_vars.known_variants : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(module:"VariantEval", branch:branch, branch_outdir:"")

    transform (".vcf.gz") to (".report") {
        // usage parameters https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_varianteval_VariantEval.php
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java ${VariantEval_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T VariantEval $VariantEval_FLAGS --eval $input -o $output
        ""","VariantEval"
    }
    forward input
}

