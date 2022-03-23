VariantEval = {
    doc title: "GATK VariantEval",
        desc:  "Variant evaluation (% in dbSNP, genotype concordance, Ti/Tv ratios, and a lot more)",
        constraints: "VariantEval is a BETA tool and is not yet ready for use in production",
        bpipe_version: "tested with 0.9.9.8.slurm",
        author: "Oliver Drechsel, modified by Frank Rühle"

    output.dir = VariantEval_vars.outdir

    def VariantEval_FLAGS = 
            (VariantEval_vars.bwa_ref        ? " -R "      + VariantEval_vars.bwa_ref        : "" ) +
            (VariantEval_vars.known_variants ? " --dbsnp " + VariantEval_vars.known_variants : "" ) +
            (VariantEval_vars.extra          ? " "         + VariantEval_vars.extra          : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix.prefix).getName())

    transform (".vcf.gz") to (".report") {

        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            gatk --java-options "${VariantEval_vars.java_flags}" VariantEval $VariantEval_FLAGS --eval $input -O $output
        ""","VariantEval"
    }
    forward input
}

