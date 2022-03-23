ValidateVariants = {
    doc title: "GATK ValidateVariants",
        desc:  "This tool validates the adherence of a file to VCF format including information contained within the fields REF, CHR_COUNTS, IDS, ALLELES.",
        constraints: "",
        bpipe_version: "tested with 0.9.9.8.slurm",
        author: "Frank Rühle"

    output.dir = ValidateVariants_vars.outdir

    def ValidateVariants_vars_FLAGS = 
            (ValidateVariants_vars.bwa_ref        ? " -R "      + ValidateVariants_vars.bwa_ref        : "" ) +
            (ValidateVariants_vars.known_variants ? " --dbsnp " + ValidateVariants_vars.known_variants : "" ) +
            (ValidateVariants_vars.extra          ? " "         + ValidateVariants_vars.extra          : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix.prefix).getName())

    transform (".vcf.gz") to (".report") {

        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            gatk --java-options "${ValidateVariants_vars.java_flags}" ValidateVariants $ValidateVariants_vars_FLAGS -V $input > $output
        ""","ValidateVariants"
    }
    forward input
}

