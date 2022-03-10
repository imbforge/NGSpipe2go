CollectVariantCallingMetrics = {
    doc title: "GATK CollectVariantCallingMetrics",
        desc:  "Collects per-sample and aggregate (spanning all samples) metrics from the provided VCF file.",
        constraints: "",
        bpipe_version: "tested with 0.9.9.8.slurm",
        author: "Frank Rühle"

    output.dir = CollectVariantCallingMetrics_vars.outdir

    def SAMPLENAME = new File(input)
    def SAMPLENAME_BASE = SAMPLENAME.getName()
    def SAMPLENAME_BASE_PRUNED = SAMPLENAME_BASE.replace(".vcf", "").replace(".gz", "") 
    println SAMPLENAME_BASE_PRUNED

    def CollectVariantCallingMetrics_vars_FLAGS = 
            (CollectVariantCallingMetrics_vars.bwa_ref        ? " -R "      + CollectVariantCallingMetrics_vars.bwa_ref        : "" ) +
            (CollectVariantCallingMetrics_vars.known_variants ? " --DBSNP " + CollectVariantCallingMetrics_vars.known_variants : "" ) +
            (CollectVariantCallingMetrics_vars.extra          ? " "         + CollectVariantCallingMetrics_vars.extra          : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix.prefix).getName())

    transform (".vcf.gz") to ("_CollectVariantCallingMetrics.done") {

        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            gatk --java-options "${CollectVariantCallingMetrics_vars.java_flags}" CollectVariantCallingMetrics $CollectVariantCallingMetrics_vars_FLAGS -I $input -O $output.dir/$SAMPLENAME_BASE_PRUNED &&
            touch $output
        ""","CollectVariantCallingMetrics"
    }
    forward input
}

