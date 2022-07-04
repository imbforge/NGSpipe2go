VariantCallHC = {
    doc title: "GATK Variant Calling HC",
        desc:  "Call variants in BAM files using GATK HaplotypeCaller",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.9.8.slurm",
        author: "Oliver Drechsel, modified by Frank Rühle"

    output.dir = VariantCallHC_vars.outdir

    def HaplotypeCaller_FLAGS =
        (VariantCallHC_vars.erc            ? " -ERC "    + VariantCallHC_vars.erc            : "" ) +
        (VariantCallHC_vars.call_region    ? " -L "      + VariantCallHC_vars.call_region    : "" ) +
        (VariantCallHC_vars.bwa_ref        ? " -R "      + VariantCallHC_vars.bwa_ref        : "" ) +
        (VariantCallHC_vars.known_variants ? " --dbsnp " + VariantCallHC_vars.known_variants : "" ) +
        (VariantCallHC_vars.extra          ? " "         + VariantCallHC_vars.extra          : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform (".bam") to (".g.vcf.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            gatk --java-options "${VariantCallHC_vars.java_flags}" HaplotypeCaller $HaplotypeCaller_FLAGS -I $input -O $output

        ""","VariantCallHC"
    }
}

