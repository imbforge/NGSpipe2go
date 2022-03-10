GenotypeGVCFs = {
    doc title: "GATK GenotypeGVCFs",
        desc:  "Perform joint genotyping on single input file containing one or more samples pre-called with HaplotypeCaller",
        constraints: "input sample file must possess genotype likelihoods produced by HaplotypeCaller with '-ERC GVCF' or '-ERC BP_RESOLUTION'",
        bpipe_version: "tested with bpipe 0.9.9.8.slurm",
        author: "Frank Rühle"

    output.dir = GenotypeGVCFs_vars.outdir

    def GenotypeGVCFs_FLAGS =
        (GenotypeGVCFs_vars.bwa_ref        ? " -R "      + GenotypeGVCFs_vars.bwa_ref        : "" ) +
        (GenotypeGVCFs_vars.extra          ? " "         + GenotypeGVCFs_vars.extra          : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1).getName())

    println "Used workspace from GenomicsDBImport: " + input.dir

    produce(GenotypeGVCFs_vars.vcf_name) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            gatk --java-options "${GenotypeGVCFs_vars.java_flags}" GenotypeGVCFs $GenotypeGVCFs_FLAGS --tmp-dir \${TMP} -V gendb://$input.dir -O $output

        ""","GenotypeGVCFs"
    }
}


