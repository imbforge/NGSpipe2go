GenomicsDBImport = {
    doc title: "GATK GenomicsDBImport",
        desc:  "Import single-sample GVCFs into GenomicsDB before joint genotyping.",
        constraints: "input GVCFs must possess genotype likelihoods containing the allele produced by HaplotypeCaller with the '-ERC GVCF' or '-ERC BP_RESOLUTION' settings. At least one interval must be provided. Input GVCFs cannot contain multiple entries for a single genomic position. If GenomicsDBImport cannot be applied use CombineGVCFs instead.",
        bpipe_version: "tested with bpipe 0.9.9.8.slurm",
        author: "Frank Rühle"

    output.dir = GenomicsDBImport_vars.outdir
    TARGETS = PROJECT + PIPELINE_ROOT + "/pipelines/DNAseq/targets.txt"

    def GenomicsDBImport_FLAGS =
        (GenomicsDBImport_vars.sample_map     ? " --sample-name-map " + GenomicsDBImport_vars.outdir + "/" + GenomicsDBImport_vars.sample_map : "" ) +
        (GenomicsDBImport_vars.call_region    ? " -L "      + GenomicsDBImport_vars.call_region    : "" ) +
        (GenomicsDBImport_vars.extra          ? " "         + GenomicsDBImport_vars.extra          : "" )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1).getName())

    produce(GenomicsDBImport_vars.workspace_name) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            rm -f ${GenomicsDBImport_vars.outdir}/${GenomicsDBImport_vars.sample_map} &&
            vcfarray=($inputs.vcf.gz) &&
            for file in "\${vcfarray[@]}"; do
              echo -e \$(echo \$(grep \$(echo \$(basename \${file}) | sed 's/\\..*//') $TARGETS | awk '{print \$1}')'\\t'\${file}) >> ${GenomicsDBImport_vars.outdir}/${GenomicsDBImport_vars.sample_map};
            done &&

            gatk --java-options "${GenomicsDBImport_vars.java_flags}" GenomicsDBImport $GenomicsDBImport_FLAGS --tmp-dir \${TMP} --genomicsdb-workspace-path $output

        ""","GenomicsDBImport"
    }

    forward GenomicsDBImport_vars.outdir + "/" + GenomicsDBImport_vars.workspace_name + "/" + "vcfheader.vcf" // must not forward dir but file
}


