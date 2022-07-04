snpEff = {
    doc title: "snpEff annotation",
        desc:  "Annotates and predicts the effects of genetic variants. Variants marked as not passing filter are excluded.",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8.slurm",
        author: "Frank Rühle"

    output.dir = snpEff_vars.outdir

    def snpEff_FLAGS =
        (snpEff_vars.config        ? " -c "      + snpEff_vars.config        : "" ) +
        (snpEff_vars.output_format ? " -o "      + snpEff_vars.output_format : "" ) +
        (snpEff_vars.extra         ? " "         + snpEff_vars.extra         : "" ) + 
        " " + snpEff_vars.genome_version

    def vcftools_FLAGS =
        (snpEff_vars.output_format=="gatk" ? "--get-INFO EFF " : "--get-INFO ANN " )

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"]) + " && " +
                   prepare_tool_env("snpEff", tools["snpEff"]["version"], tools["snpEff"]["runenv"]) + " && " +
                   prepare_tool_env("vcftools", tools["vcftools"]["version"], tools["vcftools"]["runenv"])

    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix.prefix).getName())


    transform (".vcf.gz") to (".annotated.csv", ".annotated.html", ".annotated.vcf", ".annotated.INFO") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            gatk --java-options "${GenomicsDBImport_vars.java_flags}" SelectVariants --exclude-filtered -V $input -O \${TMP}/excludeFiltered\$(basename ${input}) &&
            java "${snpEff_vars.java_flags}" -jar $TOOL_SNPEFF/snpEff.jar $snpEff_FLAGS \${TMP}/excludeFiltered\$(basename ${input}) -csvStats $output1 -s $output2  > $output3 &&
            vcftools --vcf $output3 --stdout $vcftools_FLAGS > $output4

        ""","snpEff"
    }
}


