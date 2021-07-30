miRDeep2Mapper = {
    doc title: "miRDeep2",
        desc:  "Quantification of miRNAs performed in 2 steps: (1) Processes reads and mappping to the reference genome; (2) quantification of miRNA expression.",
        constraints: "Requires mirDeep2.",
        author: "Antonio Domingues, Anke Busch"

    output.dir = miRDeep2Mapper_vars.outdir

    def MIRDEEP2MAPPER_FLAGS=
        (miRDeep2Mapper_vars.genome_ref ? " -p " + miRDeep2Mapper_vars.genome_ref : "") +
        (miRDeep2Mapper_vars.extra      ? " "    + miRDeep2Mapper_vars.extra      : "")

    def TOOL_ENV = prepare_tool_env("mirdeep2", tools["mirdeep2"]["version"], tools["mirdeep2"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".fastq.gz") to (".arf", ".fa") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            x="\${TMP}/\$(basename $input.prefix)" &&
            gzip -cd $input > \$x &&
            mapper.pl \$x $MIRDEEP2MAPPER_FLAGS -s $output2 -t $output1 &> ${output2.prefix}.mapper.log &&
            rm \$x
        ""","miRDeep2Mapper"
    }
}
