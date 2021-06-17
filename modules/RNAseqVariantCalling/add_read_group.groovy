AddRG = {
    doc title: "AddReadGroup",
        desc: "Adds reads groups to bam as part of the GATK pipeline",
        constraints: "Picard tools version >= 1.141"
        author: "Antonio Domingues"

    output.dir = AddRG_vars.outdir

    File f = new File(input1)
    def EXP = (f.getName() =~ /.bam/).replaceFirst("")

    def TOOL_ENV = prepare_tool_env("picard", tools["picard"]["version"], tools["picard"]["runenv"])
    def PREAMBLE = get_preamble(module:"AddRG", branch:branch, branch_outdir:"")

    transform(".bam") to (".rg.bam"){
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java $JAVA_FLAGS -jar \${picard} AddOrReplaceReadGroups I=$input O=$output SO=coordinate RGID=${EXP} RGLB=${EXP} RGPL=illumina RGPU=genomics RGSM=${EXP}
        ""","AddRG"
    }
}
