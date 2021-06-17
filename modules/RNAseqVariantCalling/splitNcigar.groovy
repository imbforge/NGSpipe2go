SplitNCigarReads = {
    doc title: "GATK SplitNCigarReads",
        desc: "Splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions",
        constraints: "GATK version >= 3.5",
        author: "Antonio Domingues"

    output.dir = SplitNCigarReads_vars.outdir

    def SplitNCigarReads_FLAGS =
    (SplitNCigarReads_vars.gatk_ref         ? " -R "    + SplitNCigarReads_vars.gatk_ref         : "") +
    (SplitNCigarReads_vars.read_filter_flag ? " -rf "   + SplitNCigarReads_vars.read_filter_flag : "") +
    (SplitNCigarReads_vars.map_q_from_flag  ? " -RMQF " + SplitNCigarReads_vars.map_q_from_flag  : "") +
    (SplitNCigarReads_vars.map_q_to_flag    ? " -RMQT " + SplitNCigarReads_vars.map_q_to_flag    : "") +
    (SplitNCigarReads_vars.unsafe_flag      ? " -U "    + SplitNCigarReads_vars.unsafe_flag      : "")

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, subdir:"", input:new File(input1.prefix).getName())

    transform (".duprm.bam") to (".duprm.split.bam"){
       exec """
           ${TOOL_ENV} &&
           ${PREAMBLE} &&

           java ${VariantCallHC_vars.java_flags} -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T SplitNCigarReads $SplitNCigarReads_FLAGS -I $input -o $output
       ""","SplitNCigarReads"
    }
}
