// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/splitNcigar.vars.groovy"

SplitNCigarReads = {
    doc title: "GATK SplitNCigarReads",
        desc: "Splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the intronic regions",
        constraints: "GATK version >= 3.5",
        author: "Antonio Domingues"

    output.dir = OUTDIR_STAR2ND

    def JAVA_FLAGS = "-Xmx" + SPLITNCIGAR_MAXMEM
    def SPLITCIGAR_FLAGS  = " -R " + GATK_REF +
                            " -rf " + READ_FILTER_FLAG +
                            " -RMQF " + MAP_Q_FROM_FLAG +
                            " -RMQT " + MAP_Q_TO_FLAG +
                            " -U " + UNSAFE_FLAG

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble("SplitNCigarReads")

    transform (".duprm.bam") to (".duprm.split.bam"){
       exec """
           ${TOOL_ENV} &&
           ${PREAMBLE} &&

           java $JAVA_FLAGS -jar \${gatk} -T SplitNCigarReads -I $input -o $output $SPLITCIGAR_FLAGS
       ""","SplitNCigarReads"
    }
}
