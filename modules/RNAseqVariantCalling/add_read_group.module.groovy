// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseqVariantCalling/add_read_group.vars.groovy"

AddRG = {
    doc title: "AddReadGroup",
        desc: "Adds reads groups to bam as part of the GATK pipeline",
        constraints: "Picard tools version >= 1.141"
        author: "Antonio Domingues"

    File f = new File(input1)
    def EXP = (f.getName() =~ /.bam/).replaceFirst("")

    def TOOL_ENV = prepare_tool_env("picard", tools["picard"]["version"], tools["picard"]["runenv"])
    def PREAMBLE = get_preamble("AddRG")

    transform(".bam") to (".rg.bam"){
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            echo 'VERSION INFO'  1>&2 &&
            echo \$(java -jar \${picard} AddOrReplaceReadGroups --version) 1>&2 &&
            echo '/VERSION INFO' 1>&2 &&

            PLATFORM="genomics" &&
            java $JAVA_FLAGS -jar \${picard} AddOrReplaceReadGroups I=$input O=$output SO=coordinate RGID=${EXP} RGLB=${EXP} RGPL=illumina RGPU=${PLATFORM} RGSM=${EXP}
        ""","AddRG"
    }
}
