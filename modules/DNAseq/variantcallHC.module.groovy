// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load  PIPELINE_ROOT + "/modules/DNAseq/variantcallHC.vars.groovy"

VariantCallHC = {
    doc title: "GATK Variant Calling HC",
        desc:  "Call variants in BAM files using GATK HaplotypeCaller",
        constraints: "Requires BWA ( paramteter -M ) produced BAM file, with correct chromosome order and ReadGroup attached.",
        bpipe_version: "tested with bpipe 0.9.9.3.slurm",
        author: "Oliver Drechsel"

    output.dir = RESULTS
    def GATK_FLAGS = ""

    // check if a region limit was provided
    if (GATK_CALL_REGION!=null && GATK_CALL_REGION.length()>0) {
        GATK_FLAGS = " -L " + GATK_CALL_REGION
    } else {
        GATK_FLAGS = ""
    }

    def TOOL_ENV = prepare_tool_env("java", tools["java"]["version"], tools["java"]["runenv"]) + " && " +
                   prepare_tool_env("gatk", tools["gatk"]["version"], tools["gatk"]["runenv"])
    def PREAMBLE = get_preamble("VariantCallHC")

    transform (".duprm.realigned.recalibrated.bam") to (".HC.vcf.gz") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            java -Djava.io.tmpdir=\${TMP} -jar \${gatk} -T HaplotypeCaller -nct $GATK_THREADS -R $GATK_BWA_REF --dbsnp $GATK_KNOWN_VARIANTS -I $input -o $output $GATK_FLAGS;
        ""","VariantCallHC"
    }
}

