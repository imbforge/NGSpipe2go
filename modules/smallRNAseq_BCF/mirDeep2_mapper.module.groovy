// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/smallRNAseq_BCF/mirDeep2_mapper.vars.groovy"

miRDeep2Mapper = {
    doc title: "miRDeep2",
        desc:  "Quantification of miRNAs performed in 2 steps: (1) Processes reads and mappping to the reference genome; (2) quantification of miRNA expression.",
        constraints: "Requires mirDeep2.",
        author: "Antonio Domingues, Anke Busch"

    output.dir = MIR_MAPPER_OUTDIR

    def TOOL_ENV = prepare_tool_env("mirdeep2", tools["mirdeep2"]["version"], tools["mirdeep2"]["runenv"])

    transform(".fastq.gz") to (".arf", ".fa") {
        def SAMPLENAME = input.prefix
        def OUTPUTLOG_MAIN = output2.prefix

        exec """
            ${TOOL_ENV} &&

            if [ ! -d ${TMP} ]; then
                mkdir -p ${TMP};
            fi &&

            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
            gunzip -c $input > ${TMP}/\${SAMPLENAME_BASE} &&
            cd $output.dir &&
            mapper.pl ${TMP}/\${SAMPLENAME_BASE} -e -p $GENOME_REF -s $output2 -t $output1 -h -m -i -j -o 8 &> ${OUTPUTLOG_MAIN}.mapper.log &&
            rm ${TMP}/\${SAMPLENAME_BASE}
        ""","miRDeep2Mapper"
    }
}
