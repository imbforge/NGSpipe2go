AddUMIBarcodeToFastq = {
    doc title: "Adds UMI and Barcode of to the fastq header",
        desc:  "adds UMI and barcode of the second read in MARS-Seq samples to the fastq header using umitools",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.2",
        author: "Nastasja Kreim, Frank Rühle"

    output.dir = AddUMIBarcodeToFastq_vars.outdir

    def File f = new File(input1)
    def OUTPUTFILE = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")

    def umi_tools_FLAGS =
        (AddUMIBarcodeToFastq_vars.bcpattern   ? " --bc-pattern=" + AddUMIBarcodeToFastq_vars.bcpattern   : "") +
        (AddUMIBarcodeToFastq_vars.barcodelist ? " --whitelist="  + AddUMIBarcodeToFastq_vars.barcodelist + " --filter-cell-barcode" : "") +
        (AddUMIBarcodeToFastq_vars.extra       ? " "              + AddUMIBarcodeToFastq_vars.extra       : "")

    def TOOL_ENV = prepare_tool_env("umitools", tools["umitools"]["version"], tools["umitools"]["runenv"])
    def PREAMBLE = get_preamble(module:"AddUMIBarcodeToFastq", branch:branch, branch_outdir:"")

    produce(OUTPUTFILE + ".umibarcode.fastq.gz"){
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            umi_tools extract $umi_tools_FLAGS -I $input2 --stdout \${TMP}/\$(basename ${input2.prefix}).barcode.fastq.gz --read2-in $input1 --read2-out=\${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.fastq.gz &&
            
            rm \${TMP}/\$(basename ${input2.prefix}).barcode.fastq.gz &&
            mv \${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.fastq.gz $output
        ""","AddUMIBarcodeToFastq"
    }
}

