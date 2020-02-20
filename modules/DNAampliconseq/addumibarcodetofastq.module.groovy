// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/DNAampliconseq/addumibarcodetofastq.vars.groovy"

AddUMIBarcodeToFastq = {
    doc title: "Adds UMI and Barcode of to the fastq header",
        desc:  "extracts UMI and barcode from read sequence and adds it to fastq header using umi_tools. Several Designs available via ESSENTIAL_EXPDESIGN" +
                "scRNAseq: extracts string pattern of UMI and barcode from R2 and adds them to R1 header. R2 is discarded (designed for MARS-Seq)" +
                "amplicon1: extracts regular expression of UMI and barcode from a merged fastq and adds them to header",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.2",
        author: "Nastasja Kreim, Frank Rühle"

    output.dir = AddUMIBarcodeToFastq_vars.outdir

    def File f = new File(input)
    def OUTPUTFILE = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")
    def DESIGN = ESSENTIAL_EXPDESIGN

    // create the log folder if it doesn't exists
    def umitools_logdir = new File(AddUMIBarcodeToFastq_vars.logdir)
    if (!umitools_logdir.exists()) {
        umitools_logdir.mkdirs()
    }

    def umi_tools_FLAGS =
        (AddUMIBarcodeToFastq_vars.bcpattern     ? " --bc-pattern=" + AddUMIBarcodeToFastq_vars.bcpattern   : "") +
        (AddUMIBarcodeToFastq_vars.bcpattern2    ? " --bc-pattern2=" + AddUMIBarcodeToFastq_vars.bcpattern2   : "") +
        (AddUMIBarcodeToFastq_vars.extractmethod ? " --extract-method=" + AddUMIBarcodeToFastq_vars.extractmethod : "") +
        (AddUMIBarcodeToFastq_vars.barcodelist   ? " --whitelist="  + AddUMIBarcodeToFastq_vars.barcodelist + " --filter-cell-barcode" : "") +
        (AddUMIBarcodeToFastq_vars.extra         ? " "              + AddUMIBarcodeToFastq_vars.extra       : "")

    def TOOL_ENV = prepare_tool_env("umitools", tools["umitools"]["version"], tools["umitools"]["runenv"])
    def PREAMBLE = get_preamble("AddUMIBarcodeToFastq")

    produce(OUTPUTFILE + ".umibarcode.fastq.gz"){
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

         case $DESIGN in
             scRNAseq)   
                       umi_tools extract $umi_tools_FLAGS -I ${input2.optional} --stdout=/dev/null --read2-in $input --read2-out=\${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.fastq.gz --log=\${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.log &&
                       mv \${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.fastq.gz $output &&
                       mv \${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.log $umitools_logdir/\$(basename ${OUTPUTFILE}).umibarcode.log;;

             amplicon1|amplicon2)  
                       umi_tools extract $umi_tools_FLAGS -I $input --stdout=\${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.fastq.gz --log=\${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.log &&
                       mv \${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.fastq.gz $output &&
                       mv \${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.log $umitools_logdir/\$(basename ${OUTPUTFILE}).umibarcode.log;;

             *)        echo 'error: parameter ESSENTIAL_EXPDESIGN not set properly' > $umitools_logdir/\$(basename ${OUTPUTFILE}).umibarcode.log;;
         esac;     
                                   
        ""","AddUMIBarcodeToFastq"
    }
}

