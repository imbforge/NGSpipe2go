// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

AddUMIBarcodeToFastq = {
    doc title: "Extracts UMI and barcode from read sequence and appends them to the read name",
        desc:  "Extracts UMI and barcode from read sequence and appends them to the read name using umi_tools. Several Designs available via ESSENTIAL_EXPDESIGN" +
                "scRNAseq: extracts string pattern of UMI and barcode from R2 and adds them to R1 header. R2 is discarded (designed for MARS-Seq)" +
                "amplicon1: extracts regular expression of UMI and barcode from a merged fastq file and adds them to read names" +
                "amplicon2: as amplicon1 but without UMIs. Note that the regular expression still needs an UMI segment (of length 0) for compatibility with umi_tools" +
                "amplicon3: nested call of umi_tools to extract two barcodes independently (one assembled read) but no UMIs. umi_tools_FLAGS applies to first umi_tools call and umi_tools_FLAGS_2 to 2nd call. BCs are added to read name separated by underscore." +
                "amplicon4: nested call of umi_tools to extract two barcodes independently (one of each read) but no UMIs. umi_tools_FLAGS applies only to first umi_tools call for read1, umi_tools_FLAGS_2 to R2 in 2nd call. BCs are added to read name separated by underscore.",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Nastasja Kreim, Frank Rühle"

    output.dir = AddUMIBarcodeToFastq_vars.outdir

    def File f = new File(input)
    def OUTPUTFILE = (f.getName() =~ /(.R1)*.fastq.gz/).replaceFirst("")
    def DESIGN = ESSENTIAL_EXPDESIGN
    def EXTRACTWHITELIST  = AddUMIBarcodeToFastq_vars.extractWhitelist
    def EXTRACTWHITELIST2 = AddUMIBarcodeToFastq_vars.extractWhitelist2
    def WRITEFILTEREDOUT  = AddUMIBarcodeToFastq_vars.write_filtered_out

    // create the log folder if it doesn't exists (and whitelist log folder if needed)
    def umitools_logdir = new File(AddUMIBarcodeToFastq_vars.logdir)
    if (!umitools_logdir.exists()) {
        umitools_logdir.mkdirs()
    }
    def umitools_logdirWL = new File(AddUMIBarcodeToFastq_vars.logdirWL)

    def TOOL_ENV = prepare_tool_env("umitools", tools["umitools"]["version"], tools["umitools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())


    def pattern1_FLAGS =
        (AddUMIBarcodeToFastq_vars.extractmethod ? " --extract-method=" + AddUMIBarcodeToFastq_vars.extractmethod : "") +
        (AddUMIBarcodeToFastq_vars.bcpattern     ? " --bc-pattern=" + AddUMIBarcodeToFastq_vars.bcpattern   : "") + 
        (AddUMIBarcodeToFastq_vars.write_filtered_out   ? " --filtered-out=" + "\${TMP}/${OUTPUTFILE}.filteredOut.fastq.gz"   : "")  

    def pattern2_FLAGS =
        (AddUMIBarcodeToFastq_vars.extractmethod ? " --extract-method=" + AddUMIBarcodeToFastq_vars.extractmethod : "") +
        // if PEAR, we have a single assembled read and no R2 to extract
        (AddUMIBarcodeToFastq_vars.bcpattern2    ? (RUN_PEAR ? " --bc-pattern=" : " --bc-pattern2=") + AddUMIBarcodeToFastq_vars.bcpattern2   : "") + 
        (AddUMIBarcodeToFastq_vars.write_filtered_out   ? (RUN_PEAR ? " --filtered-out=" + "\${TMP}/${OUTPUTFILE}.nested.filteredOut.fastq.gz" : "")  : "")  

    def whitelist_FLAGS =
        (AddUMIBarcodeToFastq_vars.method   ? " --method="  + AddUMIBarcodeToFastq_vars.method : "") +
        (AddUMIBarcodeToFastq_vars.knee_method   ? " --knee-method="  + AddUMIBarcodeToFastq_vars.knee_method : "") +
        (AddUMIBarcodeToFastq_vars.set_cell_number   ? " --set-cell-number="  + AddUMIBarcodeToFastq_vars.set_cell_number : "") +
        (AddUMIBarcodeToFastq_vars.expect_cells   ? " --expect-cells="  + AddUMIBarcodeToFastq_vars.expect_cells : "") +
        (AddUMIBarcodeToFastq_vars.error_correct_threshold   ? " --error-correct-threshold="  + AddUMIBarcodeToFastq_vars.error_correct_threshold : "") +
        (AddUMIBarcodeToFastq_vars.ed_above_threshold   ? " --ed-above-threshold="  + AddUMIBarcodeToFastq_vars.ed_above_threshold : "") +
        (AddUMIBarcodeToFastq_vars.extra_umitools_whitelist  ? " "    + AddUMIBarcodeToFastq_vars.extra_umitools_whitelist : "")



    def umi_tools_extract_FLAGS = pattern1_FLAGS +
        (AddUMIBarcodeToFastq_vars.barcodelist   ? " --whitelist="  + AddUMIBarcodeToFastq_vars.barcodelist + " --filter-cell-barcode" + (AddUMIBarcodeToFastq_vars.error_correct_cell ? " --error-correct-cell" : "")  : "") +
        (EXTRACTWHITELIST   ? " --whitelist="  + "\${TMP}/${OUTPUTFILE}.extracted_whitelist.tsv" + " --filter-cell-barcode" + (AddUMIBarcodeToFastq_vars.error_correct_cell ? " --error-correct-cell" : "")  : "") +
        (AddUMIBarcodeToFastq_vars.extra_umitools_extract  ? " "    + AddUMIBarcodeToFastq_vars.extra_umitools_extract : "")

    def umi_tools_extract_FLAGS_2 = pattern2_FLAGS +
        (AddUMIBarcodeToFastq_vars.barcodelist2   ? " --whitelist="  + AddUMIBarcodeToFastq_vars.barcodelist2 + " --filter-cell-barcode" + (AddUMIBarcodeToFastq_vars.error_correct_cell2 ? " --error-correct-cell" : "") : "") +
        (EXTRACTWHITELIST2   ? " --whitelist="  + "\${TMP}/${OUTPUTFILE}.extracted_whitelist2.tsv" + " --filter-cell-barcode" + (AddUMIBarcodeToFastq_vars.error_correct_cell2 ? " --error-correct-cell" : "")  : "") +
        (AddUMIBarcodeToFastq_vars.extra_umitools_extract  ? " "    + AddUMIBarcodeToFastq_vars.extra_umitools_extract : "")


    def umi_tools_whitelist_FLAGS   = pattern1_FLAGS + whitelist_FLAGS
    def umi_tools_whitelist_FLAGS_2 = pattern2_FLAGS + whitelist_FLAGS


    if (DESIGN == "amplicon4") { // for 2nd whitelist extraction: if Readpair as input (i.e. no PEAR assembly), the 2nd whitelist is extracted from R2, otherwise also from R1 (with different barcode def)
        def secondInput = input2.optional;
    } else {
        def secondInput = input;
    }   // def secondInput = input2.optional?:input // not working, blank



    produce(OUTPUTFILE + ".umibarcode.fastq.gz"){
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
      
         if [ $EXTRACTWHITELIST = true ]; then
             umi_tools whitelist $umi_tools_whitelist_FLAGS -I $input -S \${TMP}/${OUTPUTFILE}.extracted_whitelist.tsv --plot-prefix=\${TMP}/${OUTPUTFILE}.plot_extracted_whitelist --log=\${TMP}/${OUTPUTFILE}.extracted_whitelist.umibarcode.log ;
         fi &&
 
         if [ $EXTRACTWHITELIST2 = true ]; then
             umi_tools whitelist $umi_tools_whitelist_FLAGS_2 -I $secondInput -S \${TMP}/${OUTPUTFILE}.extracted_whitelist2.tsv --plot-prefix=\${TMP}/${OUTPUTFILE}.plot_extracted_whitelist2 --log=\${TMP}/${OUTPUTFILE}.extracted_whitelist.umibarcode.log ;
         fi &&
      
         case $DESIGN in
             scRNAseq)   
                       umi_tools extract $umi_tools_extract_FLAGS -I ${input2.optional} --stdout=/dev/null --read2-in $input --read2-out=\${TMP}/${OUTPUTFILE}.umibarcode.fastq.gz --log=\${TMP}/${OUTPUTFILE}.umibarcode.log ;;
 
             amplicon1|amplicon2)  
                       umi_tools extract $umi_tools_extract_FLAGS -I $input --stdout=\${TMP}/${OUTPUTFILE}.umibarcode.fastq.gz --log=\${TMP}/${OUTPUTFILE}.umibarcode.log ;;

             amplicon3)
                       umi_tools extract $umi_tools_extract_FLAGS -I $input --log=\${TMP}/${OUTPUTFILE}.umibarcode.log | 
                       umi_tools extract $umi_tools_extract_FLAGS_2 --stdout=\${TMP}/${OUTPUTFILE}.nested.umibarcode.fastq.gz --log=\${TMP}/${OUTPUTFILE}.nested.umibarcode.log ;;

             amplicon4)
                       umi_tools extract $umi_tools_extract_FLAGS -I $input --stdout=\${TMP}/${OUTPUTFILE}.umibarcode.R1.fastq.gz --read2-in ${input2.optional} --read2-out \${TMP}/${OUTPUTFILE}.umibarcode.R2.fastq.gz --log=\${TMP}/${OUTPUTFILE}.umibarcode.log && 
                       umi_tools extract $umi_tools_extract_FLAGS_2 --stdin=\${TMP}/${OUTPUTFILE}.umibarcode.R1.fastq.gz --stdout=\${TMP}/${OUTPUTFILE}.nested.umibarcode.fastq.gz --read2-in=\${TMP}/${OUTPUTFILE}.umibarcode.R2.fastq.gz --read2-out=/dev/null --log=\${TMP}/${OUTPUTFILE}.nested.umibarcode.log &&
                       rm -f \${TMP}/\$(basename ${OUTPUTFILE}).umibarcode.R?.fastq.gz ;;

             *)        echo 'error: parameter ESSENTIAL_EXPDESIGN not set properly' > $umitools_logdir/${OUTPUTFILE}.umibarcode.log;;
         esac &&
              
          if [ $EXTRACTWHITELIST = true ] || [ $EXTRACTWHITELIST2 = true ]; then
              mkdir -p $umitools_logdirWL &&
              mv -t $umitools_logdirWL/ \${TMP}/${OUTPUTFILE}.extracted_whitelist*.tsv &&
              mv -t $umitools_logdirWL/ \${TMP}/${OUTPUTFILE}*extracted_whitelist*.umibarcode.log &&
              mv -t $umitools_logdirWL/ \${TMP}/${OUTPUTFILE}.plot_extracted_whitelist* ;
          fi &&
         
          if [ $WRITEFILTEREDOUT = true ]; then
              mkdir -p ${output.dir}/filteredOut &&
              mv -t ${output.dir}/filteredOut/ \${TMP}/${OUTPUTFILE}*.filteredOut.fastq.gz ; 
          fi &&

          mv \${TMP}/${OUTPUTFILE}*.umibarcode.fastq.gz $output &&
          mv -t $umitools_logdir/ \${TMP}/${OUTPUTFILE}*umibarcode.log
                                   
        ""","AddUMIBarcodeToFastq"
    }
}



