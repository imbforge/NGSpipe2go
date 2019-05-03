//rule for task 
//desc: count the reads per gene and split based on cell barcode 
umicount = {
    doc title: "Counting reads per gene",
        desc: "Counting of mapped data and splitting accoring to cellbarcode with umi_tools",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Nastasja Kreim"

    def UMICOUNT_FLAGS =    UMICOUNT_LOG + " " +
                            UMICOUNT_PARAM + " " +
                            UMICOUNT_EXTRA
    output.dir = UMICOUNT_OUTDIR
    if(ESSENTIAL_PAIRED == "yes"){
      UMICOUNT_FLAGS = UMICOUNT_FLAGS + " --paired"
   }
    
    // run the chunk
    transform(".bam\$") to (".umicount.tsv.gz") {
        exec """
            module load umitools/${UMITOOLS_VERSION} &&
            umi_tools count $UMICOUNT_FLAGS -I $input -S $output1   
        ""","umicount"
    }
}
