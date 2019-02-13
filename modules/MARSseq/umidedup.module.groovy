//rule for task geneBodyCov from catalog RNAseq, version 1
//desc: Calculate the RNA-seq coverage over gene body
umidedup = {
    doc title: "deduplication based on UMIs",
        desc: "Deduplication of mapped data using UMIs with umi_tools",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.3",
        author: "Nastasja Kreim"

    output.dir = UMIDEDUP_OUTDIR
    def UMIDEDUP_FLAGS =    UMIDEDUP_LOG + " " +
                            UMIDEDUP_PARAM + " " +
                            UMIDEDUP_EXTRA
    
    if(ESSENTIAL_PAIRED == "yes"){
      UMIDEDUP_FLAGS = UMIDEDUP_FLAGS + " --paired"
   }
            //umi_tools dedup $UMIDEDUP_FLAGS -I $input -S $output1 -E $output2 -L $output3 --output-stats=${output1.prefix}.stats   
    
    // run the chunk
    //transform(".bam") to (".umidedup.bam", ".umidedup.err", ".umidedup.log") {
    transform(".bam") to (".umidedup.bam") {
        exec """
            module load umitools/${UMITOOLS_VERSION} &&
            umi_tools dedup $UMIDEDUP_FLAGS -I $input -S $output1 --output-stats=${output1.prefix}.stats   
        ""","umidedup"
    }
}
