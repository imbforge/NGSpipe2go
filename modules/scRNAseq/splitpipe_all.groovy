splitpipe_all = {
    doc title: "split-pipe alignment",
        desc:  "Align single/paired end reads",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    // calculate name of the sample being processed (added paired end support)
    def File f = new File(input)
    // separate removal of .fastq.gz and .R1 since possible file ending is .R1.cutadapt.fastq.gz,
    // but .cutadapt needs to be kept
    def SAMPLENAME_PRUNED = ((f.getName() =~ /.fastq.gz/).replaceFirst("") =~ /\.R1/).replaceFirst("")
    def SPLITPIPE_INPUT = " --fq1 $input1 --fq2 $input2.optional" 

    output.dir = splitpipe_all_vars.outdir + "/" + SAMPLENAME_PRUNED+ "/" 

    // split-pipe flags
    def SPLITPIPE_FLAGS =
        " --output_dir " + output.dir + 
        " --chemistry "  + splitpipe_all_vars.chemistry + 
        " --genome_dir " + splitpipe_all_vars.genome + 
        " --rseed 100"   + 
        (splitpipe_all_vars.expect_cells ? " --cell_est "      + splitpipe_all_vars.expect_cells : "") + 
        (splitpipe_all_vars.threads      ? " --nthreads "      + splitpipe_all_vars.threads      : "") + 
        (splitpipe_all_vars.extra        ? " "                 + splitpipe_all_vars.extra        : "")

    def SAMTOOLS_SORT_FLAGS = " -O bam" +
        (splitpipe_all_vars.threads ? " -@ " + splitpipe_all_vars.threads : "")


    def TOOL_ENV = prepare_tool_env("split_pipe", tools["split_pipe"]["version"], tools["split_pipe"]["runenv"]) + " && " +
                   prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    produce(SAMPLENAME_PRUNED + ".bam", SAMPLENAME_PRUNED + ".bam.bai") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            grep "$SAMPLENAME_PRUNED" ${splitpipe_all_vars.targets} | awk '{print \$1, \$5}' > $output.dir/${SAMPLENAME_PRUNED}_sample_list.txt &&
            
            split-pipe --mode all $SPLITPIPE_FLAGS $SPLITPIPE_INPUT --samp_list $output.dir/${SAMPLENAME_PRUNED}_sample_list.txt &&
            
            samtools sort --write-index -T \${TMP}/${SAMPLENAME_PRUNED}_sort $SAMTOOLS_SORT_FLAGS -o $output1##idx##$output2 ${output.dir}/process/barcode_headAligned_anno.bam

        ""","splitpipe_all"
    }
}


