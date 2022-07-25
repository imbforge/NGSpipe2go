// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

barcode_count = {
    doc title: "Deduplication and counting of cell barcodes",
        desc: "Counting cell barcode occurrences from fastq read names. If UMIs are present they are used for deduplication.",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.9.8",
        author: "Frank Rühle"

    output.dir = barcode_count_vars.outdir
    def DESIGN = ESSENTIAL_EXPDESIGN

    // create the log folder if it doesn't exists
    def barcode_count_LOGDIR = new File(barcode_count_vars.logdir)
    if (!barcode_count_LOGDIR.exists()) {
        barcode_count_LOGDIR.mkdirs()
    }


    def umicount_FLAGS = 
        (barcode_count_vars.verbose ? "--verbose=1 " : "") +
        (barcode_count_vars.method ? "--method=" + barcode_count_vars.method + " " : "") +
        (barcode_count_vars.edit_distance_threshold ? "--edit-distance-threshold=" + barcode_count_vars.edit_distance_threshold + " " : "") +
        (barcode_count_vars.spliced_is_unique ?  "--spliced-is-unique " : "") +
        (barcode_count_vars.barcode_separator ? "--barcode-separator=" + barcode_count_vars.barcode_separator + " " : "") +
        (barcode_count_vars.param   ? barcode_count_vars.param + " " : "") +
        (barcode_count_vars.extra   ? barcode_count_vars.extra : "")

    def TOOL_ENV = prepare_tool_env("umitools", tools["umitools"]["version"], tools["umitools"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    def BCSEP = barcode_count_vars.barcode_separator


    transform(".fastq.gz") to (".barcode_count.tsv") {
        def SAMPLENAME = input.prefix.prefix

        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&

         case $DESIGN in
            amplicon1)
                zcat $input | awk 'NR%4==1 {print substr(\$1,2)}' | sed -E "s/(.*[0-9]$BCSEP)([ACGTN]+)($BCSEP+)([ACGTN]+)/\\1\\4\\3\\2\\t\${SAMPLENAME_BASE}/" | 
                umi_tools count_tab $umicount_FLAGS -L ${barcode_count_LOGDIR}/\${SAMPLENAME_BASE}.barcode_count.log -E ${barcode_count_LOGDIR}/\${SAMPLENAME_BASE}.barcode_count.error |  
                sed -E -e 's/(^b'\\'')([ACGTN]+)('\\'')(.*)/\\2\\4/' -e '1s/cell/sequence/' -e '1s/gene/id/' > $output;;
            
            amplicon2)
                zcat $input | awk 'NR%4==1 {print substr(\$1,2)}' | sed -E "s/(.*[0-9]$BCSEP)([ACGTN]+)($BCSEP+)/\\2\\t\${SAMPLENAME_BASE}/" |
                awk -F "\\t" 'BEGIN{OFS = "\\t"; print "sequence", "id", "count"} {seen[\$1 "\\t" \$2]++} END{for(x in seen) print x,seen[x]}' > $output;;

            amplicon3|amplicon4)
                zcat $input | awk 'NR%4==1 {print substr(\$1,2)}' | sed -E 's/(.*[0-9]$BCSEP)([ACGTN]+)($BCSEP+)([ACGTN]+)($BCSEP)/\\2\\t\\4/' |
                awk -F "\\t" 'BEGIN{OFS = "\\t"; print "sequence", "id", "count"} {seen[\$1 "\\t" \$2]++} END{for(x in seen) print x,seen[x]}' > $output;;

            *)  echo 'error: parameter ESSENTIAL_EXPDESIGN not set properly' > ${barcode_count_LOGDIR}/\${SAMPLENAME_BASE}.barcode_count.log;;
         esac;    
            
        ""","barcode_count"
    }

    // amplicon1:
    // awk: filter fastq files for lines containing read names and skip the first "@" character
    // sed: switch UMI and CB in readnames from read_id[SEP]_CB_UMI to read_id[SEP]_UMI_CB and add the sample id in 2nd column. This format is needed for umi_tools count_tab
    // umi_tools count_tab: count remaining CB occurrences after UMI deduplication
    // sed: modify count_tab cell column format from b'BC' to BC (escape single quote) and rename columns "cell" and "gene" to "sequence" and "id"

    // amplicon2:
    // awk: filter fastq files for lines containing read names and skip the first "@" character
    // sed: extract CB in readnames and add the sample id in 2nd column (for conformity with amplicon1 design)
    // awk: write header line and count BC occurrences

    // amplicon3|amplicon4:
    // awk: filter fastq files for lines containing read names and skip the first "@" character
    // sed: extract both CB from readnames
    // awk: write header line and count occurrences of BC combinations
}
