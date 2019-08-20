// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/NGS/downsamplebam.vars.groovy"

DownsampleBAM = {
    doc title: "DownsampleBAM",
        desc:  "Call samtools tools to downsample a given bam file to roughly a given number of mapped reads",
        constraints: "Samtools tools version >= 1.3",
        bpipe_version: "tested with bpipe 0.9.9.5",
        author: "Nastasja Kreim"

    output.dir = DOWNSAMPLED

    def TOOL_ENV = prepare_tool_env("samtools", tools["samtools"]["version"], tools["samtools"]["runenv"])
    def PREAMBLE = get_preamble("DownsampleBAM")

    transform(".bam") to (".down.bam") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
    
            BASE=\$(basename $input) &&
            samtools view -F 0x04 -bh ${input} -o \${TMP}/\${BASE}_mapped.bam &&
            TOTAL_MAPPED=\$(samtools flagstat \${TMP}/\${BASE}_mapped.bam | grep mapped | head -n 1 | awk '{print \$1 }') &&
            echo mapped_info \$TOTAL_MAPPED &&
            if [[ $DOWNSAMPLE_AMOUNT > \$TOTAL_MAPPED ]]; then
                echo "Downsample amount higher than amount of mapped reads. Keeping all reads!" &&
                cp \${TMP}/\${BASE}_mapped.bam $output;
            else
                PROBABILITY=\$(echo "$DOWNSAMPLE_SEED + $DOWNSAMPLE_AMOUNT/\$TOTAL_MAPPED" | bc -l);
                echo Probability \$PROBABILITY &&
                samtools view -bs \$PROBABILITY -o $output \${TMP}/\${BASE}_mapped.bam;
            fi
        ""","DownsampleBAM"
    }
}

