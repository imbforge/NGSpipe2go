STARRPeaker_callPeak = {
    doc title: "STARRPeaker peak calling",
        desc:  "STARRPeaker 4_callPeak.py wrapper, for calling STARR-seq peaks. Performs RNA versus DNA peak calling for STARR-seq data. Can be used for both normal STARR-seq and CapSTARR-seq. This is the final peak calling step.",
        constraints: "Should be preceded by running 3_procBam.py",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Martin Oti"

    output.dir = STARRPeaker_callPeak_vars.outdir + "/$subdir" 

    def STARRPEAKER_CALLPEAK_FLAGS =
        (STARRPeaker_callPeak_vars.chromsize  ? " --chromsize " + STARRPeaker_callPeak_vars.chromsize  : "") +
        (STARRPeaker_callPeak_vars.bed  ? " --bed " + STARRPeaker_callPeak_vars.bed                    : "") +
        (STARRPeaker_callPeak_vars.cov  ? " --cov " + STARRPeaker_callPeak_vars.cov                    : "") +
        (STARRPeaker_callPeak_vars.threshold  ? " --threshold " + STARRPeaker_callPeak_vars.threshold  : "") +
        (STARRPeaker_callPeak_vars.mincov  ? " --mincov " + STARRPeaker_callPeak_vars.mincov           : "") +
        (STARRPeaker_callPeak_vars.eq  ? " --eq " + STARRPeaker_callPeak_vars.eq                       : "") +
        (STARRPeaker_callPeak_vars.mode  ? " --mode " + STARRPeaker_callPeak_vars.mode                 : "")

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("starrpeaker", tools["starrpeaker"]["version"], tools["starrpeaker"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".bam") to("_callpeak.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            
            touch $output;
            
            BASENAME=\$(basename -s .bam.bct $input) && 
            4_callPeak.py --prefix \$BASENAME --bct $input --bw \${BASENAME}.output.bw $STARRPEAKER_CALLPEAK_FLAGS ;
            
            if [ \$? -ne 0 ]; then rm $output; fi ;
            
        ""","STARRPeaker_callPeak"
    }
    
    forward input
}

