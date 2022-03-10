STARRPeaker_makeBin = {
    doc title: "STARRPeaker genomic bins file creation",
        desc:  "STARRPeaker 1_makeBin.py wrapper, for making genomic bins for STARR-seq peak calling. Can be used for both normal STARR-seq and CapSTARR-seq.",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Martin Oti"

    output.dir = STARRPeaker_makeBin_vars.outdir

    def STARRPEAKER_MAKEBIN_FLAGS =
        (STARRPeaker_makeBin_vars.chromsize  ? " --chromsize " + STARRPeaker_makeBin_vars.chromsize  : "") +
        (STARRPeaker_makeBin_vars.blacklist  ? " --blacklist " + STARRPeaker_makeBin_vars.blacklist  : "") +
        (STARRPeaker_makeBin_vars.length  ? " --length " + STARRPeaker_makeBin_vars.length           : "") +
        (STARRPeaker_makeBin_vars.step  ? " --step " + STARRPeaker_makeBin_vars.step                 : "")

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("starrpeaker", tools["starrpeaker"]["version"], tools["starrpeaker"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir)

    produce(STARRPeaker_makeBin_vars.bed) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            
            BASENAME=\$(basename -s .bin.bed $STARRPeaker_makeBin_vars.bed)  &&
            
            1_makeBin.py --prefix \$BASENAME $STARRPEAKER_MAKEBIN_FLAGS      &&
            
            if [[ ${STARRPeaker_makeBin_vars.genesbed} != ""  && -f ${STARRPeaker_makeBin_vars.genesbed} ]];
            then
                bedtools intersect -wa -a <(bedtools slop -b ${STARRPeaker_makeBin_vars.slop} ${STARRPeaker_makeBin_vars.bed}) -b ${STARRPeaker_makeBin_vars.genesbed} | sort -k1,1V > slop_${STARRPeaker_makeBin_vars.bed}  &&
                mv slop_${STARRPeaker_makeBin_vars.bed} ${STARRPeaker_makeBin_vars.bed}
            fi
        ""","STARRPeaker_makeBin"
    }
}

