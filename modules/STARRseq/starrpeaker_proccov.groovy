STARRPeaker_procCov = {
    doc title: "STARRPeaker process covariates",
        desc:  "STARRPeaker 2_procCov.py wrapper, for processing sequence covariates of genomic bins for STARR-seq peak calling. Can be used for both normal STARR-seq and CapSTARR-seq.",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Martin Oti"

    output.dir = STARRPeaker_procCov_vars.outdir

    def STARRPEAKER_PROCCOV_FLAGS =
        (STARRPeaker_procCov_vars.chromsize  ? " --chromsize " + STARRPeaker_procCov_vars.chromsize  : "") +
        (STARRPeaker_procCov_vars.blacklist  ? " --blacklist " + STARRPeaker_procCov_vars.blacklist  : "") +
        (STARRPeaker_procCov_vars.length  ? " --length " + STARRPeaker_procCov_vars.length           : "") +
        (STARRPeaker_procCov_vars.step  ? " --step " + STARRPeaker_procCov_vars.step                 : "")

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("starrpeaker", tools["starrpeaker"]["version"], tools["starrpeaker"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir)

    produce(STARRPeaker_procCov_vars.bed) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            
            BASENAME=\$(basename -s .bin.bed $STARRPeaker_procCov_vars.bed)  &&
            
            2_procCov.py --prefix \$BASENAME --cov \$(echo "$STARRPeaker_procCov_vars.covbigwigs" | tr ',' ' ')      &&
            
            if [[ ${STARRPeaker_procCov_vars.linfoldtsv} != ""  && -f ${STARRPeaker_procCov_vars.linfoldtsv} ]];
            then
                paste ${STARRPeaker_procCov_vars.bed} ${STARRPeaker_procCov_vars.linfoldtsv} > lf_${STARRPeaker_procCov_vars.bed}  &&
                mv lf_${STARRPeaker_procCov_vars.bed} ${STARRPeaker_procCov_vars.bed}
            fi
        ""","STARRPeaker_procCov"
    }
}

