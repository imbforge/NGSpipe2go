STARRPeaker_calcFoldingEnergy = {
    doc title: "STARRPeaker utility script to calculate RNA secondary structure in (genomewide) bins",
        desc:  "STARRPeaker calcFoldingEnergy.py wrapper, for calculating RNA secondary structure in (genomewide) bins, to use as a covariate for STARR-seq peak calling.",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Martin Oti"

    output.dir = STARRPeaker_calcFoldingEnergy_vars.outdir

    def STARRPEAKER_CALCFOLDINGENERGY_FLAGS =
        (STARRPeaker_calcFoldingEnergy_vars.bed ? " --bed " + STARRPeaker_calcFoldingEnergy_vars.bed           : "") +
        (STARRPeaker_calcFoldingEnergy_vars.out ? " --out " + STARRPeaker_calcFoldingEnergy_vars.out           : "") +
        (STARRPeaker_calcFoldingEnergy_vars.genome ? " --genome " + STARRPeaker_calcFoldingEnergy_vars.genome  : "") +
        (STARRPeaker_calcFoldingEnergy_vars.linearfold ? " --linearfold " + STARRPeaker_calcFoldingEnergy_vars.linearfold : "") +
        (STARRPeaker_calcFoldingEnergy_vars.cpus ? " --cpus " + STARRPeaker_calcFoldingEnergy_vars.cpus                           : "")

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("starrpeaker", tools["starrpeaker"]["version"], tools["starrpeaker"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir)

    produce(STARRPeaker_calcFoldingEnergy_vars.out) {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
            
            calcFoldingEnergy.py $STARRPEAKER_CALCFOLDINGENERGY_FLAGS ;
            
        ""","STARRPeaker_calcFoldingEnergy"
    }
}

