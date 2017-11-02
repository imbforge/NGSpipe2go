PlotSensorCoverage = {
    doc title: "PlotSensorCoverage",
    desc: "Plots coverage at a small RNA sensor. See Albuquerque et al. for details.",
    constraints: "none",
    author: "António Domingues"

    output.dir = PLOT_SENSOR_COVERAGE_OUTDIR

    if(NORMALIZATION_TO_NONSTRUCT == "yes") {
        FACS = TRACKS + "/normalization_factors.txt"
    } else {
        FACS = ""
    }

    from("*22G.minus.cov") produce(
        "coverage.replicates.pdf",
        "coverage.SD.pdf"
        ){
        
        exec """
            module load R/${R_VERSION} &&
            
            cd $PLOT_SENSOR_COVERAGE_OUTDIR &&
            Rscript $PLOT_SENSOR_COVERAGE_TOOL_PATH $ESSENTIAL_SENSOR_REF $FACS $inputs

      ""","PlotSensorCoverage"
    }
}
