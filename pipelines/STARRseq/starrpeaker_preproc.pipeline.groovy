PIPELINE="STARRPeaker_preproc"
PIPELINE_VERSION="1.0"
PIPELINE_ROOT="./NGSpipe2go/"

load PIPELINE_ROOT + "/pipelines/STARRseq/essential.vars.groovy"
load PIPELINE_ROOT + "/pipelines/STARRseq/tools.groovy"
load PIPELINE_ROOT + "/config/preambles.groovy"
load PIPELINE_ROOT + "/config/bpipe.config.groovy"

load PIPELINE_ROOT + "/modules/STARRseq/starrpeaker_makebin.header"
load PIPELINE_ROOT + "/modules/STARRseq/starrpeaker_calcfoldingenergy.header"
load PIPELINE_ROOT + "/modules/STARRseq/starrpeaker_proccov.header"

//MAIN PIPELINE TASK

Bpipe.run {
    STARRPeaker_makeBin + STARRPeaker_calcFoldingEnergy + STARRPeaker_procCov
}


