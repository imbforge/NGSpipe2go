STARRPeaker_procBam = {
    doc title: "STARRPeaker BAM file processing",
        desc:  "STARRPeaker 3_procBam.py wrapper, for calling STARR-seq peaks. Performs RNA versus DNA peak calling for STARR-seq data. Can be used for both normal STARR-seq and CapSTARR-seq. This step processes the RNA & DNA BAM files prior to peak calling.",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Martin Oti"

    output.dir = STARRPeaker_procBam_vars.outdir + "/$subdir" 

    def STARRPEAKER_PROCBAM_FLAGS =
        (STARRPeaker_procBam_vars.chromsize  ? " --chromsize " + STARRPeaker_procBam_vars.chromsize : "") +
        (STARRPeaker_procBam_vars.bed  ? " --bed " + STARRPeaker_procBam_vars.bed                  : "") +
        (STARRPeaker_procBam_vars.cpus  ? " --cpus " + STARRPeaker_procBam_vars.cpus               : "") +
        (STARRPeaker_procBam_vars.strand ? " --strand " + STARRPeaker_procBam_vars.strand           : "") +
        (STARRPeaker_procBam_vars.min ? " --min " + STARRPeaker_procBam_vars.min                    : "") +
        (STARRPeaker_procBam_vars.max ? " --max " + STARRPeaker_procBam_vars.max                    : "") +
        (STARRPeaker_procBam_vars.se  ? " --se "                                                    : "") +
        (STARRPeaker_procBam_vars.readstart  ? " --readstart "                                      : "")

    def TOOL_ENV = prepare_tool_env("bedtools", tools["bedtools"]["version"], tools["bedtools"]["runenv"]) + " && " +
                   prepare_tool_env("starrpeaker", tools["starrpeaker"]["version"], tools["starrpeaker"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

    transform(".bam") to("_procbam.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            touch $output;
            
            if [ ! -e ${STARRPeaker_procBam_vars.targets} ]; then
                echo "Targets file ${STARRPeaker_procBam_vars.targets} doesn't exist" >> $output &&
                exit 0;
            fi;

            BAM=\$(basename $input) &&           
            extension="\${BAM#*.}" &&
            BAM="\${BAM%%.*}" &&
            grep "^\$BAM" ${STARRPeaker_procBam_vars.targets} | while read -r TARGET; do
                IP=\$(       echo $TARGET | tr '\t' ' ' | cut -f1 -d" ").\$extension &&
                IPname=\$(   echo $TARGET | tr '\t' ' ' | cut -f2 -d" ") &&
                INPUT=\$(    echo $TARGET | tr '\t' ' ' | cut -f3 -d" ").\$extension &&
                INPUTname=\$(echo $TARGET | tr '\t' ' ' | cut -f4 -d" ");

                if [ "\$BAM" != "\$INPUT" ]; then
                    
                    starrpeaker 3_procBam.py --prefix \$BAM --input ${STARRPeaker_procBam_vars.mapped}/\$INPUT --output ${STARRPeaker_procBam_vars.mapped}/\$IP  $STARRPEAKER_PROCBAM_FLAGS ;
                    
                    if [ \$? -ne 0 ]; then rm $output; fi ;
                fi;
            done
        ""","STARRPeaker_procBam"
    }
    
    forward input
}

