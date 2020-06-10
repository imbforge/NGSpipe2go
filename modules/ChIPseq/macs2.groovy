macs2 = {
    doc title: "MACS2",
        desc:  "MACS2 wrapper",
        constraints: "Only performs treatment versus control peak calling",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = macs2_vars.outdir

    def MACS2_FLAGS =
        (macs2_vars.gsize  ? " -g "           + macs2_vars.gsize  : "") +
        (macs2_vars.bwidth ? " --bw "         + macs2_vars.bwidth : "") +
        (macs2_vars.minlen ? " --min-length " + macs2_vars.minlen : "") +
        (macs2_vars.paired ? " --format BAMPE"                    : "") +
        (macs2_vars.extra  ? " " + macs2_vars.extra               : "")

    def TOOL_ENV = prepare_tool_env("macs2", tools["macs2"]["version"], tools["macs2"]["runenv"])
    def PREAMBLE = get_preamble("macs2")

    transform(".bam") to("_macs2.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            touch $output;
            if [ ! -e ${macs2_vars.targets} ]; then
                echo "Targets file ${macs2_vars.targets} doesn't exist" >> $output &&
                exit 0;
            fi;

            BAM=\$(basename $input) &&
            grep \$BAM ${macs2_vars.targets} | while read -r TARGET; do
                IP=\$(       echo $TARGET | tr '\t' ' ' | cut -f1 -d" ") &&
                IPname=\$(   echo $TARGET | tr '\t' ' ' | cut -f2 -d" ") &&
                INPUT=\$(    echo $TARGET | tr '\t' ' ' | cut -f3 -d" ") &&
                INPUTname=\$(echo $TARGET | tr '\t' ' ' | cut -f4 -d" ");

                if [ "\$BAM" != "\$INPUT" ]; then
                    echo "\${IPname} vs \${INPUTname}" >> $output &&
                    macs2 callpeak -t ${macs2_vars.mapped}/\$IP -c ${macs2_vars.mapped}/\$INPUT -n \${IPname}.vs.\${INPUTname}_macs2 $MACS2_FLAGS &&
                    if [ \$? -ne 0 ]; then rm $output; fi &&
                    mv \${IPname}.vs.\${INPUTname}_macs2* $output.dir;
                fi;
            done
        ""","macs2"
    }

    forward input
}

