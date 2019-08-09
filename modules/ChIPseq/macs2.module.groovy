// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/macs2.vars.groovy"

macs2 = {
    doc title: "MACS2",
        desc:  "MACS2 wrapper",
        constraints: "Only performs treatment versus control peak calling",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = RESULTS + "/macs2"

    def MACS2_FLAGS= MACS2_GSIZE  + " " +
                     MACS2_BWIDTH + " " +
                     MACS2_MINLEN + " " +
                     MACS2_EXTRA
    if(MACS2_PAIRED == "yes") {
        MACS2_FLAGS = MACS2_FLAGS + " " + "--format BAMPE"
    }

    def TOOL_ENV = prepare_tool_env("macs2", tools["macs2"]["version"], tools["macs2"]["runenv"])

    transform(".bam") to("_macs2.done") {
        exec """
            ${TOOL_ENV} &&

            touch $output;
            if [ ! -e $MACS2_TARGETS ]; then
                echo "Targets file $MACS2_TARGETS doesn't exist" >> $output &&
                exit 0;
            fi;

            BAM=\$(basename $input) &&
            grep \$BAM $MACS2_TARGETS | while read -r TARGET; do
                IP=\$(       echo $TARGET | tr '\t' ' ' | cut -f1 -d" ") &&
                IPname=\$(   echo $TARGET | tr '\t' ' ' | cut -f2 -d" ") &&
                INPUT=\$(    echo $TARGET | tr '\t' ' ' | cut -f3 -d" ") &&
                INPUTname=\$(echo $TARGET | tr '\t' ' ' | cut -f4 -d" ");

                if [ "\$BAM" != "\$INPUT" ]; then
                    echo "\${IPname} vs \${INPUTname}" >> $output &&
                    macs2 callpeak -t $MACS2_MAPPED/\$IP -c $MACS2_MAPPED/\$INPUT -n \${IPname}.vs.\${INPUTname}_macs2 $MACS2_FLAGS &&
                    if [ \$? -ne 0 ]; then rm $output; fi &&
                    mv \${IPname}.vs.\${INPUTname}_macs2* $output.dir;
                fi;
            done
        ""","macs2"
    }

    forward input
}

