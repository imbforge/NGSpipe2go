// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/ChIPseq/ipstrength.vars.groovy"

ipstrength = {
    doc title: "IPstrength plot",
        desc:  "IPstrength",
        constraints: "install the right reference BSgenome",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols"

    output.dir = QC + "/ipstrength"

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])

    transform(".bam") to("_ipstrength.done") {
        exec """
            ${TOOL_ENV} &&

            touch $output;
            if [ ! -e $IPSTRENGTH_TARGETS ]; then
                echo "Targets file $IPSTRENGTH_TARGETS doesn't exist" >> $output &&
                exit 0;
            fi;

            BAM=\$(basename $input) &&
            grep \$BAM $IPSTRENGTH_TARGETS | while read -r TARGET; do
                IP=\$(       echo $TARGET | tr '\t' ' ' | cut -f1 -d" ") &&
                IPname=\$(   echo $TARGET | tr '\t' ' ' | cut -f2 -d" ") &&
                INPUT=\$(    echo $TARGET | tr '\t' ' ' | cut -f3 -d" ") &&
                INPUTname=\$(echo $TARGET | tr '\t' ' ' | cut -f4 -d" ");

                if [ "\$BAM" != "\$INPUT" ]; then
                    echo "\${IPname} vs \${INPUTname}" >> $output ;
                    Rscript ${PIPELINE_ROOT}/tools/ENCODEqc/IPstrength.R $IPSTRENGTH_MAPPED/\$IP \$IPname $IPSTRENGTH_MAPPED/\$INPUT \$INPUTname \${IPname}.vs.\${INPUTname}_ipstrength $IPSTRENGTH_BSGENOME;
                    if [ \$? -ne 0 ]; then rm $output; fi;
                    mv \${IPname}.vs.\${INPUTname}_ipstrength* $output.dir;
                fi;
            done
        ""","ipstrength"
    }

    forward input
}

