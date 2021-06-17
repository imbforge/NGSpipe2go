ipstrength = {
    doc title: "IPstrength plot",
        desc:  "IPstrength",
        constraints: "install the right reference BSgenome",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Frank Rühle"

    var subdir : ""
    output.dir = ipstrength_vars.outdir + "/$subdir" 

    def TOOL_ENV = prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, subdir:subdir, input:new File(input1.prefix).getName())

    transform(".bam") to("_ipstrength.done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&

            touch $output;
            if [ ! -e ${ipstrength_vars.targets} ]; then
                echo "Targets file ${ipstrength_vars.targets} doesn't exist" >> $output &&
                exit 0;
            fi;

            BAM=\$(basename $input) &&           
            extension="\${BAM#*.}" &&
            BAM="\${BAM%%.*}" &&
            grep "^\$BAM" ${ipstrength_vars.targets} | while read -r TARGET; do
                IP=\$(       echo $TARGET | tr '\t' ' ' | cut -f1 -d" ").\$extension &&
                IPname=\$(   echo $TARGET | tr '\t' ' ' | cut -f2 -d" ") &&
                INPUT=\$(    echo $TARGET | tr '\t' ' ' | cut -f3 -d" ").\$extension &&
                INPUTname=\$(echo $TARGET | tr '\t' ' ' | cut -f4 -d" ");

                if [ "\$BAM" != "\$INPUT" ]; then
                    echo "\${IPname} vs \${INPUTname}" >> $output ;
                    Rscript ${PIPELINE_ROOT}/tools/ENCODEqc/IPstrength.R ${ipstrength_vars.mapped}/\$IP \$IPname ${ipstrength_vars.mapped}/\$INPUT \$INPUTname $subdir\${IPname}.vs.\${INPUTname}_ipstrength ${ipstrength_vars.bsgenome};
                    if [ \$? -ne 0 ]; then rm $output; fi;
                    find . -maxdepth 1 -name "$subdir\${IPname}.vs.\${INPUTname}_ipstrength*" -exec sh -c 'mv "\$1" "$output.dir/\${1#./$subdir}"' _ {} \\;;
                fi;
            done
        ""","ipstrength"
    }

    forward input
}

