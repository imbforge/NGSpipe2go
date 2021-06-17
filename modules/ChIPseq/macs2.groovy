macs2 = {
    doc title: "MACS2",
        desc:  "MACS2 wrapper",
        constraints: "Performs treatment versus control peak calling. If no input control sample is available, the INPUT field in the targets.txt file must be given as none.",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Sergi Sayols, Frank Rühle"

    var subdir : ""
    output.dir = macs2_vars.outdir + "/$subdir" 

    def MACS2_FLAGS =
        (macs2_vars.gsize  ? " -g "           + macs2_vars.gsize  : "") +
        (macs2_vars.minlen ? " --min-length " + macs2_vars.minlen : "") +
        (macs2_vars.broad  ? " --broad "                          : "") +
        (macs2_vars.paired ? " --format BAMPE"                    : "") +
        (macs2_vars.extra  ? " " + macs2_vars.extra               : "")

    def TOOL_ENV = prepare_tool_env("macs2", tools["macs2"]["version"], tools["macs2"]["runenv"])
    def PREAMBLE = get_preamble(module:"macs2", branch:branch, branch_outdir:subdir)

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
            extension="\${BAM#*.}" &&
            BAM="\${BAM%%.*}" &&
            grep "^\$BAM" ${macs2_vars.targets} | while read -r TARGET; do
                IP=\$(       echo $TARGET | tr '\t' ' ' | cut -f1 -d" ").\$extension &&
                IPname=\$(   echo $TARGET | tr '\t' ' ' | cut -f2 -d" ") &&
                INPUT=\$(    echo $TARGET | tr '\t' ' ' | cut -f3 -d" ").\$extension &&
                INPUTname=\$(echo $TARGET | tr '\t' ' ' | cut -f4 -d" ");

                if [ "\$BAM" != "\$INPUT" ]; then
                    if [ "\$INPUT" != "none.\$extension" ]; then
                        NAMEoutput="\${IPname}.vs.\${INPUTname}" &&
                        INPUTflag="-c ${macs2_vars.mapped}/\$INPUT";
                    else
                        NAMEoutput="\${IPname}" &&
                        INPUTflag="";
                    fi &&   
                    echo "\$NAMEoutput" >> $output &&
                    macs2 callpeak -t ${macs2_vars.mapped}/\$IP \$INPUTflag -n $subdir\${NAMEoutput}_macs2 $MACS2_FLAGS &&
                    if [ \$? -ne 0 ]; then rm $output; fi &&
                    find . -maxdepth 1 -name "$subdir\${NAMEoutput}_macs2*" -exec sh -c 'mv "\$1" "$output.dir/\${1#./$subdir}"' _ {} \\;;
                fi;
            done
        ""","macs2"
    }

    forward input
}

