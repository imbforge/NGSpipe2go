// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/rmats.vars.groovy"

rMATS = {
    doc title: "rMats Wrapper for alternative splicing events",
        desc:  "Differential expression analysis for alternative splicing events",
        constraints: "",
        author: "Nastasja Kreim"

    output.dir = rMATS_vars.outdir
    def RMATS_FLAGS =
        (rMATS_vars.gtf     ? " --gtf "        +  rMATS_vars.gtf     : "" ) +
        (rMATS_vars.length  ? " --readLength " +  rMATS_vars.length  : "" ) +
        (rMATS_vars.threads ? " --nthread "    +  rMATS_vars.threads : "" ) + 
        (rMATS_vars.extra   ? " "              +  rMATS_vars.extra   : "" ) 

    if(rMATS_vars.paired) {
        RMATS_FLAGS = RMATS_FLAGS + " -t paired"
    } else {
        RMATS_FLAGS = RMATS_FLAGS + " -t single"
    }
    if(rMATS_vars.stranded == "no") {
        RMATS_FLAGS = "--libType fr-unstranded " + RMATS_FLAGS
    }
    else if (rMATS_vars.stranded == "yes") {
        RMATS_FLAGS = "--libType fr-firststrand " + RMATS_FLAGS
    }
    else {
        RMATS_FLAGS = "--libType fr-secondstrand " + RMATS_FLAGS
    }

    def TOOL_ENV = prepare_tool_env("rmats", tools["rmats"]["version"], tools["rmats"]["runenv"])
    def PREAMBLE = get_preamble("rMATS")

    // run the chunk
    transform (".txt") to (".done") {
        exec """
            ${TOOL_ENV} &&
            ${PREAMBLE} &&
    
            input_var=`basename $input`;
            suffix=${rMATS_vars.suffix};
            groups=\${input_var%\$suffix};
            echo "groups after suffix removal" \$groups;
            sep=${rMATS_vars.sep};
            groups=(\${groups//\$sep/ });
            echo "groups 0 " \${groups[0]};
            gcol=\$(head -n 1 $input | awk 'BEGIN{FS="\t"}{ for(fn=1; fn<=NF;fn++){ if(\$fn == "group") print fn;}}' ) &&
            fcol=\$(head -n 1 $input | awk 'BEGIN{FS="\t"}{ for(fn=1; fn<=NF;fn++){ if(\$fn == "file") print fn;}}' ) &&
            mkdir $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS &&
            bamgroup0=`awk -v g="\${groups[0]}" -v M="$MAPPED" -v fcol=\$fcol -v gcol=\$gcol 'BEGIN{OFS=""} {if (\$gcol == g) print M "/" \$fcol}' $input | paste -sd, -` &&
            echo \$bamgroup0 > $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS/\${groups[0]}_samples.txt &&
            bamgroup1=`awk -v g="\${groups[1]}" -v M="$MAPPED" -v fcol=\$fcol -v gcol=\$gcol 'BEGIN{ OFS=""} {if (\$gcol == g) print M "/" \$fcol}' $input | paste -sd, -` &&
            echo \$bamgroup1 > $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS/\${groups[1]}_samples.txt &&
            rmats.py --b1 $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS/\${groups[0]}_samples.txt --b2 $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS/\${groups[1]}_samples.txt $RMATS_FLAGS --od $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS &&
            touch $output;
        ""","rMATS"
    }
}

