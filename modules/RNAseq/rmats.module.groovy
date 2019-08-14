// Notes:
//  * Indentation is important in this file. Please, use 4 spaces for indent. *NO TABS*.

load PIPELINE_ROOT + "/modules/RNAseq/rmats.vars.groovy"

rMATS = {
    doc title: "rMats Wrapper for alternative splicing events",
        desc:  "Differential expression analysis for alternative splicing events",
        constraints: "",
        author: "Nastasja Kreim"

    output.dir = RMATS_OUTDIR
    def RMATS_FLAGS = RMATS_GTF +
                      RMATS_LENGTH +
                      RMATS_THREADS + 
                      RMATS_EXTRA
    if(RMATS_PAIRED == "yes") {
        RMATS_FLAGS = RMATS_FLAGS + " -t paired"
    } else {
        RMATS_FLAGS = RMATS_FLAGS + " -t single"
    }
    if(RMATS_STRANDED == "no") {
        RMATS_FLAGS = "--libType fr-unstranded " + RMATS_FLAGS
    }
    else if (RMATS_STRANDED == "yes") {
        RMATS_FLAGS = "--libType fr-firststrand " + RMATS_FLAGS
    }
    else {
        RMATS_FLAGS = "--libType fr-secondstrand " + RMATS_FLAGS
    }

    def TOOL_ENV = prepare_tool_env("rmats", tools["rmats"]["version"], tools["rmats"]["runenv"])

    // run the chunk
    transform (".txt") to (".done") {
        exec """
            ${TOOL_ENV} &&
    
            if [ -n "\$SLURM_JOBID" ]; then
                export TMPDIR=/jobdir/\${SLURM_JOBID};
            fi;
            input_var=`basename $input`;
            postfix=$RMATS_POSTFIX;
            groups=\${input_var%\$postfix};
            echo "groups after postfix removal" \$groups;
            sep=$RMATS_SEP;
            groups=(\${groups//\$sep/ });
            echo "groups 0 " \${groups[0]};
            gcol=\$(head -n 1 $input | awk 'BEGIN{FS="\t"}{ for(fn=1; fn<=NF;fn++){ if(\$fn == "group") print fn;}}' ) &&
            fcol=\$(head -n 1 $input | awk 'BEGIN{FS="\t"}{ for(fn=1; fn<=NF;fn++){ if(\$fn == "file") print fn;}}' ) &&
            mkdir $output.dir/\${input_var%$PREMATS_POSTFIX}_rMATS &&
            bamgroup0=`awk -v g="\${groups[0]}" -v M="$MAPPED" -v fcol=\$fcol -v gcol=\$gcol 'BEGIN{OFS=""} {if (\$gcol == g) print M "/" \$fcol}' $input | paste -sd, -` &&
            echo \$bamgroup0 > $output.dir/\${input_var%$PREMATS_POSTFIX}_rMATS/\${groups[0]}_samples.txt &&
            bamgroup1=`awk -v g="\${groups[1]}" -v M="$MAPPED" -v fcol=\$fcol -v gcol=\$gcol 'BEGIN{ OFS=""} {if (\$gcol == g) print M "/" \$fcol}' $input | paste -sd, -` &&
            echo \$bamgroup1 > $output.dir/\${input_var%$PREMATS_POSTFIX}_rMATS/\${groups[1]}_samples.txt &&
            rmats.py --b1 $output.dir/\${input_var%$PREMATS_POSTFIX}_rMATS/\${groups[0]}_samples.txt --b2 $output.dir/\${input_var%$PREMATS_POSTFIX}_rMATS/\${groups[1]}_samples.txt $RMATS_FLAGS --od $output.dir/\${input_var%$PREMATS_POSTFIX}_rMATS &&
            touch $output;
        ""","rMATS"
    }
}

