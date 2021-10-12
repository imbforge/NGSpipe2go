rMATS = {
    doc title: "rMats Wrapper for alternative splicing events",
        desc:  "Differential expression analysis for alternative splicing events",
        constraints: "",
        author: "Nastasja Kreim, Modified by Sivarajan Karunanithi"

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

    def MASER_FLAGS = 
        (maser_vars.gtf      ? " gtf="      +   maser_vars.gtf      : "" ) +
        (maser_vars.db       ? " db="       +   maser_vars.db       : "" ) +
        (maser_vars.ftype    ? " ftype="    +   maser_vars.ftype    : "" ) +
        (maser_vars.mincov   ? " mincov="   +   maser_vars.mincov   : "" ) +
        (maser_vars.fdr      ? " fdr="      +   maser_vars.fdr      : "" ) + 
        (maser_vars.dpsi     ? " dpsi="     +   maser_vars.dpsi     : "" )

    def TOOL_ENV = prepare_tool_env("rmats", tools["rmats"]["version"], tools["rmats"]["runenv"]) + " && " +
                   prepare_tool_env("R", tools["R"]["version"], tools["R"]["runenv"])
    def PREAMBLE = get_preamble(stage:stageName, outdir:output.dir, input:new File(input1.prefix).getName())

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
            add_extensions="";
            if [[ "$RUN_CUTADAPT" == "true" ]]; then
                echo "Will add cutadapt extension";
                add_extensions=".cutadapt";
            fi;
            echo \${add_extensions};
            mkdir -p $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS &&
            bamgroup0=`awk -v g="\${groups[0]}" -v M="$MAPPED" -v P="${rMATS_vars.prefix}" -v C="\${add_extensions}" 'BEGIN{OFS=""} {if (\$3 == g) print M "/" P \$1 C ".bam"}' $input | paste -sd, -` &&
            echo \$bamgroup0 > $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS/\${groups[0]}_samples.txt &&
            bamgroup1=`awk -v g="\${groups[1]}" -v M="$MAPPED" -v P="${rMATS_vars.prefix}" -v C="\${add_extensions}" 'BEGIN{OFS=""} {if (\$3 == g) print M "/" P \$1 C ".bam"}' $input | paste -sd, -` &&
            echo \$bamgroup1 > $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS/\${groups[1]}_samples.txt &&
            run_rmats --b1 $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS/\${groups[0]}_samples.txt --b2 $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS/\${groups[1]}_samples.txt $RMATS_FLAGS --od $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS --tmp $output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS_tmp && 
            Rscript ${PIPELINE_ROOT}/tools/maser/createMaserPlots.R $MASER_FLAGS rmats_dir=$output.dir/\${input_var%${rMATS_vars.suffix}}_rMATS scripts_dir=${PIPELINE_ROOT}/tools/maser group1=\${groups[0]} group2=\${groups[1]} &&
            touch $output;
        ""","rMATS"
    }
}

