//rule for task rMats from catalog miscellaneous, version 1
//desc: preform RMats analysis for the contrast target files as produced by PREMATS
rMATS = {
    doc title: "rMats Wrapper for alternative splicing evencts",
        desc:  "Differential expression analysis for alternative splicing events",
        constraints: "",
        author: "Nastasja Kreim"

    output.dir = RMATS_OUTDIR
    
    // run the chunk
    transform (".txt") to (".done") {
        exec """
            module load rmats/${RMATS_VERSION} &&
			if [ -n "\$SLURM_JOBID" ]; then
				export TMPDIR=/jobdir/\${SLURM_JOBID};
			fi;
				input_var=`basename $input`;
				postfix=$RMATS_POSTFIX
				groups=\${input_var%\$postfix};
				echo "groups after postfix removal" \$groups;
				sep=$RMATS_SEP
				groups=(\${groups//\$sep/ });
				echo "groups 0 " \${groups[0]};
				group0=`awk -v g="\${groups[0]}" -v M="$MAPPED" 'BEGIN{OFS=""} {if (\$2 == g) print M "/" \$1}' $input | paste -sd, -`;
				group1=`awk -v g="\${groups[1]}" -v M="$MAPPED" 'BEGIN{ OFS=""} {if (\$2 == g) print M "/" \$1}' $input | paste -sd, -`;
				PAIRED="single";
				PAIRED_ARGS="";
				if [ $RMATS_PAIRED == "yes" ]; then
					PAIRED="paired";
					insert_group0=`awk -v g="\${groups[0]}" 'BEGIN{OFS=""} {if (\$2 == g) print \$5}' $input | paste -sd"," -`;
					insert_group1=`awk -v g="\${groups[1]}" 'BEGIN{OFS=""} {if (\$2 == g) print \$5}' $input | paste -sd"," -`;
					sd_group0=`awk -v g="\${groups[0]}" 'BEGIN{OFS=""} {if (\$2 == g) print \$6}' $input | paste -sd"," -`;
					sd_group1=`awk -v g="\${groups[1]}" 'BEGIN{OFS=""} {if (\$2 == g) print \$6}' $input | paste -sd"," -`;
					PAIRED_ARGS="-r1 \$insert_group0 -r2 \$insert_group1 -sd1 \$sd_group0 -sd2 \$sd_group1";
				fi;
				python ${TOOL_RMATS}/RNASeq-MATS.py -b1 \$group0 -b2 \$group1 -t \$PAIRED -gtf $RMATS_GTF -a $RMATS_ANCHOR -c $RMATS_CUTOFF -analysis $RMATS_TYPE -len $RMATS_LEN  \$PAIRED_ARGS -o $output.dir/\${input_var%_targets_rMats.txt}_rMats;
				touch $output;
				echo "python ${TOOL_RMATS}/RNASeq-MATS.py -b1 \$group0 -b2 \$group1 -t \$PAIRED -gtf $RMATS_GTF -a $RMATS_ANCHOR -c $RMATS_CUTOFF -analysis $RMATS_TYPE -len $RMATS_LEN  \$PAIRED_ARGS -o $output.dir/\${input_var%_targets_rMats.txt}_rMats" > $output;
				""","rMATS"
        }
}

