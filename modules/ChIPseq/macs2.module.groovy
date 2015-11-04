macs2 = {
	doc title: "MACS2",
		desc:  "MACS2 wrapper",
		constraints: "Only performs treatment versus control peakcalling",
		bpipe_version: "tested with bpipe 0.9.8.7",
		author: "Sergi Sayols"

	output.dir = RESULTS + "/macs2"
	MACS2_FLAGS= " -m "   + MACS2_MFOLD  +
	             " -g "   + MACS2_GSIZE  +
				 " --bw " + MACS2_BWIDTH +
				 MACS2_OTHER

	transform(".bam") to("_macs2.done") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES  &&
			source ${TOOL_MACS2}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi;

			touch $output;
			if [ ! -e $MACS2_TARGETS ]; then
				echo "Targets file $MACS2_TARGETS doesn't exist" >> $output &&
				exit 0;
			fi;
			
			echo 'VERSION INFO'  1>&2 &&
			${TOOL_MACS2}/bin/macs2 --version &&
			echo '/VERSION INFO' 1>&2 &&
			
			BAM=\$(basename $input) &&
			grep \$BAM $MACS2_TARGETS | while read -r TARGET; do
				IP=\$(       echo \$TARGET | cut -f1 -d" ") &&
				IPname=\$(   echo \$TARGET | cut -f2 -d" ") &&
				INPUT=\$(    echo \$TARGET | cut -f3 -d" ") &&
				INPUTname=\$(echo \$TARGET | cut -f4 -d" ");

				if [ "\$BAM" != "\$INPUT" ]; then
					echo "\${IPname} vs \${INPUTname}" >> $output ;
					${TOOL_MACS2}/bin/macs2 callpeak -t $MACS2_MAPPED/\$IP -c $MACS2_MAPPED/\$INPUT -n \${IPname}.vs.\${INPUTname}_macs2 $MACS2_FLAGS;
					if [ \$? -ne 0 ]; then rm $output; fi;
					mv \${IPname}.vs.\${INPUTname}_macs2* $output.dir;
				fi;
			done
		""","macs2"
	}

	forward input
}

