//rule for task ipstrength from catalog ChIPseq, version 1
//desc: IPstrength
ipstrength = {
	doc title: "IPstrength plot",
		desc:  "IPstrength",
		constraints: "install the right reference BSgenome",
		author: "Sergi Sayols"

	output.dir = QC + "/ipstrength"

	transform(".bam") to("_ipstrength.done") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES  &&
			source ${TOOL_R}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi;

			touch $output;
			if [ ! -e $IPSTRENGTH_TARGETS ]; then
				echo "Targets file $IPSTRENGTH_TARGETS doesn't exist" >> $output &&
				exit 0;
			fi;
			
			echo 'VERSION INFO'  1>&2 &&
			echo \$(${TOOL_R}/bin/Rscript --version) 1>&2 &&
			echo '/VERSION INFO' 1>&2 &&
			
			
			BAM=\$(basename $input) &&
			grep \$BAM $IPSTRENGTH_TARGETS | while read -r TARGET; do
				IP=\$(       echo $TARGET | cut -f1 -d" ") &&
				IPname=\$(   echo $TARGET | cut -f2 -d" ") &&
				INPUT=\$(    echo $TARGET | cut -f3 -d" ") &&
				INPUTname=\$(echo $TARGET | cut -f4 -d" ");

				if [ "\$BAM" != "\$INPUT" ]; then
					echo "\${IPname} vs \${INPUTname}" >> $output ;
					${TOOL_R}/bin/Rscript ${TOOL_ENCODEqc}/IPstrength.R $IPSTRENGTH_MAPPED/\$IP \$IPname $IPSTRENGTH_MAPPED/\$INPUT \$INPUTname \${IPname}.vs.\${INPUTname}_ipstrength $IPSTRENGTH_BSGENOME;
					if [ \$? -ne 0 ]; then rm $output; fi;
					mv \${IPname}.vs.\${INPUTname}_ipstrength* $output.dir;
				fi;
			done
		""","ipstrength"
	}

	forward input
}

