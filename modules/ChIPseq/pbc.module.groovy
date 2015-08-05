//rule for task pbc from catalog ChIPseq, version 1
//desc: PCR Bottleneck Coefficient
pbc = {
	doc title: "PBC",
		desc:  "PCR Bottleneck Coefficient",
		constraints: "",
		author: "Sergi Sayols"

	output.dir = QC + "/pbc"

	transform(".bam") to("_PBC.csv") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES  &&
			source ${TOOL_R}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi;

			${TOOL_R}/bin/Rscript ${TOOL_DEPENDENCIES}/imb-forge/encodeChIPqc/PBC.R $input && mv ${input.prefix}_PBC.csv $output.dir
		""","pbc"
	}
	
	forward input

}