//rule for task qualimap from catalog NGS, version 1
//desc: get the amount of exonic intronic and intergenic reads to test for DNA contamination 
qualimap = {
	doc title: "Qualimap",
		desc:  "Call qualimap to do rnaseq qualitycontrol",
		author: "NastasjaKreim"

	output.dir = QUALIMAP_OUTDIR
	// no|yes|reverse
	if(QUALIMAP_STRANDED == "no") {
		QUALIMAP_STRANDED = "non-strand-specific"
	}
	else if (QUALIMAP_STRANDED == "yes") {
		QUALIMAP_STRANDED = "strand-specific-forward"
	}
	else {
		QUALIMAP_STRANDED = "strand-specific-reverse"
	}
	if(QUALIMAP_PAIRED == "yes"){
		QUALIMAP_EXTRA = QUALIMAP_EXTRA + " -pe"
	}

	transform(".bam") to("_counts.txt") {
		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_QUALIMAP}/env.sh &&
			if [ -n "\$LSB_JOBID" ]; then
				export TMPDIR=/jobdir/\${LSB_JOBID};
			fi &&
			unset DISPLAY;
			echo $output.prefix;
			qualimap rnaseq -bam $input -outdir ${output.prefix}_qualimap -outformat html $QUALIMAP_GENESGTF -oc $output -p $QUALIMAP_STRANDED $QUALIMAP_EXTRA
		""","qualimap"
	}

	forward input
}

