//rule for task qualimap from catalog NGS, version 1
//desc: get the amount of exonic intronic and intergenic reads to test for DNA contamination 
qualimap = {
	doc title: "Qualimap",
		desc:  "Call qualimap to do rnaseq qualitycontrol",
		author: "NastasjaKreim"

	output.dir = QUALIMAP_OUTDIR
	// no|yes|reverse
	if(QUALIMAP_STRANDED == "no") {
		QUALIMAP_PROTOCOL = "non-strand-specific"
	}
	else if (QUALIMAP_STRANDED == "yes") {
		QUALIMAP_PROTOCOL = "strand-specific-forward"
	}
	else {
		QUALIMAP_PROTOCOL = "strand-specific-reverse"
	}
	if(QUALIMAP_PAIRED == "yes"){
		QUALIMAP_EXTRA = QUALIMAP_EXTRA + " -pe"
	}

	transform(".bam") to("_counts.txt") {
		exec """
			module load qualimap/${QUALIMAP_VERSION} &&
			unset DISPLAY;
			echo $output.prefix;
			qualimap rnaseq -bam $input -outdir ${output.prefix}_qualimap -outformat html $QUALIMAP_GENESGTF -oc $output -p $QUALIMAP_PROTOCOL $QUALIMAP_EXTRA
		""","qualimap"
	}

	forward input
}

