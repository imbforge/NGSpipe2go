GO_Enrichment = {

	doc title: "GO_Enrichment",
		desc: "Gene Ontology enrichment analysis",
		constraints: "",
		bpipe_version: "",
		author: ""

	output.dir = GO_Enrichment_OUTDIR.replaceFirst("out=", "")
	def GO_Enrichment_FLAGS = GO_Enrichment_RDATA.replaceFirst("out=", "") + "/DE_DESeq2.RData " + 
				GO_Enrichment_LOG2FOLD + " " + 
				GO_Enrichment_PVALUE + " " +
				GO_Enrichment_ORGDB + " " +
				GO_Enrichment_UNIV + " " +
				GO_Enrichment_TYPE + " " +
                GO_Enrichment_CATEGORY + " " +
				GO_Enrichment_OUTDIR + " " +
				GO_Enrichment_EXTRA

	transform(".RData") to("_GO.done") {
		exec """
		module load R &&

		touch $output; 
		echo 'VERSION INFO' 1>&2 ;
		echo \$(${TOOL_GO}/bin/Rscript --version 2>&1 | cut -d' '-f5) 1>&2 ;
		echo '/VERSION INFO'  1>&2 ;

			Rscript ${TOOL_GO}/GO_Enrichment.R $GO_Enrichment_FLAGS;
			if [ \$? -ne 0 ]; then rm $output; fi;
				
	""","GO_Enrichment"
	}

}
