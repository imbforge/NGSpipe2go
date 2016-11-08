TrimUMIs = {
	doc title: "Trim UMIs",
	desc:  """Trims random barcodes that help in the identification of PCR duplicates and are in adapter-removed reads: NNNN-insert-NNNN.""",
      	constraints: "Requires seqtk.",
	author: "Antonio Domingues, Anke Busch"

	output.dir = TRIM_OUTDIR
	def SEQTK_FLAGS = 	" -l " + LEFT_TRIM + 
				" -b " + RIGHT_TRIM

	transform(".deduped.fastq.gz") to (".deduped.trimmed.fastq.gz") {

		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_SEQTK}/env.sh &&

		        echo 'VERSION INFO'  1>&2 &&
		        echo \$(seqtk 2>&1 | grep 'Version' | cut -d' ' -f2) 1>&2 &&
		        echo '/VERSION INFO' 1>&2 &&

		        seqtk trimfq -b ${LEFT_TRIM} -e ${RIGHT_TRIM} $input | gzip > $output

		""","TrimUMIs"
	}
}
