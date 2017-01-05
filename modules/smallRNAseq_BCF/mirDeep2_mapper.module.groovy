miRDeep2Mapper = {
	doc title: "miRDeep2",
	desc:  """Quantification of miRNAs performed in 2 steps: (1) Processes reads and mappping to the reference genome; (2) quantification of miRNA expression.""",
      	constraints: "Requires mirDeep2.",
	author: "Antonio Domingues, Anke Busch"

	output.dir = MIR_MAPPER_OUTDIR

	transform(".fastq.gz") to (".arf", ".fa") {

		def SAMPLENAME = input.prefix
		def OUTPUTLOG_MAIN = output2.prefix

		exec """
			export TOOL_DEPENDENCIES=$TOOL_DEPENDENCIES &&
			source ${TOOL_MIRDEEP2}/env.sh &&

			if [ ! -d ${TMP} ]; then
                                mkdir -p ${TMP};
                        fi &&

			SAMPLENAME_BASE=\$(basename ${SAMPLENAME}) &&
			gunzip -c $input > ${TMP}/\${SAMPLENAME_BASE} &&

			cd $output.dir &&

			mapper.pl ${TMP}/\${SAMPLENAME_BASE} -e -p $GENOME_REF -s $output2 -t $output1 -h -m -i -j -o 8 &> ${OUTPUTLOG_MAIN}.mapper.log &&

			rm ${TMP}/\${SAMPLENAME_BASE}
			
		""","miRDeep2Mapper"
	}
}
