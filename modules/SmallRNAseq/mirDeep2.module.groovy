MirDeep2 = {
	doc title: "mirDeep2",
		desc:  """Quantification of miRNAs performed in 2 steps: (1) Processes reads and mappping to the reference genome; (2) quantification of miRNA expression.""",
      constraints: "Requires mirDeep2.",
		author: "Antonio Domingues"

	output.dir = MIR_OUTDIR

	transform(".deduped_barcoded.trimmed.fastq") to (".arf", ".fa") {


      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi                                          &&

         export PATH=${TOOL_DEPENDENCIES}:$PATH &&
         export PATH=${TOOL_MIRDEEP2}:$PATH &&

         ${TOOL_MIRDEEP2}/mapper.pl $input -e -p $GENOME_REF -s $output2.fa -t $output1 -h -m -i -j -o 8 1>&2 ${output2.prefix}.mapper.log

		""","MirDeep2"
	}
}
