MirDeep2 = {
	doc title: "mirDeep2",
		desc:  """Quantification of miRNAs performed in 2 steps: (1) Processes reads and mappping to the reference genome; (2) quantification of miRNA expression. This is step 2""",
      constraints: "Requires mirDeep2.",
		author: "Antonio Domingues"

	output.dir = MIR_OUTDIR

	transform(".arf", ".fa") to (".csv", ".html") {


      exec """
         if [ -n "\$LSB_JOBID" ]; then
            export TMPDIR=/jobdir/\${LSB_JOBID};
         fi                                          &&

         export PATH=${TOOL_DEPENDENCIES}:$PATH &&
         export PATH=${TOOL_VIENNA}:$PATH &&
         export PATH=${TOOL_RANDFOLD}:$PATH &&
         export PERL5LIB=$PERL5LIB:/fsimb/groups/imb-kettinggr/common_bin/mirdeep2_0_0_7/lib/ &&
         export PATH=${TOOL_MIRDEEP2}:$PATH &&

         ${TOOL_MIRDEEP2}/miRDeep2.pl $input2 $GENOME_SEQ $input1 $MATURE_MIRNA none $HAIRPIN_MIRNA -t zebrafish -c -r $input1.prefix 2> $output.dir/report.log

		""","MirDeep2"
	}
}
