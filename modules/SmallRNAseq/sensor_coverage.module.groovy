SensorCoverage = {
	doc title: "SensorCoverage",
		desc:  "Coverage at a small RNA sensor construct.",
		author: "António Domingues"

	output.dir = SENSOR_COVERAGE_OUTDIR

	transform(".bam") to(".plus.cov", ".minus.cov") {
		exec """

			module load bedtools/${BEDTOOLS_VERSION} &&

			genomeCoverageBed -ibam $input -d -strand "+" | grep "21U_mCherry_sensor" > $output1 &&
			genomeCoverageBed -ibam $input -d -strand "-" | grep "21U_mCherry_sensor" > $output2

		""","SensorCoverage"
	}
}

