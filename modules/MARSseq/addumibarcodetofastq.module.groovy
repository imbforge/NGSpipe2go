AddUMIBarcodeToFastq = {
	doc title: "Adds UMI and Barcode of to the fastq header",
		desc:  "adds UMI and barcode of the second read in MARS-Seq samples to the fastq header using umitools",
		constraints: "",
		bpipe_version: "tested with bpipe 0.9.9.2",
		author: "Nastasja Kreim"

		output.dir = ADDUMIBARCODE_OUTDIR
    def OUTPUTFILE = input1
    int path_index = OUTPUTFILE.lastIndexOf("/")
    OUTPUTFILE = OUTPUTFILE.substring(path_index+1)
    OUTPUTFILE = (OUTPUTFILE =~ /.R1.fastq.gz/).replaceFirst("")
    def ADDUMIBARCODE_FLAGS = ADDUMIBARCODE_BCPATTERN + 
                              " " +  ADDUMIBARCODE_BARCODELIST +
                              " " +  ADDUMIBARCODE_EXTRA

  produce(OUTPUTFILE + ".umibarcode.fastq.gz"){
    exec """
     module load umitools/${UMITOOLS_VERSION} &&
     umi_tools extract $ADDUMIBARCODE_FLAGS -I $input2 --stdout ${input2.prefix}.barcode.fastq.gz --read2-in $input1 --read2-out=$output &&
     rm ${input2.prefix}.barcode.fastq.gz
  ""","AddUMIBarcodeToFastq"
}
}

