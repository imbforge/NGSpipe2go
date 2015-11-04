//rule for task BWA mem from DNAseq pipeline
//desc: Align single and paired end reads
BWA_pe = {
	doc title: "BWA PE alignment",
	desc:  "Align paired end reads",
	constraints: "Only works with compressed input. Set all global vars.",
	bpipe_version: "tested with bpipe 0.9.8.7",
	author: "Oliver Drechsel"

	// check Simon Sadedin solution
	//	Hi Oliver,
	//
	//How I'd suggest to do it is by using "transform" - but, with two arguments:
	//
	//transform(".read_1.fastq.gz",".read_2.fastq.gz") to(".bam") {
	//    exec """
	//        cat $input1 $input2 > $output.bam
	//    """
	//    }
	//
	//
	//Make sure you update to the latest (beta) version as I remember there may
	//have been a bug fix related to something around this.
	//
	//Cheers / hope this helps!
	//
	//Simon
	//
	//On 1/09/2015 2:43 am, "Oliver Drechsel" <o.drechsel@imb-mainz.de> wrote:
	//
	//>Dear Simon,
	//>
	//>out group is happy user of bpipe.
	//>Currently, i try to build up a new pipeline, which requires paired input
	//>data.
	//>
	//>The input file names will look like this:
	//>NA12878.read_1.fastq.gz
	//>NA12878.read_2.fastq.gz
	//>
	//>
	//>My current usage looks like this:
	//>
	//>transform('.fastq.gz') to('.bam') {
	//>bwa mem ... $input1 $input2 > $ouput
	//>}
	//>
	//>Unfortunately, the output will be formatted like NA12878.read_1.bam, but
	//>my aim would be NA12878.bam for the downstream steps.
	//>
	//>I tried several options to reach this file name to no avail:
	//># defining several produces and transforms
	//> transform('.read*.fastq.gz') to('.bam') {...}
	//> transform('.read%.fastq.gz') to('.bam') {...}
	//> transform(.bam) {...}
	//> produce(.bam) {...} # also using various * and % constructs
	//>
	//># trying to modify $output
	//> bwa mem ... > ${output($output.prefix.prefix)}
	//>
	//># playing around with
	//> ${input.prefix} a lot
	//>
	//>All i tried now, didn't yield my expected output name. Neither did the
	//>example for paired input bwa module.
	//>
	//>Could you please point me to my misunderstanding and mistake, so the
	//>pipeline could feature the naming scheme i'd like to use?
	
	transform(".read_1.fastq.gz",".read_2.fastq.gz") to(".bam") {
		
        exec """
        
            echo "input1:"$input1 &&
			echo "input2:"$input2 &&
			echo "output:"$output &&
			echo "output1:"$output1 &&
			echo "output2:"$output2 &&
			
			
        ""","BWA_pe"
    }
}
