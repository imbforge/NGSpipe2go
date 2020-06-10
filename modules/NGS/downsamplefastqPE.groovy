DownsamplefastqPE = {
    doc title: "downsample",
        desc:    "downsample wrapper for fastq files (paired end)",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = DownsamplefastqPE_vars.outdir
    def OUTPUTFILES = new ArrayList()
    inputs.eachWithIndex { item, index -> 
        File f = new File(item)
        OUTPUTFILES.add((f.getName() =~ /.fastq.gz/).replaceFirst(".down.fastq.gz"))
        println OUTPUTFILES[index]
    }

    def PREAMBLE = get_preamble("DownsamplefastqPE")

    produce(OUTPUTFILES){
        exec """
            paste <(zcat $input1) <(zcat $input2) | awk '{ printf("%s",\$0); n++; if(n%4==0) { printf("\\n");} else { printf("\\t\\t");} }' | shuf    | head -n ${DownsamplefastqPE_vars.amount} | sed 's/\\t\\t/\\n/g' | awk -v r1=$output1.prefix -v r2=$output2.prefix 'BEGIN {FS="\\t"}{print \$1 >r1; print \$2>r2 }' &&
            gzip $output1.prefix &&
            gzip $output2.prefix;
        ""","Downsamplefastq"
    }
}

