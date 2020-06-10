DownsamplefastqSE = {
    doc title: "downsample",
        desc: "downsample wrapper for fastq files (single end)",
        constraints: "",
        bpipe_version: "tested with bpipe 0.9.8.7",
        author: "Nastasja Kreim"

    output.dir = DownsamplefastqSE_vars.outdir
    def OUTPUTFILES = new ArrayList()
    inputs.eachWithIndex { item, index -> 
        File f = new File(item)
        OUTPUTFILES.add((f.getName() =~ /.fastq.gz/).replaceFirst(".down.fastq.gz"))
        OUTPUTFILES.add(item)
        println OUTPUTFILES[index]
    }

    def PREAMBLE = get_preamble("DownsamplefastqSE")

    produce(OUTPUTFILES) {
        exec """
            paste <(zcat $input) | awk '{ printf("%s",\$0); n++; if(n%4==0) { printf("\\n");} else { printf("\\t\\t");} }' | shuf | head -n ${DownsamplefastqSE_vars.amount} | sed 's/\\t\\t/\\n/g' | awk -v r1=$output1.prefix 'BEGIN {FS="\\t"}{print \$1 >r1}' &&
            gzip $output1.prefix
    ""","Downsamplefastq"
    }
}

