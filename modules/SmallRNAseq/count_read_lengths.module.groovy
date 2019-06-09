CountReadLengths = {
    doc title: "CountReadLengths",
    desc: "Determines the sequence length distribution in a fastq file",
    constraints: "Gzipped fastq",
    author: "António Domingues"

    output.dir = READ_LENGTH_DIR

    transform(".fastq.gz") to (".readlength.txt") {
        exec """

            zcat $input | awk '{if(NR%4==2) print length(\$1)}' | sort -n | uniq -c > $output
        ""","CountReadLengths"
   }
   forward input
}
