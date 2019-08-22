filbowtie2unique_vars=[
    mapped          : MAPPED,
    samtools_mapq   : "3", // MAPQ >=3 should exclude "true multireads", multi mapped reads within the window of insert size
    samtools_threads: Integer.toString(ESSENTIAL_THREADS)
]
