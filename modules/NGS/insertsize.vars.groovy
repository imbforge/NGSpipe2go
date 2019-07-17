//vars for task MarkDups from catalog NGS, version 1
INSERTSIZE_MAXMEM="-Xmx5000m" //set the java heap size
INSERTSIZE_OUTDIR= QC + "/insertsize" //location of the OUTPUT Dir
INSERTSIZE_OTHER="ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT" //sometimes the sorted flag is not set and we should not care if we have reads which overhang chromosomes
