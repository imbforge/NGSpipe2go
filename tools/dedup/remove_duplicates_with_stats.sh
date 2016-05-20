#!/usr/bin/env bash

input=${1}
output=${input/.highQ.fastq.gz/.deduped_barcoded.fastq.gz}
highQ=${input%.highQ.fastq.gz}".highQ"
unique=${input%.highQ.fastq.gz}".unique"

zcat ${input} | paste -d, - - - - | tee >(awk -v var="$highQ" 'END {print NR,var}' >> dedup.stats.txt) | sort -u -t, -k2,2 | tee >(awk -v var="$uniq" 'END {print NR,var}' >> dedup.stats.txt) | tr ',' '\\n' | gzip > ${output}
