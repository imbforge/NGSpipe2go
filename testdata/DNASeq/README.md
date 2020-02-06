Run this bash code to generate 5M 101bp of PE reads for a reference (eg. Yeast), 2 replicates x 2 conditions.

```sh
#!/bin/bash
# 
set -euo pipefail

export READ_LEN=101
export FRAG_LEN=500
export NUM_READS=5000000
export ERR_RATE=0.001
export MUT_RATE=0.0001
export INDEL_RATE=0.15
export INDEL_EXTEND_RATE=0.3

export REF="/fsimb/groups/imb-bioinfocf/common-data/GATK_resources/gatk_bundle_hg38_v0/Homo_sapiens_assembly38.fasta"
export RAWDATA="./rawdata"

CORES=4

function f {
  REF=$1
  BASE=$2
  REPL=${RAWDATA}/${BASE}_$3
  SEED=$3

  echo "replicate $REPL using ref $REF"
  wgsim \
    -1${READ_LEN} \
    -2${READ_LEN} \
    -d${FRAG_LEN} \
    -N${NUM_READS} \
    -e${ERR_RATE} \
    -r${MUT_RATE} \
    -R${INDEL_RATE} \
    -X${INDEL_EXTEND_RATE} \
    -S${SEED} \
    ${REF} ${REPL}.R1.fastq ${REPL}.R2.fastq | gzip > ${REPL}_sim.txt.gz

  gzip ${REPL}.R1.fastq ${REPL}.R2.fastq
}
export -f f

# Generate 2 replicates for 2 conditions. Both from the same reference, thus they'll have no differences.
# Use a modified reference for treated, if you wanna introduce changes
parallel --xapply -j $CORES "f {1} {2} {3}" ::: "$REF" ::: control treated ::: 1 1 2 2
```
