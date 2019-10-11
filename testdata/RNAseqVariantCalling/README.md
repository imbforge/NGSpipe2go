To simulate paired end RNA-seq libs for variant calling, first we'll generate a fake reference transcriptome
with the sequence (exonic+intronic) of 5 genes from chr19.

This part is done with tis R script. Dependencies are the `BSgenome.Hsapiens.UCSC.hg38` and the `ShortRead` packages:

```R
#############################
##
## Generate a fasta file with sequences of the following genes:
## Coordinates where retrieved from Biomart
##
## Gene stable ID	Gene start (bp)	Gene end (bp)	Chromosome/scaffold name	Strand	Gene name
## ENSG00000141837	13206442	13633025	19	-1	CACNA1A
## ENSG00000074181	15159038	15200995	19	-1	NOTCH3
## ENSG00000105679	35533455	35545319	19	1	GAPDHS
## ENSG00000130203	44905791	44909393	19	1	APOE
## ENSG00000121410	58345178	58353499	19	-1	A1BG
##
#############################
library(BSgenome.Hsapiens.UCSC.hg38)
library(ShortRead)

# retrieve sequences
genes <- with(data.frame(seqnames=c("chr19", "chr19", "chr19", "chr19", "chr19"),
                         start   =c(13206442, 15159038, 58345178, 44905791, 35533455),
                         end     =c(13633025, 15200995, 58353499, 44909393, 35545319),
                         strand  =strand(c("-", "-", "-", "+", "+")),
                         name    =c("CACNA1A", "NOTCH3", "A1BG", "APOE", "GAPDHS")),
              GRanges(seqnames, IRanges(start, end), strand, name))
seqs <- as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, genes))
names(seqs) <- genes$name

# save fasta
dir.create("ref")
ShortRead::writeFasta(DNAStringSet(seqs), "ref/ref.fa")
system("samtools faidx ref/ref.fa")
system("java -jar $picard CreateSequenceDictionary R=ref/ref.fa O=ref/ref.dict")
system("STAR --runMode genomeGenerate --genomeDir ref --genomeFastaFiles ref/ref.fa")
```

The second part, simulating reads from this reference, is done using samtools' `wgsim`:

```sh
#!/bin/bash
set -euo pipefail

export READ_LEN=101
export FRAG_LEN=500
export NUM_READS=1000000
export ERR_RATE=0.001
export MUT_RATE=0.0001
export INDEL_RATE=0.15
export INDEL_EXTEND_RATE=0.3

export REF="./ref/ref.fa"
export RAWDATA="./rawdata"

CORES=3

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
    ${REF} ${REPL}_R1.fastq ${REPL}_R2.fastq | gzip > ${REPL}_sim.txt.gz

  gzip ${REPL}_R1.fastq ${REPL}_R2.fastq
}
export -f f

parallel --xapply -j $CORES "f {1} {2} {3}" ::: "$REF" ::: sample ::: 1 2 3
```

The VCF with know sites was extracted from the gatk_bundle_b37, and only those variants within the coordinates of the 5 genes where selected.
Then, coordinates were made relative to the start of the gene (with R).

