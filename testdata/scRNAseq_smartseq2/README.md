To simulate paired end RNA-seq libs, first we'll generate a fake reference transcriptome
with the sequence (exonic+intronic) of all genes from chr19.

This part is done with tis R script. Dependencies are the `BSgenome.Hsapiens.UCSC.hg38` and the `ShortRead` packages:

```R
#############################
##
## Generate a fasta file with sequences of all genes from chr19.
## Get coordinates from Biomart
##
#############################
library(BSgenome.Hsapiens.UCSC.hg38)
library(ShortRead)
library(biomaRt)

# retrieve sequences
genes <- getBM(attributes=c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand", "transcript_appris"),
                  filter="chromosome_name", values="19",
                  mart=useMart("ensembl", dataset="hsapiens_gene_ensembl"))
genes <- genes[genes$transcript_appris == "principal1", ]  # get only principal isoforms
genes <- with(genes, GRanges(chromosome_name, IRanges(start_position, end_position), strand, ensembl_gene_id=ensembl_gene_id))
seqlevelsStyle(BSgenome.Hsapiens.UCSC.hg38) <- "Ensembl"
seqs <- as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38, genes))
names(seqs) <- genes$ensembl_gene_id

# save fasta
dir.create("ref")
ShortRead::writeFasta(DNAStringSet(seqs), "ref/genes_chr19.fa")

# duplicate (5x) 100 genes and pretend they're there multiple times, showing thus 5x expression
set.seed(666)
seqs.de <- c(seqs, rep(seqs[sample(1:length(seqs), 100)], 5))
ShortRead::writeFasta(DNAStringSet(seqs.de), "ref/genes_chr19.de.fa")
```

The second part, simulating reads from this reference, is done using samtools' `wgsim`:
Run this bash code to generate 20k 50bp of SE reads for a reference (eg. chr19 Human), 96 replicates

```sh
#!/bin/bash
set -euo pipefail

export READ_LEN=50
export FRAG_LEN=200
export NUM_READS=20000
export ERR_RATE=0.001
export MUT_RATE=0.0001
export INDEL_RATE=0.15
export INDEL_EXTEND_RATE=0.3

export REF="ref/chr19.fa"
export REF_RNA="ref/genes_chr19.fa"
export REF_RNA2="ref/genes_chr19.de.fa"
export RAWDATA="./rawdata"

CORES=16

# get references and annotations
wget -qO- ftp://ftp.ensembl.org/pub/release-98/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.19.fa.gz | gzip -cd > $REF
samtools faidx $REF
ml star/2.6.1b_debian9 && STAR --runMode genomeGenerate --runThreadN $CORES --genomeDir ref --genomeFastaFiles $REF
wget -qO- ftp://ftp.ensembl.org/pub/release-98/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz | gzip -cd | grep "^19" > ref/genes_chr19.gtf
wget -qO- http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz | gzip -cd | awk '{print $3, "\t", $5, "\t", $6, "\t", $2, "\t0\t", $4, "\t", $5, "\t", $6, "\t0\t", $9, "\t", $10, "\t", $11;}' | sed 's/^chr//' | grep "^19" > ref/genes_chr19.bed

function f {
  REPL=${RAWDATA}/sample_$1
  $REF=$2
  SEED=$1

  wgsim \
    -1${READ_LEN} \
    -d${FRAG_LEN} \
    -N${NUM_READS} \
    -e${ERR_RATE} \
    -r${MUT_RATE} \
    -R${INDEL_RATE} \
    -X${INDEL_EXTEND_RATE} \
    -S${SEED} \
    ${REF} ${REPL}.R1.fastq ${REPL}.R2.fastq > ${REPL}_sim.txt

  gzip -c ${REPL}.R1.fastq > ${REPL}.fastq.gz
  rm ${REPL}.R1.fastq ${REPL}.R2.fastq ${REPL}_sim.txt
}
export -f f

parallel --xapply -j $CORES "f {1} {2}" ::: `seq 1 96` ::: $REF_RNA $REF_RNA2
```

