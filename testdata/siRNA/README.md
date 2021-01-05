## Download test data

```bash
PROJECT='/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline'
mkdir -p ${PROJECT}; cd ${PROJECT}
ml sratoolkit bowtie

# get rawdata
mkdir -p rawdata/celegans
parallel --xapply -j4 "fastq-dump -X 10000000 --stdout {1} | gzip > rawdata/celegans/{2}.fastq.gz" ::: SRR6002308 SRR6002311 SRR6002314 SRR6854649 SRR6854650 SRR6854651 ::: N2_rep1 N2_rep2 N2_rep3 RppH_rep1 RppH_rep2 RppH_rep3

## also with the sensor
mkdir -p rawdata/celegans_sensor
parallel --xapply -j4 "fastq-dump -X 10000000 --stdout {1} | gzip > rawdata/celegans_sensor/{2}.fastq.gz" ::: SRR11396359 SRR11396358 SRR11396351 SRR11396357 ::: pid1_rep1 pid1_rep2 prg1_pid2_rep1 prg1_pid2_rep2


# get worm reference & annotation
mkdir ref
wget -qO- ftp://ftp.ensembl.org/pub/release-98/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz | gzip -cd > ref/cel.fa
bowtie-build --threads 8 ref/cel.fa ref/cel

wget -qO- http://ftp.ensembl.org/pub/release-98/gtf/caenorhabditis_elegans/Caenorhabditis_elegans.WBcel235.98.gtf.gz | gzip -cd > ref/cel.gtf

## transposons
wget -qO- ftp://ftp.wormbase.org/pub/wormbase/releases/WS235/species/c_elegans/c_elegans.WS235.annotations.gff3.gz | gzip -cd | grep "WBTrans" | sed 's/CHROMOSOME_//' | sed 's/transposable_element/exon/'  | sed 's/transposon=/gene_id /' | sed 's/$/; gene_biotype "transposon";/' > ref/cel.transposons.gtf

cat ref/cel.gtf ref/cel.transposons.gtf | sort -k1,1 -k4,4n > ref/cel.genes_and_transposons.gtf

# rrna reference
ml R
Rscript -e 'biomaRt::exportFASTA(biomaRt::getBM(filters="biotype", values="rRNA", attributes=c("gene_exon_intron", "ensembl_gene_id"), mart=biomaRt::useEnsembl("ensembl", dataset="celegans_gene_ensembl")), "ref/rrna.fa")'
bowtie-build --threads 8 ref/rrna.fa ref/rrna

```

This will create the GTF annotations needed to filter small RNA classes:

```R
library("rtracklayer")
gtf_file <- "ref/cel.gtf"
gtf <- import(gtf_file)
unique(elementMetadata(gtf)$gene_biotype)

## Transposons
transposons <- import("ref/cel.transposons.gtf")
unique(elementMetadata(transposons)$gene_biotype )

all <- c(gtf, transposons)
unique(elementMetadata(all)$gene_biotype )

## piRNAs/21Us
anno_pirna <- gtf[elementMetadata(gtf)$gene_biotype %in% c("piRNA")]
export(anno_pirna, "ref/piRNA.pipeline.gtf")

## 22Gs
anno_22G <- gtf[elementMetadata(gtf)$gene_biotype %in% c("protein_coding", "pseudogene", "lincRNA", "transposon")]
export(anno_22G, "ref/proteincoding_pseudogenes_lincRNA_transposons.pipeline.gtf")

## 26Gs
anno_26G <- gtf[elementMetadata(gtf)$gene_biotype %in% c("protein_coding", "pseudogene", "lincRNA")]
export(anno_26G, "ref/proteincoding_pseudogenes_lincRNA.pipeline.gtf")

## miRNAs
anno_miRNA <- gtf[elementMetadata(gtf)$gene_biotype %in% c("protein_coding", "pseudogene", "lincRNA")]
export(anno_miRNA, "ref/miRNA.pipeline.gtf")

## structural RNAs
anno_structural <- gtf[elementMetadata(gtf)$gene_biotype %in% c("rRNA", "tRNA", "snoRNA", "snRNA")]
export(anno_structural, "ref/structural.pipeline.gtf")

```


## Run the pipeline

```bash
PROJECT='/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline'
cd ${PROJECT}

rsync --exclude='/.git' -r -t -x -v --progress -u -l -z -s /fsimb/groups/imb-kettinggr/adomingues/projects/NGSpipe2go/ ${PROJECT}/NGSpipe2go/

ln -s NGSpipe2go/pipelines/SmallRNAseq/* .

ml bpipe/0.9.9.8.slurm
bpipe run siRNA.pipeline.groovy rawdata/celegans/*.fastq.gz
bpipe run siRNA_sensor.pipeline.groovy rawdata/celegans_sensor/*.fastq.gz
```
