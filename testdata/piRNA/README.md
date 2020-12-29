## Download test data

```bash
PROJECT='/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline'
mkdir -p ${PROJECT}; cd ${PROJECT}
ml sratoolkit bowtie

# get rawdata
mkdir -p rawdata/drerio
parallel --xapply -j4 "fastq-dump -X 1000000 --stdout {1} | gzip > rawdata/celegans/{2}.fastq.gz" ::: SRR3231361 SRR3231362 SRR3231358 SRR3231359 ::: input_1 input_2 Tdrd6a_IP_1 Tdrd6a_IP_2

# get worm reference & annotation
mkdir ref
wget -qO- ftp://ftp.ensembl.org/pub/release-98/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.toplevel.fa.gz | gzip -cd > ref/GRCz11.fa
bowtie-build --threads 8 ref/GRCz11.fa ref/GRCz11

wget -qO- ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gzip -cd > ref/GRCz11.gtf

## transposons
wget -qO- ftp://ftp.wormbase.org/pub/wormbase/releases/WS235/species/c_elegans/c_elegans.WS235.annotations.gff3.gz | gzip -cd | grep "WBTrans" | sed 's/CHROMOSOME_//' | sed 's/transposable_element/exon/'  | sed 's/transposon=/gene_id /' | sed 's/$/; gene_biotype "transposon";/' > ref/cel.transposons.gtf

cat ref/cel.gtf ref/cel.transposons.gtf | sort -k1,1 -k4,4n > ref/cel.genes_and_transposons.gtf

# rrna reference
ml R
Rscript -e 'biomaRt::exportFASTA(biomaRt::getBM(filters="biotype", values="rRNA", attributes=c("gene_exon_intron", "ensembl_gene_id"), mart=biomaRt::useEnsembl("ensembl", dataset="celegans_gene_ensembl")), "ref/rrna.fa")'
bowtie-build --threads 8 ref/rrna.fa ref/rrna

```


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

```


```bash
PROJECT='/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline'
cd ${PROJECT}

rsync --exclude='/.git' -r -t -x -v --progress -u -l -z -s /fsimb/groups/imb-kettinggr/adomingues/projects/NGSpipe2go/ ${PROJECT}/NGSpipe2go/
ln -s NGSpipe2go/pipelines/SmallRNAseq/smallrnaseq_ce.pipeline.groovy ./
ln -s NGSpipe2go/pipelines/SmallRNAseq/bpipe.config ./
ln -s NGSpipe2go/modules/SmallRNAseq/essential.vars.groovy ./


ml bpipe/0.9.9.8.slurm
bpipe run smallrnaseq_ce.pipeline.groovy rawdata/celegans/*.fastq.gz
```


