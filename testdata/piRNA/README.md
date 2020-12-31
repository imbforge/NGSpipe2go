## Download test data

```bash
PROJECT='/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline'
mkdir -p ${PROJECT}; cd ${PROJECT}
ml sratoolkit bowtie

## get rawdata
mkdir -p rawdata/drerio
parallel --xapply -j4 "fastq-dump -X 1000000 --stdout {1} | gzip > rawdata/drerio/{2}.fastq.gz" ::: SRR3231361 SRR3231362 SRR3231358 SRR3231359 ::: input_1 input_2 Tdrd6a_IP_1 Tdrd6a_IP_2

# get zebrafish reference & annotation
mkdir -p ref/drerio
wget -qO- ftp://ftp.ensembl.org/pub/release-98/fasta/danio_rerio/dna/Danio_rerio.GRCz11.dna.primary_assembly.fa.gz | gzip -cd > ref/drerio/GRCz11.fa

bowtie-build --threads 8 ref/drerio/GRCz11.fa ref/drerio/GRCz11

## rrna reference
ml R
Rscript -e 'biomaRt::exportFASTA(biomaRt::getBM(filters="biotype", values="rRNA", attributes=c("gene_exon_intron", "ensembl_gene_id"), mart=biomaRt::useEnsembl("ensembl", dataset="drerio_gene_ensembl")), "ref/drerio/rrna.fa")'
bowtie-build --threads 8 ref/drerio/rrna.fa ref/drerio/rrna


## Gene annotations
wget -qO- ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gzip -cd > ref/drerio/GRCz11.gtf

## repeats and transposons
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/danRer11/database/rmsk.txt.gz ref/drerio/

zcat ref/drerio/rmsk.txt.gz > ref/drerio/GRCz11.repeat_masker.txt
cat ref/drerio/GRCz11.repeat_masker.txt | \
    awk -v OFS='\t' '{print $6, $7, $8, $13, $11, $10, $12}' | \
    sort -k1,1 -k2,2n > ref/drerio/GRCz11.repeat_masker.bed 

## sanity check
egrep "LINE|SINE|DNA" ref/drerio/GRCz11.repeat_masker.bed | cut -f 7 | sort | uniq

## convert to Ensembl contig names:
##conda activate base 
##conda install -c bioconda cvbio 
curl https://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCz11_UCSC2ensembl.txt > ref/drerio/GRCz11_UCSC2ensembl.txt

cvbio UpdateContigNames \
    -i ref/drerio/GRCz11.repeat_masker.bed \
    -o ref/drerio/GRCz11.repeat_masker.ens.bed \
    -m ref/drerio/GRCz11_UCSC2ensembl.txt \
    --comment-chars '#' \
    --columns 0 \
    --skip-missing true

## it's essential that column 7 has the class name
egrep "LINE|SINE|DNA" ref/drerio/GRCz11.repeat_masker.ens.bed > ref/drerio/GRCz11.transposons.bed

## create repRenrich reference
## here only DNA, SINE, LINE elements were used
awk -v FS='\t' -v OFS='\t' '{print $1,$2,$3,$5,$7,$4}' ref/drerio/GRCz11.transposons.bed > ref/drerio/GRCz11.transposons_repEnrich.bed

module load bedtools/2.27.1_debian9 
module load samtools/1.5_debian9 
module load bowtie/0.12.8 
module load RepEnrich/1.2_py3

RepEnrich_setup.py ref/drerio/GRCz11.transposons_repEnrich.bed ref/drerio/GRCz11.fa ref/drerio/RepEnrich --is_bed TRUE
```


```R
library("rtracklayer")
gtf_file <- "ref/drerio/cel.gtf"
gtf <- import(gtf_file)
unique(elementMetadata(gtf)$gene_biotype)

## Transposons
transposons <- import("ref/drerio/cel.transposons.gtf")
unique(elementMetadata(transposons)$gene_biotype )

all <- c(gtf, transposons)
unique(elementMetadata(all)$gene_biotype )

## piRNAs/21Us
anno_pirna <- gtf[elementMetadata(gtf)$gene_biotype %in% c("piRNA")]
export(anno_pirna, "ref/drerio/piRNA.pipeline.gtf")

## 22Gs
anno_22G <- gtf[elementMetadata(gtf)$gene_biotype %in% c("protein_coding", "pseudogene", "lincRNA", "transposon")]
export(anno_22G, "ref/drerio/proteincoding_pseudogenes_lincRNA_transposons.pipeline.gtf")

## 26Gs
anno_26G <- gtf[elementMetadata(gtf)$gene_biotype %in% c("protein_coding", "pseudogene", "lincRNA")]
export(anno_26G, "ref/drerio/proteincoding_pseudogenes_lincRNA.pipeline.gtf")

## miRNAs
anno_miRNA <- gtf[elementMetadata(gtf)$gene_biotype %in% c("protein_coding", "pseudogene", "lincRNA")]
export(anno_miRNA, "ref/drerio/miRNA.pipeline.gtf")

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


