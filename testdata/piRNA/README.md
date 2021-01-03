## Download test data

```bash
PROJECT='/fsimb/groups/imb-kettinggr/adomingues/projects/test_pipeline'
mkdir -p ${PROJECT}; cd ${PROJECT}
ml sratoolkit bowtie

## get rawdata
mkdir -p rawdata/drerio
parallel --xapply -j4 "fastq-dump -X 1000000 --stdout {1} | gzip > rawdata/drerio/{2}.fastq.gz" ::: SRR6345658 SRR6345657 ::: tdrd6a_het tdrd6a_mut

# get zebrafish reference & annotation
mkdir -p ref/drerio
wget -qO- ftp://hgdownload.soe.ucsc.edu/goldenPath/danRer10/bigZips/danRer10.fa.gz | gzip -cd > ref/drerio/danRer10.fa

sbatch --job-name=bt_index  --partition=long --time=10:00:00 --nodes=1 --cpus-per-task=8 --mem-per-cpu=8G --wrap="bowtie-build --threads 8 ref/drerio/danRer10.fa ref/drerio/danRer10"

## rrna reference
ml R
Rscript -e 'biomaRt::exportFASTA(biomaRt::getBM(filters="biotype", values="rRNA", attributes=c("gene_exon_intron", "ensembl_gene_id"), mart=biomaRt::useEnsembl("ensembl", dataset="drerio_gene_ensembl")), "ref/drerio/rrna.fa")'
bowtie-build --threads 8 ref/drerio/rrna.fa ref/drerio/rrna


## Gene annotations
wget -qO- ftp://ftp.ensembl.org/pub/release-87/gtf/danio_rerio/Danio_rerio.GRCz10.87.chr.gtf.gz | gzip -cd > ref/drerio/danRer10.gtf

# convert gene names to UCSC:
wget --no-check-certificate -qO- http://raw.githubusercontent.com/dpryan79/ChromosomeMappings/master/GRCz10_ensembl2UCSC.txt \
   | awk '{if($1!=$2) print "s/^"$1"/"$2"/g"}' > remap.sed
cat remap.sed
cat ref/drerio/danRer10.gtf | sed -f remap.sed  > ref/drerio/danRer10.chr.gtf

cut -f1 ref/drerio/danRer10.chr.gtf | sort | uniq -c


## repeats and transposons
wget -qO- http://www.repeatmasker.org/genomes/danRer10/RepeatMasker-rm406-dfam2.0/danRer10.fa.out.gz | gzip -cd | egrep -v 'Simple_repeat|Low_complexity' > ref/drerio/danRer10.RepeatMasker.fa.out

module load bedtools/2.27.1_debian9 
module load samtools/1.5_debian9 
module load bowtie/0.12.8 
module load RepEnrich/1.2_py3

sbatch --job-name=bt_index  --partition=long --time=10:00:00 --nodes=1 --cpus-per-task=8 --mem-per-cpu=8G --wrap="RepEnrich_setup.py ref/drerio/danRer10.RepeatMasker.fa.out ref/drerio/danRer10.fa ref/drerio/RepEnrich"

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


