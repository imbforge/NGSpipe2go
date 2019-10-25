Get mouse microRNA experiment from [GSE57138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57138).
Groenendyk J, Fan X, Peng Z, Ilnytskyy Y et al. Genome-wide analysis of thapsigargin-induced microRNAs and their targets in NIH3T3 cells. Genom Data 2014 Dec;2:325-7. PMID: 26484121

```bash
ml sratoolkit bowtie

# get rawdata
mkdir rawdata
parallel --xapply -j4 "fastq-dump --stdout {1} | gzip > rawdata/{2}.fastq.gz" ::: SRR1269676 SRR1269677 SRR1269678 SRR1269679 ::: control1 control2 thapsigargin1 thapsigargin2

# get mouse reference & annotation
mkdir ref
wget -qO- ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.chromosome.19.fa.gz | gzip -cd > ref/mmu_chr19.fa
bowtie-build --threads 8 ref/mmu_chr19.fa ref/mmu_chr19
wget -qO- ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz | gzip -cd | grep "^19"> ref/mmu_chr19.gtf

# rrna reference
Rscript -e 'biomaRt::exportFASTA(biomaRt::getBM(filters="biotype", values="rRNA", attributes=c("gene_exon_intron", "ensembl_gene_id"), mart=biomaRt::useEnsembl("ensembl", dataset="mmusculus_gene_ensembl")), "ref/rrna.fa")'
bowtie-build --threads 8 ref/rrna.fa ref/rrna
```
