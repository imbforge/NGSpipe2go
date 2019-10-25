Get mouse microRNA experiment from [GSE57138](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57138).
Groenendyk J, Fan X, Peng Z, Ilnytskyy Y et al. Genome-wide analysis of thapsigargin-induced microRNAs and their targets in NIH3T3 cells. Genom Data 2014 Dec;2:325-7. PMID: 26484121

```bash
ml sratoolkit bowtie

# get rawdata
mkdir rawdata
parallel --xapply -j4 "fastq-dump --stdout {1} | gzip > rawdata/{2}.fastq.gz" ::: SRR1269676 SRR1269677 SRR1269678 SRR1269679 ::: control1 control2 thapsigargin1 thapsigargin2

# get mouse reference & annotation
mkdir ref
wget -qO- ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/ncrna/Mus_musculus.GRCm38.ncrna.fa.gz | gzip -cd | sed 's/ /_/g' > ref/mmu.fa
bowtie-build --threads 8 ref/mmu.fa ref/mmu
wget -qO- ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz | gzip -cd > ref/mmu.gtf

# rrna reference
Rscript -e 'biomaRt::exportFASTA(biomaRt::getBM(filters="biotype", values="rRNA", attributes=c("gene_exon_intron", "ensembl_gene_id"), mart=biomaRt::useEnsembl("ensembl", dataset="mmusculus_gene_ensembl")), "ref/rrna.fa")'
bowtie-build --threads 8 ref/rrna.fa ref/rrna

# mirna coordinate annotation (mirbase)
wget -qO- ftp://mirbase.org/pub/mirbase/CURRENT/genomes/mmu.gff3 > ref/mmu_mirna.gff3
Rscript -e 'rtracklayer::export(rtracklayer::import("ref/mmu_mirna.gff3", format="gff3"), "ref/mmu_mirna.gtf", format="gtf")'

# mirna hairpin & mature forms (from mirbase, needed for mirdeep2)
wget -qO- ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz | gzip -cd | perl -lne 'if($_ =~ "^>") { $_ =~ tr/ /_/ } else { $_ =~ tr/actguACTGU/N/c } print $_;' > ref/hairpin.fa
wget -qO- ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz | gzip -cd | perl -lne 'if($_ =~ "^>") { $_ =~ tr/ /_/ } else { $_ =~ tr/actguACTGU/N/c } print $_;' > ref/mature.fa
```
