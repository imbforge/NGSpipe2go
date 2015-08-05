---
title: "SHINYREPS_PROJECT"
output:
  html_document:
    toc: yes
---



## Quality Control RNAseq ##

### Sequencing quality ###
The sequencing quality of the run was good, and the read distribution over the libraries was good. See tables the SAV quality tables:

 
Desc | Value
:--: | :---:
Machine | HWI-ST558
Run Folder | 141010-SN558-0265-AC49HNACXX
Chip ID | unknown 

Lane | Lane Yield (kbases) | Clusters (raw) | Clusters (PF) | First Cycle Int (PF) | \% intensity after 20 cycles (PF) | \% PF Clusters
:--: | :-----------------: | :------------: | :-----------: | :------------------: | :--------------------------: | :------------:
1  | 9650991 | 2201886 $\pm$ 330864 | 1971199 $\pm$ 243511 | 1131 $\pm$ 126 | 79.50 $\pm$ 1.41 | 89.90 $\pm$ 3.11 
2  | 8901946 | 2826947 $\pm$ 152352 | 1856893 $\pm$ 491548 | 949 $\pm$ 210 | 75.97 $\pm$ 14.31 | 66.08 $\pm$ 18.49 
3  | 8314272 | 1836538 $\pm$ 270076 | 1698176 $\pm$ 215295 | 1002 $\pm$ 113 | 79.78 $\pm$ 5.05 | 92.74 $\pm$ 2.04 
4  | 10036281 | 2356506 $\pm$ 334192 | 2049894 $\pm$ 260168 | 1282 $\pm$ 189 | 66.48 $\pm$ 6.95 | 87.59 $\pm$ 7.99 
5  | 10837353 | 2652785 $\pm$ 175894 | 2213511 $\pm$ 84079 | 1030 $\pm$ 105 | 78.68 $\pm$ 2.15 | 83.66 $\pm$ 3.96 
6  | 9773962 | 2222935 $\pm$ 439671 | 1996315 $\pm$ 314381 | 1096 $\pm$ 110 | 76.58 $\pm$ 1.84 | 90.49 $\pm$ 3.65 
7  | 10894107 | 2540112 $\pm$ 527322 | 2225103 $\pm$ 341898 | 1069 $\pm$ 107 | 77.01 $\pm$ 1.55 | 88.62 $\pm$ 5.39 
8  | 10914079 | 2536972 $\pm$ 501254 | 2229182 $\pm$ 326316 | 1069 $\pm$ 101 | 77.08 $\pm$ 1.52 | 88.76 $\pm$ 4.86 

_ | **Average** | 2396835 | 2030034 | 1079 | 76.39 | 85.98 

Lane | Clusters(tile $\mu$,raw) | \% Phas | \% Prephas | \% Retained(raw) | Cyc2-4 $\mu$ Int(raw,PF) | Cyc2-10 $\mu$ \% Loss(filt,PF) | Cyc10-20 $\mu$ \% Loss(filt,PF)
:--: | :----------------------: | :-----: | :--------: | :--------------: | :------------------------: | :----------------------------: | :-----------------------------:
1 | 2201886.01 | 0.2412 | 0.0944 | 89.90 | 1056.01 $\pm$ 111.56 | 1.49 $\pm$ 0.19 | 0.81 $\pm$ 0.08 
2 | 2826947.72 | 0.2412 | 0.0944 | 66.08 | 916.64 $\pm$ 107.86 | -4.29 $\pm$ 17.52 | 0.92 $\pm$ 0.19 
3 | 1836538.88 | 0.2412 | 0.0944 | 92.74 | 937.78 $\pm$ 100.02 | 0.63 $\pm$ 0.16 | 0.77 $\pm$ 0.16 
4 | 2356506.43 | 0.2412 | 0.0944 | 87.59 | 996.38 $\pm$ 105.28 | 0.85 $\pm$ 0.11 | 0.75 $\pm$ 0.12 
5 | 2652785.95 | 0.2412 | 0.0944 | 83.66 | 959.75 $\pm$ 104.07 | 0.43 $\pm$ 0.10 | 0.87 $\pm$ 0.14 
6 | 2222935.44 | 0.2412 | 0.0944 | 90.49 | 991.34 $\pm$ 104.88 | -0.01 $\pm$ 0.13 | 1.32 $\pm$ 0.10 
7 | 2540112.88 | 0.2412 | 0.0944 | 88.62 | 976.06 $\pm$ 105.42 | 0.25 $\pm$ 0.18 | 0.81 $\pm$ 0.12 
8 | 2536972.71 | 0.2412 | 0.0944 | 88.76 | 977.67 $\pm$ 98.56 | 0.25 $\pm$ 0.09 | 0.82 $\pm$ 0.10 

### Mapping
Mapping was done using STAR. The program version and genome assembly are described in the following table. The following parameters were used to get only uniquely mapping reads:

|                       |                                   parms|
|:----------------------|---------------------------------------:|
|version                |                      STAR_2.3.1z13_r470|
|runMode                |                              alignReads|
|limitGenomeGenerateRAM |                             31000000000|
|limitIObufferSize      |                               150000000|
|genomeDir              | /fsimb/groups/imb-bi...tar_genomes/mm9/|
|runThreadN             |                                       8|
|outFilterMismatchNmax  |                                       2|
|outFilterMultimapNmax  |                                      10|
|genomeLoad             |                          NoSharedMemory|
|alignIntronMin         |                                      21|
|outStd                 |                                     SAM|
|outSAMattributes       |                                Standard|
|outSJfilterReads       |                                  Unique|
|sjdbGTFfile            | /fsimb/groups/imb-bi.../Genes/genes.gtf|
|readFilesCommand       |                                    zcat|
|versionGenome          |                                   20101|
|genomeFastaFiles       |                                  mm9.fa|
|genomeSAindexNbases    |                                      14|
|genomeChrBinNbits      |                                      18|
|genomeSAsparseD        |                                       1|
|sjdbOverhang           |                                     100|
|sjdbFileChrStartEnd    |                                       -|

The mapping statistics show the total number of reads delivered to the aligner, the number of uniquely mapped reads, and the number of reads discarded because of mapping to multiple positions or not mapping at all to the reference genome.

|                       | input_reads|      uniq_mapped|      multimapped| unmapped|
|:----------------------|-----------:|----------------:|----------------:|--------:|
|12_mESCs_shMed12_rep_1 |     4839552| 2927429 (60.49%)| 1912123 (39.51%)|       0%|
|19_mESCs_shNMC_rep_2   |     5678147| 3582292 (63.09%)| 2095855 (36.91%)|       0%|
|21_mESCs_shMed12_rep_2 |     4321806| 2656616 (61.47%)| 1665190 (38.53%)|       0%|
|28_mESCs_shNMC_rep_1   |     5135863| 3284015 (63.94%)| 1851848 (36.06%)|       0%|
|37_mESCs_shNMC_rep_3   |     4295755| 2744295 (63.88%)| 1551460 (36.12%)|       0%|
|39_mESCs_shMed12_rep_3 |     5464407| 3386132 (61.97%)| 2078275 (38.03%)|       0%|

### Raw reads qualities, sequence bias and duplication
|                              |Duplication                                                                                                                                                                       |Read qualities                                                                                                                                                                  |Sequence bias                                                                                                                                                                            |
|:-----------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|12_mESCs_shMed12_rep_1_fastqc |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_12_mESCs_shMed12_rep_1_fastqc/Images/duplication_levels.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_12_mESCs_shMed12_rep_1_fastqc/Images/per_base_quality.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_12_mESCs_shMed12_rep_1_fastqc/Images/per_base_sequence_content.png) |
|19_mESCs_shNMC_rep_2_fastqc   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_19_mESCs_shNMC_rep_2_fastqc/Images/duplication_levels.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_19_mESCs_shNMC_rep_2_fastqc/Images/per_base_quality.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_19_mESCs_shNMC_rep_2_fastqc/Images/per_base_sequence_content.png)   |
|21_mESCs_shMed12_rep_2_fastqc |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_21_mESCs_shMed12_rep_2_fastqc/Images/duplication_levels.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_21_mESCs_shMed12_rep_2_fastqc/Images/per_base_quality.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_21_mESCs_shMed12_rep_2_fastqc/Images/per_base_sequence_content.png) |
|28_mESCs_shNMC_rep_1_fastqc   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_28_mESCs_shNMC_rep_1_fastqc/Images/duplication_levels.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_28_mESCs_shNMC_rep_1_fastqc/Images/per_base_quality.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_28_mESCs_shNMC_rep_1_fastqc/Images/per_base_sequence_content.png)   |
|37_mESCs_shNMC_rep_3_fastqc   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_37_mESCs_shNMC_rep_3_fastqc/Images/duplication_levels.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_37_mESCs_shNMC_rep_3_fastqc/Images/per_base_quality.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_37_mESCs_shNMC_rep_3_fastqc/Images/per_base_sequence_content.png)   |
|39_mESCs_shMed12_rep_3_fastqc |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_39_mESCs_shMed12_rep_3_fastqc/Images/duplication_levels.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_39_mESCs_shMed12_rep_3_fastqc/Images/per_base_quality.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/fastqc/Sample_imb_richly_2014_05_39_mESCs_shMed12_rep_3_fastqc/Images/per_base_sequence_content.png) |

### PCR duplication assessment
Measuring the fraction of duplicated reads is a common quality control step for NGS data sets, which can hint towards too low library complexity due to too many PCR cycles during library preparation or other problems in the sequencing workflow. However in RNA-Seq duplicate reads arise naturally due to highly expressed genes, which renders the overall read duplication rate useless.

|                                                                                                                                                                  |                                                                                                                                                                  |                                                                                                                                                                  |                                                                                                                                                                |
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------|
|![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/dupRadar/Sample_imb_richly_2014_05_12_mESCs_shMed12_rep_1_duprm_dupRadar.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/dupRadar/Sample_imb_richly_2014_05_19_mESCs_shNMC_rep_2_duprm_dupRadar.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/dupRadar/Sample_imb_richly_2014_05_21_mESCs_shMed12_rep_2_duprm_dupRadar.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/dupRadar/Sample_imb_richly_2014_05_28_mESCs_shNMC_rep_1_duprm_dupRadar.png) |
|12_mESCs_shMed12_rep_1_duprm                                                                                                                                      |19_mESCs_shNMC_rep_2_duprm                                                                                                                                        |21_mESCs_shMed12_rep_2_duprm                                                                                                                                      |28_mESCs_shNMC_rep_1_duprm                                                                                                                                      |
|![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/dupRadar/Sample_imb_richly_2014_05_37_mESCs_shNMC_rep_3_duprm_dupRadar.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/dupRadar/Sample_imb_richly_2014_05_39_mESCs_shMed12_rep_3_duprm_dupRadar.png) |                                                                                                                                                                  |                                                                                                                                                                |
|37_mESCs_shNMC_rep_3_duprm                                                                                                                                        |39_mESCs_shMed12_rep_3_duprm                                                                                                                                      |                                                                                                                                                                  |                                                                                                                                                                |

### RNA types
Check for an enrichment of specific RNA types. This plots helps in detecting if the extraction protocol worked well, and the expected RNA types were sequenced.

|V1                                                                                                                    |V2                                                                                                                    |V3                                                                                                                    |
|:---------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------|
|![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/RNAtypes/RNAtypes.counts.per.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/RNAtypes/RNAtypes.counts.raw.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/RNAtypes/RNAtypes.counts.rpk.png) |

### Gene body coverage
Check if reads coverage is uniform and if there is any 5' or 3' bias:

|                                                                                                                                                                                    |                                                                                                                                                                                    |                                                                                                                                                                                    |                                                                                                                                                                                  |
|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/geneBodyCov/Sample_imb_richly_2014_05_12_mESCs_shMed12_rep_1_duprm.geneBodyCoverage.curves.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/geneBodyCov/Sample_imb_richly_2014_05_19_mESCs_shNMC_rep_2_duprm.geneBodyCoverage.curves.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/geneBodyCov/Sample_imb_richly_2014_05_21_mESCs_shMed12_rep_2_duprm.geneBodyCoverage.curves.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/geneBodyCov/Sample_imb_richly_2014_05_28_mESCs_shNMC_rep_1_duprm.geneBodyCoverage.curves.png) |
|12_mESCs_shMed12_rep_1_duprm                                                                                                                                                        |19_mESCs_shNMC_rep_2_duprm                                                                                                                                                          |21_mESCs_shMed12_rep_2_duprm                                                                                                                                                        |28_mESCs_shNMC_rep_1_duprm                                                                                                                                                        |
|![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/geneBodyCov/Sample_imb_richly_2014_05_37_mESCs_shNMC_rep_3_duprm.geneBodyCoverage.curves.png)   |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/qc/geneBodyCov/Sample_imb_richly_2014_05_39_mESCs_shMed12_rep_3_duprm.geneBodyCoverage.curves.png) |                                                                                                                                                                                    |                                                                                                                                                                                  |
|37_mESCs_shNMC_rep_3_duprm                                                                                                                                                          |39_mESCs_shMed12_rep_3_duprm                                                                                                                                                        |                                                                                                                                                                                    |                                                                                                                                                                                  |

### Batch effect
The MDS plot shows the sample relations based on multidimensional scaling (MDS). Distances in the plot can be interpreted in terms of leading log2 fold change. Overall, we see sample replicates that show little distances within them, while groups show higher distance between them. The ratio between/within groups distances is high which means there is no batch effect (or technical variation) if the experiment was well randomized:

![plot of chunk unnamed-chunk-9](figure/unnamed-chunk-9-1.png) 

Clustering the samples based on the distance calculated on the expression levels of the 50 most variable genes, one can have an idea which samples have similar gene expression profile. This, together with the MDS plot, enables us to detect if the replicates behave similar in terms of expression and spot batches of samples related by other characteristics rather than condition:

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 

Finally, sample to sample correlation based on the distance calculated on the expression profile gives us a measure of sample similarity:

![plot of chunk unnamed-chunk-11](figure/unnamed-chunk-11-1.png) 

## Differential Expression Analysis

###  Med12vsNMC

There are 329 
genes highlighted in the MA plot with corrected pvalue (Benjamini and Hochberg, FDR) < 
0.01
![plot of chunk unnamed-chunk-12](figure/unnamed-chunk-12-1.png) 

The table contains the top 25 differentially expressed genes in terms of fold change and Pvalue:

|              | logFC| logCPM|     LR| PValue| FDR|
|:-------------|-----:|------:|------:|------:|---:|
|Gprc6a        |  5.61|   0.89|  20.48|      0|   0|
|Dcn           |  5.26|   1.98|  36.16|      0|   0|
|Srgn          | -5.20|   6.77| 412.31|      0|   0|
|Tac2          | -5.16|   2.02|  42.05|      0|   0|
|Mettl7b       |  4.79|   2.28|  55.60|      0|   0|
|Nav3          |  3.78|   2.46|  50.57|      0|   0|
|Eya4          |  3.68|   4.62| 167.49|      0|   0|
|Col6a1        |  3.58|   4.95| 204.39|      0|   0|
|Tmcc3         |  3.47|   8.88| 621.95|      0|   0|
|AI317395      |  3.35|   2.33|  36.43|      0|   0|
|Rfx6          |  3.32|   2.64|  42.40|      0|   0|
|Best3         | -3.26|   2.18|  25.64|      0|   0|
|Fam162b       |  3.20|   1.99|  28.23|      0|   0|
|Slc5a4a       | -3.20|   1.59|  19.73|      0|   0|
|Gm15915       | -3.13|   2.08|  30.81|      0|   0|
|Frk           |  3.12|   4.94| 107.04|      0|   0|
|Nts           |  3.11|   1.70|  22.63|      0|   0|
|Syt1          |  3.02|   5.51| 150.38|      0|   0|
|Mettl24       |  3.01|   2.23|  31.42|      0|   0|
|Hkdc1         | -3.00|   4.36| 102.84|      0|   0|
|Tmem200a      |  2.97|   4.94|  95.01|      0|   0|
|C920009B18Rik | -2.91|   1.68|  20.83|      0|   0|
|Hal           | -2.80|   5.91| 202.23|      0|   0|
|Cdh23         | -2.69|   5.54| 123.34|      0|   0|
|Rgs17         |  2.67|   5.05|  54.15|      0|   0|
