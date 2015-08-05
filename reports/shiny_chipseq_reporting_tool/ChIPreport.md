---
title: "chipseq"
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
Mapping was done using Bowtie. The mapping statistics show the total number of reads delivered to the aligner, the number of mapped reads, the reads that failed to align, and the number of reads discarded because of mapping to multiple positions:

|         | input_reads|           mapped|      failed|   discarded|     duplicates|
|:--------|-----------:|----------------:|-----------:|-----------:|--------------:|
|ChIP_R1  |     1333473| 1333085 (99.97%)| 338 (0.03%)|  50 (0.00%)| 118586 (8.89%)|
|ChIP_R2  |     1491016| 1490554 (99.97%)| 423 (0.03%)|  39 (0.00%)| 103528 (6.94%)|
|Input_R1 |     2635762| 2635005 (99.97%)| 676 (0.03%)|  81 (0.00%)|  36698 (1.39%)|
|Input_R2 |     1664522| 1663537 (99.94%)| 872 (0.05%)| 113 (0.01%)|  20950 (1.26%)|

### PCR bottleneck coefficient
The PBC (PCR bottleneck coefficient) is an approximate measure of library complexity. Provisionally, 0-0.5 is severe bottlenecking, 0.5-0.8 is moderate bottlenecking, 0.8-0.9 is mild bottlenecking, while 0.9-1.0 is no bottlenecking.

Very low values can indicate a technical problem, such as PCR bias, or a biological finding, such as a very rare genomic feature. Nuclease-based assays (DNase, MNase) detecting features with base-pair resolution (transcription factor footprints, positioned nucleosomes) are expected to recover the same read multiple times, resulting in a lower PBC score for these assays. Note that the most complex library, random DNA, would approach 1.0, thus the very highest values can indicate technical problems with libraries. Some common numbers from ENCODE datasets are: 82% for TF ChIP, 89% for His ChIP, 77% for DNase, 98% for FAIRE, and 97% for control ENCODE datasets have no or mild bottlenecking.

|                                |       PBC|
|:-------------------------------|---------:|
|Sample_imb_gcf_2014_07_ChIP_R1  | 0.9078274|
|Sample_imb_gcf_2014_07_ChIP_R2  | 0.9277504|
|Sample_imb_gcf_2014_07_Input_R1 | 0.9823947|
|Sample_imb_gcf_2014_07_Input_R2 | 0.9850511|

### Raw reads qualities, sequence bias and duplication
|                |Duplication                                                                                                                                                              |Read qualities                                                                                                                                                         |Sequence bias                                                                                                                                                                   |
|:---------------|:------------------------------------------------------------------------------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|ChIP_R1_fastqc  |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_ChIP_R1_fastqc/Images/duplication_levels.png)  |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_ChIP_R1_fastqc/Images/per_base_quality.png)  |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_ChIP_R1_fastqc/Images/per_base_sequence_content.png)  |
|ChIP_R2_fastqc  |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_ChIP_R2_fastqc/Images/duplication_levels.png)  |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_ChIP_R2_fastqc/Images/per_base_quality.png)  |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_ChIP_R2_fastqc/Images/per_base_sequence_content.png)  |
|Input_R1_fastqc |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_Input_R1_fastqc/Images/duplication_levels.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_Input_R1_fastqc/Images/per_base_quality.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_Input_R1_fastqc/Images/per_base_sequence_content.png) |
|Input_R2_fastqc |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_Input_R2_fastqc/Images/duplication_levels.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_Input_R2_fastqc/Images/per_base_quality.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/fastqc/Sample_imb_gcf_2014_07_Input_R2_fastqc/Images/per_base_sequence_content.png) |

### IPstrength
To estimate the IP strength, we attempt to decompose the population of IP reads into two distinct components: those pulled down by the antibody, and background.

|                                                                                                                                           |                                                                                                                                           |                                                                                                                                           |                                                                                                                                           |
|:------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------|:------------------------------------------------------------------------------------------------------------------------------------------|
|![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/ipstrength/chip_r1.vs.input_r1_ipstrength.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/ipstrength/chip_r1.vs.input_r2_ipstrength.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/ipstrength/chip_r2.vs.input_r1_ipstrength.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/ipstrength/chip_r2.vs.input_r2_ipstrength.png) |
|chip_r1.vs.input_r1                                                                                                                        |chip_r1.vs.input_r2                                                                                                                        |chip_r2.vs.input_r1                                                                                                                        |chip_r2.vs.input_r2                                                                                                                        |

The plot is useful to detect several forms of biases:
* The distance between the two curves represents the enrichment of the IP versus the input. Close curves means weak enrichment.
* Undersequencing is shown by the curves close to the x-axis for a large percentage of bins, representing no reads were aligned on this bins. 
* Sequencing bias towards a small set of specific regions is shown by the curves close to the y-axis for a large percentage of tags.

### Cross-correlation analysis

**Cross-correlation**
A measure of enrichment derived without dependence on prior determination of enriched regions. Forward and reverse strand read coverage signal tracks are computed (number of unique mapping read starts at each base in the genome on the + and - strand counted separately). The forward and reverse tracks are shifted towards and away from each other by incremental distances and for each shift, the Pearson correlation coefficient is computed. In this way, a cross-correlation profile is computed, representing the correlation between forward and reverse strand coverage at different shifts. The highest cross-correlation value is obtained at a strand shift equal to the predominant fragment length in the dataset as a result of clustering/enrichment of relative fixed-size fragments around the binding sites of the target factor or feature.

**Normalized Strand Cross-correlation coefficient (NSC):**
The NSC is the ratio of the maximal cross-correlation value (which occurs at strand shift equal to fragment length) divided by the background cross-correlation (minimum cross-correlation value over all possible strand shifts). Higher values indicate more enrichment, values less than 1.1 are relatively low NSC scores, and the minimum possible value is 1 (no enrichment). This score is sensitive to technical effects; for example, high-quality antibodies such as H3K4me3 and CTCF score well for all cell types and ENCODE production groups, and variation in enrichment in particular IPs is detected as stochastic variation. This score is also sensitive to biological effects; narrow marks score higher than broad marks (H3K4me3 vs H3K36me3, H3K27me3) for all cell types and ENCODE production groups, and features present in some individual cells, but not others, in a population are expected to have lower scores.

**Relative Strand Cross-correlation coefficient (RSC):**
The RSC is the ratio of the fragment-length cross-correlation value minus the background cross-correlation value, divided by the phantom-peak cross-correlation value minus the background cross-correlation value. The minimum possible value is 0 (no signal), highly enriched experiments have values greater than 1, and values much less than 1 may indicate low quality.

|                                                                                                                                                              |                                                                                                                                                              |                                                                                                                                                               |                                                                                                                                                               |
|:-------------------------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------------------|
|![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/phantompeak/Sample_imb_gcf_2014_07_ChIP_R1_duprm_phantompeak.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/phantompeak/Sample_imb_gcf_2014_07_ChIP_R2_duprm_phantompeak.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/phantompeak/Sample_imb_gcf_2014_07_Input_R1_duprm_phantompeak.png) |![alt text](/fsimb/groups/imb-bioinfocf/projects/cfb_internal/imbforge-code/test/chipseq/qc/phantompeak/Sample_imb_gcf_2014_07_Input_R2_duprm_phantompeak.png) |
|ChIP_R1_duprm                                                                                                                                                 |ChIP_R2_duprm                                                                                                                                                 |Input_R1_duprm                                                                                                                                                 |Input_R2_duprm                                                                                                                                                 |

## Peaks called

A summary of the peaks called in the different comparisons is shown here:


|     | chip_r1 vs. input_r1| chip_r1 vs. input_r2| chip_r2 vs. input_r1| chip_r2 vs. input_r2|
|:----|--------------------:|--------------------:|--------------------:|--------------------:|
|chr1 |                  139|                  134|                  107|                  106|
