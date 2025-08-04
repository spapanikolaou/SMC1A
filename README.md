# SMC1A-in-Inflammation

In this git we provide the source code for all analyses performed in our paper:

**"SMC1A, a sex-biased chromatin modifier, acquires specific regulatory function in lupus shaping inflammatory pathways that promote autoimmunity"**

_Despoina Kosmara*, Sofia Papanikolaou*, Chrysoula Stathopoulou, Giannis Vatsellas, Aggelos Banos, Prodromos Sidiropoulos, Matthieu D Lavigne, Panayotis Verginis, Dimitrios Boumpas, Charalampos Spilianakis, Dimitris Konstantopoulos, Christoforos Nikolaou, George Bertsias_

**Abstract**  

A strong female bias is characteristic of systemic lupus erythematosus (SLE), the prototypic systemic autoimmune disease. Through an unbiased analysis of sexually differentially expressed genes, we demonstrate that among various blood immune cells, monocytes show exaggerated female-biased expression of the cohesin complex subunit SMC1A in SLE compared to healthy individuals and those with ankylosing spondylitis, a non-sex-biased autoimmune disorder. SMC1A downregulation impacted the expression of several immune genes, especially those responsive to lupus-relevant stimuli. Integration of the genome-wide landscape of SMC1A binding, chromatin activity and accessibility in lupus-like monocytes, revealed extensive SMC1A redistribution to active gene enhancers of immune/inflammatory pathways, culminating to their transcriptional induction. Analysis of monocyte transcriptomes from female and male SLE patients demonstrated significant enrichment of SMC1A targets among female-biased genes involved in immune/inflammatory response, corroborated by the increased abundance of secreted cytokines by female counterparts. Notably, SMC1A-regulated genes such as IL6 and IL1A exhibited higher expression in female lupus-like monocytes compared to males, with increased SMC1A binding at associated enhancers. Collectively, our study highlights SMC1A as a female-biased chromatin modifier that acquires specific regulatory function during lupus, accentuating the inflammatory pathways and offering mechanistic insights into why females are predisposed to SLE and other autoimmune diseases.

Source Code (in R) is provided [here](https://github.com/christoforos-nikolaou/SMC1A-in-Inflammation/blob/main/final_fig_script.R).

Data Files may be found [here](https://github.com/christoforos-nikolaou/SMC1A-in-Inflammation/tree/main/datafiles).

**System requirements**
The code was developed and tested on the following system:
platform: x86_64-pc-linux-gnu         
arch: x86_64                      
os: linux-gnu                   
system: x86_64, linux-gnu 
R version: 4.3.0 (2023-04-21)

> packageVersion("dplyr")
[1] ‘1.1.4’
> packageVersion("ggplot2")
[1] ‘3.5.2’
> packageVersion("ggrepel")
[1] ‘0.9.6’
> packageVersion("ggvenn")
[1] ‘0.1.10’
> packageVersion("clusterProfiler")
[1] ‘4.8.3’
> packageVersion("org.Hs.eg.db")
[1] ‘3.17.0’
> packageVersion("ggpubr")
[1] ‘0.6.1’
> packageVersion("RColorBrewer")
[1] ‘1.1.3’
> packageVersion("ChIPseeker")
[1] ‘1.36.0’
> packageVersion("TxDb.Hsapiens.UCSC.hg19.knownGene")
[1] ‘3.2.2’
> packageVersion("Homo.sapiens")
[1] ‘1.3.1’
> packageVersion("eulerr")
[1] ‘7.0.2’
> packageVersion("reshape2")
[1] ‘1.4.4’
> packageVersion("pheatmap")
[1] ‘1.0.13’
> packageVersion("readxl")
[1] ‘1.4.5’
> packageVersion("tidyverse")
[1] ‘2.0.0’

**Installation Guide**
Download datafiles and R script
Install R and required packages
Run the R script

