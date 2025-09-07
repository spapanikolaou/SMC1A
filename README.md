# SMC1A-in-Inflammation

In this git we provide the source code for all analyses performed in our paper:

**"The X-linked chromatin modifier SMC1A is deregulated in lupus female monocytes and promotes autoimmunity"**

_Despoina Kosmara*, Sofia Papanikolaou*, Chrysoula Stathopoulou, Giannis Vatsellas, Aggelos Banos, Prodromos Sidiropoulos, Matthieu D Lavigne, Panayotis Verginis, Dimitrios Boumpas, Charalampos Spilianakis, Dimitris Konstantopoulos, Christoforos Nikolaou, George Bertsias_

**Abstract**  

Biological sex underlines the onset of autoimmunity yet the underlying molecular mechanisms remain ill-defined. Female susceptibility is prominent in systemic lupus erythematosus (SLE), an autoimmune disease of substantial burden. Herein we demonstrate that SMC1A, a subunit of the cohesin complex that escapes X-chromosome inactivation, exhibits female-biased expression in monocytes from SLE patients and those cultured under lupus-inducing stimulus. By altering SMC1A dosage, the expression of immune-related genes, particularly those responsive to lupus environment, was affected. Active enhancers of immune and inflammatory pathways showed widespread SMC1A redistribution in lupus-like monocytes, as revealed by integrated analysis of SMC1A binding, chromatin activity and accessibility. Our data support SMC1A as a sex-biased chromatin modifier that enhances inflammatory responses in lupus.

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

Download datafiles and R script,
Install R and required packages,
Run the R script.

