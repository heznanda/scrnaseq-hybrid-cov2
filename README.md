# Single Cell RNA Sequencing with Hybrid Pipeline
### Published as Single cell resolution of SARS-CoV-2 tropism, antiviral responses, and susceptibility to therapies in primary human airway epithelium, Langlois et al. 2020

##### The 10x Chromium Single Cell sequencing raw data of experimental samples were reference mapped using Cell Ranger Count (v3.0.1) for alignment to a hybrid human+SARS-CoV-2 viral genome reference (Homo_sapiens.GRCh38; MN985325):
* Concatenating nucleotide and assembly's annotation and gtf files
* Creating Cell Ranger reference
* Alignment and Count using Cell Ranger Count for each of the sample

The hybrid genome reference, annotation file analysis script: 

## Analysis Method:

##### The resulting raw count matrix for each experimental data set was imported into an R pipeline using Seurat (v2.3.4)
where the filter criteria for empty droplets are minimum 1000 genes per cell, for genes are minimum three cells and for mitochondria genes percentage is no more than 30%. The default parameters were used for subsequent data normalization and scaling. The Seurat default method was used to select variable genes for PCA and clustering analysis, which the possible viral and mitochondria genes were excluded from final var gene lists for all samples. The Seurat utilities were used for subsequent cluster marker gene selection and results visualization.
