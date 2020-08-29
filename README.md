# Single Cell RNA Sequencing with Hybrid Pipeline
## Published as Single cell resolution of SARS-CoV-2 tropism, antiviral responses, and susceptibility to therapies in primary human airway epithelium, Langlois et al. 2020
The journal link: <a class="link-gray" href="#">Publishing...</a>

### The 10x Chromium Single Cell sequencing raw data of experimental samples were reference mapped using Cell Ranger Count (v3.0.1)
The alignment is done in a hybrid manner with human+SARS-CoV-2 viral genome reference (Homo_sapiens.GRCh38; MN985325)

The hybrid genome bash script: <a class="link-gray" href="https://github.com/heznanda/scrnaseq-hybrid-cov2/blob/master/hybrid_pipeline_bash_script.sh">hybrid_pipeline_bash_script.sh</a>

The raw data can be found in NCBI GEO website: <a class="link-gray" href="#">Publishing...</a>

This NCBI GEO website includes raw count matrix for each experimental data set.

###### The overall hybrid pipeline bash steps:
* Concatenating nucleotide and assembly's annotation and gtf files
* Creating Cell Ranger reference
* Alignment and Count using Cell Ranger Count for each of the sample (the filter criteria for empty droplets are minimum 1000 genes per cell)

### Analysis Method for the raw count matrix was done using R (v3.6.3) with Seurat (v2.3.4)

The Seurat hybrid pipeline: <a class="link-gray" href="https://github.com/heznanda/scrnaseq-hybrid-cov2/blob/master/hybrid_pipeline_seurat_R_script.R">hybrid_pipeline_seurat_R_script.R</a>

The CoV2 cells analysis (an extension from the Seurat hybrid pipeline: <a class="link-gray" href="https://github.com/heznanda/scrnaseq-hybrid-cov2/blob/master/CoV2_cells_analysis_seurat_R_script.R">CoV2_cells_analysis_seurat_R_script.R</a>

###### The overall hybrid pipeline Seurat steps:
* QC and filtering on metadata (filter for genes with minimum three cells and for mitochondria genes percentage less than 30%)
* Data normalization, scaling, and detection of variable genes
* Dimensionality Reduction using PCA and clustering analysis (mitochondria and viral genes excluded)
* Find all gene markers and key genes results visualization

###### The overall CoV2 cells analysis steps:
* Subset data for lower threshold of 0.0000001, 0.000001, 0.00001, 0.0001, 0.001, 0.01, 0.1% CoV2 cells
* Plot TSNE dimensionality reduction with specific threshold virus cells highlighted
