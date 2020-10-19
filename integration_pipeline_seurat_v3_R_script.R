# R Script using Seurat (v3.2.2) for Hybrid Human CoV2 Pipeline
# Created by Hezkiel Nanda and advised by Prof. Steve Shen
# Last update: 19-Oct-2020

# 0. Prerequisite
# ---------------
# Setting working dir
setwd("path/to/project/folder/")
# Library usage
library(Seurat)
library(ggplot2)


# 1. Data loading and seurat object update to v3
# ----------------------------------------------
load("robj/Proj033-matureRNA-lane1-pca.rdata",verbose=T)

rm(d4,d7,d8)
d1 <- UpdateSeuratObject(d1)
# 33276 features across 1604 samples within 1 assay 
# Active assay: RNA (33276 features, 2559 variable features)
d2 <- UpdateSeuratObject(d2)
# 33276 features across 2136 samples within 1 assay 
# Active assay: RNA (33276 features, 2567 variable features)
d3 <- UpdateSeuratObject(d3)
# 33276 features across 1698 samples within 1 assay 
# Active assay: RNA (33276 features, 1670 variable features)
d5 <- UpdateSeuratObject(d5)
# 33276 features across 1724 samples within 1 assay 
# Active assay: RNA (33276 features, 2388 variable features)
d6 <- UpdateSeuratObject(d6)
# 33276 features across 1485 samples within 1 assay 
# Active assay: RNA (33276 features, 1697 variable features)


# 2. Merging and object exploration
# -------------------------------------
d <- merge(d1, c(d2,d3,d5,d6), add.cell.ids=c("d1","d2","d3","d5","d6"))
# 33276 features across 8647 samples within 1 assay 

table(d$orig.ident)
# S1   S2   S3   S5   S6   
# 1604 2136 1698 1724 1485
table(Idents(d))
# 0    1    2    3    4    5 
# 3704 2186 1421  796  511   29
table(Idents(d), d$orig.ident)
#   S1  S2  S3  S5  S6
# 0 628 916 711 634 815
# 1 297 692 401 468 328
# 2 270 424 340 241 146
# 3 246  66 122 231 131
# 4 163  38  95 150  65
# 5   0   0  29   0   0


# 3. Integration (with Anchor)
# -------------------------------
d <- SplitObject(d, split.by = "orig.ident")
d <- d[c("S1","S2","S3","S5","S6")]
# replace variable d with anchor set
d <- FindIntegrationAnchors(object.list = d, dims = 1:16)
# replace variable d with integration data
d <- IntegrateData(anchorset = d, dims = 1:16)

DefaultAssay(d) <- "RNA"

dim(GetAssayData(object = d))
# [1] 33276  8647


# 4. Data normalization and scalling
# ----------------------------------
d <- NormalizeData(object = d, normalization.method = "LogNormalize", scale.factor = 10000)
d <- FindVariableFeatures(object = d, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.1,8), dispersion.cutoff = c(0.1,Inf))
d <- ScaleData(object = d, vars.to.regress = c("nCount_RNA", "percent.mt"))


# 5. Detection of host genes across the integration object
# --------------------------------------------------------
d.mito.gene.names <- grep("^MT-", x = rownames(x = d), value=TRUE)
d.cov2.gene.names <- grep("^cov2-", x = rownames(x = d), value=TRUE) 

d.int.vhostgenes <- setdiff(VariableFeatures(object = d), d.cov2.gene.names)
d.int.vhostgenes <- setdiff(d.int.vhostgenes, d.mito.gene.names)

length(d.int.vhostgenes)
# [1] 1989


# 6. Dimensionality Reduction: PCA with host genes
# --------------------------------------------------------------
d1 <- ScaleData(object = d1, vars.to.regress = c("nUMI", "percent.mito"))
d2 <- ScaleData(object = d2, vars.to.regress = c("nUMI", "percent.mito"))
d3 <- ScaleData(object = d3, vars.to.regress = c("nUMI", "percent.mito"))
d5 <- ScaleData(object = d5, vars.to.regress = c("nUMI", "percent.mito"))
d6 <- ScaleData(object = d6, vars.to.regress = c("nUMI", "percent.mito"))


# 7. Dimensionality Reduction: PCA
# --------------------------------
d <- RunPCA(d, features = d.int.vhostgenes, verbose = FALSE)
Idents(d) <- d$RNA_snn_res.0.15
levels(d)
# [1] "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"


# 8. Determine statistically significant principal components
# -----------------------------------------------------------
d <- JackStraw(object = d, num.replicate = 100)
d <- ScoreJackStraw(d, dims=1:16)


# 9. Cluster the cells (resolution 0.15)
# --------------------------------------
d <- FindNeighbors(d)
d <- FindClusters(object = d, reduction = "pca", dims = 1:16, resolution = 0.15)


# 10. Non-linear Dimensionality Reduction: UMAP
# ---------------------------------------------
d <- RunUMAP(object = d)
Idents(d) <- d$orig.ident

pdf(file="output/integration2/integration_umap_res_0_15_per_sample.pdf")
print(DimPlot(d, reduction="umap",label=F, pt.size = 0.5) + labs(caption="Integration res 0.15"))
dev.off()


# 11. Check for IFN genes
# ------------------------------------
ifn.gene.names <- grep("^IFNL[1-4]|^IFNB1", x = rownames(x = d), value=TRUE) 
# [1] "IFNB1" "IFNL3" "IFNL2" "IFNL1" same for all 5 samples
d.ifn.norm <- log2(Matrix::colSums(GetAssayData(object = d, slot = "counts")[ifn.gene.names, ])/Matrix::colSums(GetAssayData(object = d, slot = "counts")) * 10000 + 1)
d[["norm.ifn.reads_log2"]] <- d.ifn.norm

mean(d$norm.ifn.reads_log2) #[1] 0.02036072

table(d$orig.ident)
# S1   S2   S3   S5   S6 
# 1604 2136 1698 1724 1485 

difn <- subset(d, subset=norm.ifn.reads_log2 > 0.00000001)
table(difn$orig.ident)
# S2 S3 S5 S6 
# 11  5 85  2
table(difn$RNA_snn_res.0.15,difn$orig.ident)
#   S2 S3 S5 S6
# 0  0  0  0  1
# 1  2  1 15  0
# 2  6  2 10  0
# 3  0  0  3  0
# 4  0  0  0  0
# 5  0  0 23  0
# 6  0  0  0  1
# 7  3  2 32  0
# 8  0  0  2  0
# 9  0  0  0  0


# 12. Plot UMAP with IFN cells highlighted (red)
# ----------------------------------------------
dcifn <- WhichCells(d, expression = norm.ifn.reads_log2 > 0.00000001)

pdf(file="output/integration2/integration_umap_IFN_highlight_percentages_res_0_15_red.pdf")
print(DimPlot(d, cells.highlight = dcifn, reduction="umap", cols.highlight = "red",label=T, pt.size = 0.5) + NoLegend() + labs(caption="Integration res 0.15 cells >0.00000001% IFN"))
dev.off()

pdf(file="output/integration2/integration_umap_split_IFN_highlight_percentages_res_0_15_red.pdf", width=15, height=8 )
print(DimPlot(d, cells.highlight = dcifn, cols.highlight = "red", reduction="umap", label=F, pt.size = 0.5, split.by = "RNA_snn_res.0.15") + NoLegend() + labs(caption="Integration res 0.15 cells >0.00000001% IFN"))
print(DimPlot(d, cells.highlight = dcifn, cols.highlight = "red", reduction="umap", label=T, pt.size = 0.5, split.by = "orig.ident") + NoLegend() + labs(caption="Integration res 0.15 cells >0.00000001% IFN"))
dev.off()


# 13. Plot UMAP with CoV2 cells (>0.001%) highlighted (cadetblue)
# ---------------------------------------------------------------
ds0.001 <- subset(d, subset = percent.cov2 > 0.001)

dc0.001 <- WhichCells(d, expression = percent.cov2 > 0.001)

pdf(file="output/integration2/integration_umap_cov2_highlight_percentages_res_0_15_cadetblue.pdf")
print(DimPlot(d, cells.highlight = dc0.001, reduction="umap", cols.highlight = "cadetblue",label=T, pt.size = 0.5) + NoLegend() + labs(caption="Integration res 0.15 cells >0.001% SARS-CoV2"))
dev.off()

pdf(file="output/integration2/integration_umap_split_cov2_highlight_percentages_res_0_15_cadetblue.pdf", width=15, height=8 )
print(DimPlot(d, cells.highlight = dc0.001, cols.highlight = "cadetblue", reduction="umap", label=F, pt.size = 0.5, split.by = "RNA_snn_res.0.15") + NoLegend() + labs(caption="Integration res 0.15 cells >0.001% SARS-CoV2"))
print(DimPlot(d, cells.highlight = dc0.001, cols.highlight = "cadetblue", reduction="umap", label=T, pt.size = 0.5, split.by = "orig.ident") + NoLegend() + labs(caption="Integration res 0.15 cells >0.001% SARS-CoV2"))
dev.off()