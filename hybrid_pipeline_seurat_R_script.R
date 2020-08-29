# R Script using Seurat (v2.3.4) for Hybrid Human CoV2 Pipeline
# Created by Prof. Steve Shen and organized by Hezkiel Nanda
# Last update: 29-Aug-2020

# 0. Prerequisite
# ---------------
# Setting working dir
setwd("path/to/project/folder/")

# Library usage
library(Seurat)

# 1. Data loading and seurat object creation
# ------------------------------------------
S1 <- Read10X(data.dir = "./S1a")
S2 <- Read10X(data.dir = "./S2a")
S3 <- Read10X(data.dir = "./S3a")
S5 <- Read10X(data.dir = "./S5a")
S6 <- Read10X(data.dir = "./S6a")

d1 <- CreateSeuratObject(raw.data = S1, project = "S1")
d2 <- CreateSeuratObject(raw.data = S2, project = "S2")
d3 <- CreateSeuratObject(raw.data = S3, project = "S3")
d5 <- CreateSeuratObject(raw.data = S5, project = "S5")
d6 <- CreateSeuratObject(raw.data = S6, project = "S6")


# 2. Raw data exploration 
# -----------------------
dim(d1@raw.data)
# [1] 33276  2020
dim(d2@raw.data)
# [1] 33276  2902
dim(d3@raw.data)
# [1] 33276  1881
dim(d5@raw.data)
# [1] 33276  2163
dim(d6@raw.data)
# [1] 33276  1780


# 3. QC and filtering on metadata
# -------------------------------
# Get gene names
mito.gene.names <- grep("^MT-", x = rownames(x = d1@data), value=TRUE) #case sensitive MT
cov2.gene.names <- grep("^cov2-", x = rownames(x = d1@data), value=TRUE) 

# Get TSS normalized mitochodrial counts
# Cell total expression
d1.total.counts <- Matrix::colSums(d1@raw.data) 
d2.total.counts <- Matrix::colSums(d2@raw.data)
d3.total.counts <- Matrix::colSums(d3@raw.data)
d5.total.counts <- Matrix::colSums(d5@raw.data) 
d6.total.counts <- Matrix::colSums(d6@raw.data)

sum(d1.total.counts)
# [1] 70595432
sum(d2.total.counts)
# [1] 71829322
sum(d3.total.counts)
# [1] 78541729
sum(d5.total.counts)
# [1] 56014354
sum(d6.total.counts)
# [1] 71512188

# Scale each count by the cell total
d1.mito.percent <- Matrix::colSums(d1@raw.data[mito.gene.names, ])/d1.total.counts
d2.mito.percent <- Matrix::colSums(d2@raw.data[mito.gene.names, ])/d2.total.counts 
d3.mito.percent <- Matrix::colSums(d3@raw.data[mito.gene.names, ])/d3.total.counts 
d5.mito.percent <- Matrix::colSums(d5@raw.data[mito.gene.names, ])/d5.total.counts 
d6.mito.percent <- Matrix::colSums(d6@raw.data[mito.gene.names, ])/d6.total.counts 

# Cov2 percent
d1.cov2.percent <- Matrix::colSums(d1@raw.data[cov2.gene.names, ])/d1.total.counts 
d2.cov2.percent <- Matrix::colSums(d2@raw.data[cov2.gene.names, ])/d2.total.counts 
d3.cov2.percent <- Matrix::colSums(d3@raw.data[cov2.gene.names, ])/d3.total.counts 
d5.cov2.percent <- Matrix::colSums(d5@raw.data[cov2.gene.names, ])/d5.total.counts 
d6.cov2.percent <- Matrix::colSums(d6@raw.data[cov2.gene.names, ])/d6.total.counts 

# Add to seurat object as a metadata
d1 <- AddMetaData(object = d1, metadata = d1.mito.percent, col.name = "percent.mito")
d2 <- AddMetaData(object = d2, metadata = d2.mito.percent, col.name = "percent.mito")
d3 <- AddMetaData(object = d3, metadata = d3.mito.percent, col.name = "percent.mito")
d5 <- AddMetaData(object = d5, metadata = d5.mito.percent, col.name = "percent.mito")
d6 <- AddMetaData(object = d6, metadata = d6.mito.percent, col.name = "percent.mito")

# cov2
d1 <- AddMetaData(object = d1, metadata = d1.cov2.percent, col.name = "percent.cov2")
d2 <- AddMetaData(object = d2, metadata = d2.cov2.percent, col.name = "percent.cov2")
d3 <- AddMetaData(object = d3, metadata = d3.cov2.percent, col.name = "percent.cov2")
d5 <- AddMetaData(object = d5, metadata = d5.cov2.percent, col.name = "percent.cov2")
d6 <- AddMetaData(object = d6, metadata = d6.cov2.percent, col.name = "percent.cov2")

# Note: comparison between mature RNA mapping and all RNA mapping drives this conclusion, premature RNA mapping yields more cell but less genes.
# focusing on preRNA mapping strategy for production mode, 06/16/2020

# Filter max complexity
d1 <- FilterCells(object = d1, subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(12000, 0.3))
d2 <- FilterCells(object = d2, subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(12000, 0.3))
d3 <- FilterCells(object = d3, subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(12000, 0.3))
d5 <- FilterCells(object = d5, subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(12000, 0.3))
d6 <- FilterCells(object = d6, subset.names = c("nGene", "percent.mito"), low.thresholds = c(1000, -Inf), high.thresholds = c(12000, 0.3))

# low 1000
dim(d1@data)
# [1] 33276  1604
dim(d2@data)
# [1] 33276  2136
dim(d3@data)
# [1] 33276  1698
dim(d5@data)
# [1] 33276  1724
dim(d6@data)
# [1] 33276  1485


# 4. Data normalization
# ---------------------
d1 <- NormalizeData(object = d1, normalization.method = "LogNormalize", scale.factor = 10000)
d2 <- NormalizeData(object = d2, normalization.method = "LogNormalize", scale.factor = 10000)
d3 <- NormalizeData(object = d3, normalization.method = "LogNormalize", scale.factor = 10000)
d5 <- NormalizeData(object = d5, normalization.method = "LogNormalize", scale.factor = 10000)
d6 <- NormalizeData(object = d6, normalization.method = "LogNormalize", scale.factor = 10000)


# 5. Detection of variable genes across the single cells
# ------------------------------------------------------
d1 <- FindVariableGenes(object = d1, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.25, y.high.cutoff = Inf)
d2 <- FindVariableGenes(object = d2, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.1, y.high.cutoff = Inf)
d3 <- FindVariableGenes(object = d3, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.1, y.high.cutoff = Inf)
d5 <- FindVariableGenes(object = d5, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.1, y.high.cutoff = Inf)
d6 <- FindVariableGenes(object = d6, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.1, x.high.cutoff = 8, y.cutoff = 0.1, y.high.cutoff = Inf)

length(x = d1@var.genes)
# [1] 2559
length(x = d2@var.genes)
# [1] 2567
length(x = d3@var.genes)
# [1] 1670
length(x = d5@var.genes)
# [1] 2388
length(x = d6@var.genes)
# [1] 1697

# exclude cov2 genes and mito genes (no influence to host)
d1.vhostgenes <- setdiff(d1@var.genes, cov2.gene.names); d1.vhostgenes <- setdiff(d1.vhostgenes, mito.gene.names)
d2.vhostgenes <- setdiff(d2@var.genes, cov2.gene.names); d2.vhostgenes <- setdiff(d2.vhostgenes, mito.gene.names)
d3.vhostgenes <- setdiff(d3@var.genes, cov2.gene.names); d3.vhostgenes <- setdiff(d3.vhostgenes, mito.gene.names)
d5.vhostgenes <- setdiff(d5@var.genes, cov2.gene.names); d5.vhostgenes <- setdiff(d5.vhostgenes, mito.gene.names)
d6.vhostgenes <- setdiff(d6@var.genes, cov2.gene.names); d6.vhostgenes <- setdiff(d6.vhostgenes, mito.gene.names)

length(d1.vhostgenes)
# [1] 2550
length(d2.vhostgenes)
# [1] 2553
length(d3.vhostgenes)
# [1] 1664
length(d5.vhostgenes)
# [1] 2376
length(d6.vhostgenes)
# [1] 1689


# 6. Scaling the data and removing unwanted sources of variation
# --------------------------------------------------------------
d1 <- ScaleData(object = d1, vars.to.regress = c("nUMI", "percent.mito"))
d2 <- ScaleData(object = d2, vars.to.regress = c("nUMI", "percent.mito"))
d3 <- ScaleData(object = d3, vars.to.regress = c("nUMI", "percent.mito"))
d5 <- ScaleData(object = d5, vars.to.regress = c("nUMI", "percent.mito"))
d6 <- ScaleData(object = d6, vars.to.regress = c("nUMI", "percent.mito"))


# 7. Dimensionality Reduction: PCA
# --------------------------------
d1 <- RunPCA(object = d1, pc.genes = d1.vhostgenes, do.print = FALSE)
d2 <- RunPCA(object = d2, pc.genes = d2.vhostgenes, do.print = FALSE)
d3 <- RunPCA(object = d3, pc.genes = d3.vhostgenes, do.print = FALSE)
d5 <- RunPCA(object = d5, pc.genes = d5.vhostgenes, do.print = FALSE)
d6 <- RunPCA(object = d6, pc.genes = d6.vhostgenes, do.print = FALSE)

# scoring all gene (including non-var genes) in the data based on their correlation with the calculated components.
d1 <- ProjectPCA(object = d1, do.print = FALSE) 
d2 <- ProjectPCA(object = d2, do.print = FALSE)
d3 <- ProjectPCA(object = d3, do.print = FALSE)
d5 <- ProjectPCA(object = d5, do.print = FALSE)
d6 <- ProjectPCA(object = d6, do.print = FALSE)


# 8. Determine statistically significant principal components
# -----------------------------------------------------------
d1 <- JackStraw(object = d1, num.replicate = 100)
d2 <- JackStraw(object = d2, num.replicate = 100)
d3 <- JackStraw(object = d3, num.replicate = 100)
d5 <- JackStraw(object = d5, num.replicate = 100)
d6 <- JackStraw(object = d6, num.replicate = 100)

JackStrawPlot(object = d1, PCs = 1:20)
JackStrawPlot(object = d2, PCs = 1:20)
JackStrawPlot(object = d3, PCs = 1:20)
JackStrawPlot(object = d5, PCs = 1:20)
JackStrawPlot(object = d6, PCs = 1:20)

PCElbowPlot(object = d1)
PCElbowPlot(object = d2)
PCElbowPlot(object = d3)
PCElbowPlot(object = d5)
PCElbowPlot(object = d6)


# 9. Cluster the cells (phenograph method and SNN-cliq)
# -----------------------------------------------------
d1 <- FindClusters(object = d1, reduction.type = "pca", dims.use = 1:16, resolution = 0.4, print.output = 0, force.recalc = TRUE, save.SNN = TRUE)
d2 <- FindClusters(object = d2, reduction.type = "pca", dims.use = 1:16, resolution = 0.4, print.output = 0, force.recalc = TRUE, save.SNN = TRUE) 
d3 <- FindClusters(object = d3, reduction.type = "pca", dims.use = 1:16, resolution = 0.4, print.output = 0, force.recalc = TRUE, save.SNN = TRUE) 
d5 <- FindClusters(object = d5, reduction.type = "pca", dims.use = 1:16, resolution = 0.4, print.output = 0, force.recalc = TRUE, save.SNN = TRUE)
d6 <- FindClusters(object = d6, reduction.type = "pca", dims.use = 1:16, resolution = 0.4, print.output = 0, force.recalc = TRUE, save.SNN = TRUE) 

# to check parameters:
#PrintFindClustersParams(object = d)


# 10. non-linear dimensional reduction
# ------------------------------------
d1 <- RunTSNE(object = d1, dims.use = 1:16, do.fast = TRUE)
d2 <- RunTSNE(object = d2, dims.use = 1:16, do.fast = TRUE)
d3 <- RunTSNE(object = d3, dims.use = 1:16, do.fast = TRUE)
d5 <- RunTSNE(object = d5, dims.use = 1:16, do.fast = TRUE)
d6 <- RunTSNE(object = d6, dims.use = 1:16, do.fast = TRUE)

t1 <- TSNEPlot(object = d1, do.label = TRUE, do.return = TRUE, pt.size = 0.5)
t2 <- TSNEPlot(object = d2, do.label = TRUE, do.return = TRUE, pt.size = 0.5)
t3 <- TSNEPlot(object = d3, do.label = TRUE, do.return = TRUE, pt.size = 0.5)
t5 <- TSNEPlot(object = d5, do.label = TRUE, do.return = TRUE, pt.size = 0.5)
t6 <- TSNEPlot(object = d6, do.label = TRUE, do.return = TRUE, pt.size = 0.5)

plot_grid(t1, t2, t3, t5, t6, ncol=3)

dev.print(device=pdf, file="lane1/plot_grid.pdf")


# 11. Finding differentially expressed genes (cluster biomarkers)
# ---------------------------------------------------------------
# find clusters example:
# find all markers of cluster 1
#cluster1.markers <- FindMarkers(object = d, ident.1 = 1, min.pct = 0.25)

# find all markers distinguishing cluster 5 from clusters 0 and 3
#cluster5.markers <- FindMarkers(object = d, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)

# find markers for every cluster compared to all remaining cells, report
# only the positive ones
#d.markers <- FindAllMarkers(object = d, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

d1.all.markers <- FindAllMarkers(d1, only.pos=FALSE, min.pct=0.25, thresh.use=0.25)
d2.all.markers <- FindAllMarkers(d2, only.pos=FALSE, min.pct=0.25, thresh.use=0.25)
d3.all.markers <- FindAllMarkers(d3, only.pos=FALSE, min.pct=0.25, thresh.use=0.25)
d5.all.markers <- FindAllMarkers(d5, only.pos=FALSE, min.pct=0.25, thresh.use=0.25)
d6.all.markers <- FindAllMarkers(d6, only.pos=FALSE, min.pct=0.25, thresh.use=0.25)

write.table(d1.all.markers, file="lane1/s1_allmarkers.xls", sep="\t", quote=F)
write.table(d2.all.markers, file="lane1/s2_allmarkers.xls", sep="\t", quote=F)
write.table(d3.all.markers, file="lane1/s3_allmarkers.xls", sep="\t", quote=F)
write.table(d5.all.markers, file="lane1/s5_allmarkers.xls", sep="\t", quote=F)
write.table(d6.all.markers, file="lane1/s6_allmarkers.xls", sep="\t", quote=F)


# 12. Plot Gene Markers:
# ----------------------
# Gene markers: "IFNB","IFNL1","IFNL2","IFNL3","ACE2","TMPRSS2","IL6","IL1B","IL18","CXCL10","IFITM3","LY6E"

key.genes1 <- c("IFNB","IFNL1","IFNL2","IFNL3","ACE2","TMPRSS2","IL6","IL1B","IL18","CXCL10","IFITM3","LY6E")

# plotting violin plot
svg(filename="lane1/Proj033_s1-lane1-keygenes1-vol.svg", width=16, height=12, pointsize=1)
VlnPlot(d1, features = key.genes1, point.size.use=0)
dev.off()
svg(filename="lane1/Proj033_s2-lane1-keygenes1-vol.svg", width=16, height=12, pointsize=1)
VlnPlot(d2, features = key.genes1, point.size.use=0)
dev.off()
svg(filename="lane1/Proj033_s3-lane1-keygenes1-vol.svg", width=16, height=12, pointsize=1)
VlnPlot(d3, features = key.genes1, point.size.use=0)
dev.off()
svg(filename="lane1/Proj033_s5-lane1-keygenes1-vol.svg", width=16, height=12, pointsize=1)
VlnPlot(d5, features = key.genes1, point.size.use=0)
dev.off()
svg(filename="lane1/Proj033_s6-lane1-keygenes1-vol.svg", width=16, height=12, pointsize=1)
VlnPlot(d6, features = key.genes1, point.size.use=0)
dev.off()

# plotting feature plot
svg(filename="lane1/Proj033_s1-lane1-keygenes1.svg", width=16, height=12, pointsize=1)
FeaturePlot(d1, features = key.genes1, cols=c("grey", "red"), min.cutoff="q9")
dev.off()
svg(filename="lane1/Proj033_s2-lane1-keygenes1.svg", width=16, height=12, pointsize=1)
FeaturePlot(d2, features = key.genes1, cols=c("grey", "red"), min.cutoff="q9")
dev.off()
svg(filename="lane1/Proj033_s3-lane1-keygenes1.svg", width=16, height=12, pointsize=1)
FeaturePlot(d3, features = key.genes1, cols=c("grey", "red"), min.cutoff="q9")
dev.off()
svg(filename="lane1/Proj033_s5-lane1-keygenes1.svg", width=16, height=12, pointsize=1)
FeaturePlot(d5, features = key.genes1, cols=c("grey", "red"), min.cutoff="q9")
dev.off()
svg(filename="lane1/Proj033_s6-lane1-keygenes1.svg", width=16, height=12, pointsize=1)
FeaturePlot(d6, features = key.genes1, cols=c("grey", "red"), min.cutoff="q9")
dev.off()

# note: for sample 1
key.genes1 <- c("ACE2","TMPRSS2","IL6","IL1B","IL18","IFITM3","LY6E") #missing "IFNB"; no value "IFNL1","IFNL2","IFNL3","CXCL10",

# note:for sample 6
key.genes1 <- c("ACE2","TMPRSS2","IL6","IL1B","IL18","CXCL10","IFITM3","LY6E")#no value,"IFNL1", "IFNL2","IFNL3",
