# R Script using Seurat (v2.3.4) for CoV2 Viral Cells Analysis
# Extension from Hybrid Human CoV2 Pipeline (hybrid_pipeline_seurat_R_script.R)
# Created by Prof. Steve Shen and organized by Hezkiel Nanda
# Last update: 29-Aug-2020

# 1. Preliminary check
# --------------------
# check for how many cells contain virus reads after filtering
table(Matrix::colSums(d1@data[cov2.gene.names, ])>0)
# FALSE 
#   851 
table(Matrix::colSums(d3@data[cov2.gene.names, ])>0)
# FALSE  TRUE 
#  1447    10 
 
# 2. CoV2 virus infected cells information
# ----------------------------------------
# for virus reads > 0 using a temporary variable
temp1 <- SubsetData(d1, Matrix::colSums(d1@data[cov2.gene.names, ]) !=0)
temp2 <- SubsetData(d2, Matrix::colSums(d2@data[cov2.gene.names, ]) !=0)
temp3 <- SubsetData(d3, Matrix::colSums(d3@data[cov2.gene.names, ]) !=0)
temp5 <- SubsetData(d5, Matrix::colSums(d5@data[cov2.gene.names, ]) !=0)
temp6 <- SubsetData(d6, Matrix::colSums(d6@data[cov2.gene.names, ]) !=0)

dim(temp1@data)
# [1] 33276     2
dim(temp2@data)
# [1] 33276  1981
dim(temp3@data)
# [1] 33276   335
dim(temp5@data)
# [1] 33276  1723
dim(temp6@data)
# [1] 33276    27


# summary
#				S1		S2 		S3 		S5 		S6
# at 0.0000001	2		1981	335		1723	27
# at 0.000001	2		1981	335		1723	27
# at 0.00001	2		1980	312		1723	26
# at 0.0001		0		1195	193		1712	7
# at 0.001		0		469		118		566		0
# at 0.01		0		391		7		463		0
# at 0.1		0		105		0		333		0


# for virus at perc accept.low=0.0000001
lowerthreshold <- 0.0000001
temp1 <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
temp2 <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
temp3 <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
temp5 <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
temp6 <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

dim(temp1@data)
# [1] 33276     2
dim(temp2@data)
# [1] 33276  1981
dim(temp3@data)
# [1] 33276   335
dim(temp5@data)
# [1] 33276  1723
dim(temp6@data)
# [1] 33276    27

# for virus at perc accept.low=0.000001
lowerthreshold <- 0.000001
temp1 <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
temp2 <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
temp3 <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
temp5 <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
temp6 <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

dim(temp1@data)
# [1] 33276     2
dim(temp2@data)
# [1] 33276  1981
dim(temp3@data)
# [1] 33276   335
dim(temp5@data)
# [1] 33276  1723
dim(temp6@data)
# [1] 33276    27

# for virus at perc accept.low=0.00001
lowerthreshold <- 0.00001
temp1 <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
temp2 <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
temp3 <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
temp5 <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
temp6 <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

dim(temp1@data)
# [1] 33276     2
dim(temp2@data)
# [1] 33276  1980
dim(temp3@data)
# [1] 33276   312
dim(temp5@data)
# [1] 33276  1723
dim(temp6@data)
# [1] 33276    26

# for virus at perc accept.low=0.0001
lowerthreshold <- 0.0001
temp1 <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
temp2 <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
temp3 <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
temp5 <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
temp6 <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

dim(temp1@data)
# NULL
dim(temp2@data)
# [1] 33276  1195
dim(temp3@data)
# [1] 33276   193
dim(temp5@data)
# [1] 33276  1712
dim(temp6@data)
# [1] 33276     7

# for virus at perc accept.low=0.001
lowerthreshold <- 0.001
temp1 <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
temp2 <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
temp3 <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
temp5 <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
temp6 <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

dim(temp1@data)
# NULL
dim(temp2@data)
# [1] 33276   469
dim(temp3@data)
# [1] 33276   118
dim(temp5@data)
# [1] 33276   566
dim(temp6@data)
# NULL

# for virus at perc accept.low=0.01
lowerthreshold <- 0.01
temp1 <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
temp2 <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
temp3 <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
temp5 <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
temp6 <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

dim(temp1@data)
# NULL
dim(temp2@data)
# [1] 33276   391
dim(temp3@data)
# [1] 33276     7
dim(temp5@data)
# [1] 33276   463
dim(temp6@data)
# NULL

# for virus at perc accept.low=0.1
lowerthreshold <- 0.1
temp1 <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
temp2 <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
temp3 <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
temp5 <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
temp6 <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

dim(temp1@data)
# NULL
dim(temp2@data)
# [1] 33276   105
dim(temp3@data)
# [1] 33276     0
dim(temp5@data)
# [1] 33276   333
dim(temp6@data)
# [1] 33276     0


# 3. Plot Virus Cell
# -------------------
# summary percentages:
# at0.00001		0.177	93.093	19.642	99.217	99.948	1.536	0.000	1.902	0.00001
# at0.0001		0.000	57.237	11.438	80.331	99.424	0.320	0.000	0.207	0.0001
# at0.001		0.000	22.921	5.604	24.456	35.288	0.128	0.000	0.000	0.001

# at 0.00001
lowerthreshold <- 0.00001
d1x <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
d2x <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
d3x <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
d5x <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
d6x <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

# at 0.0001
lowerthreshold <- 0.0001
d1y <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
d2y <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
d3y <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
d5y <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
d6y <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

# at 0.001
lowerthreshold <- 0.001
d1z <- SubsetData(d1, subset.name="percent.cov2", accept.low=lowerthreshold)
d2z <- SubsetData(d2, subset.name="percent.cov2", accept.low=lowerthreshold)
d3z <- SubsetData(d3, subset.name="percent.cov2", accept.low=lowerthreshold)
d5z <- SubsetData(d5, subset.name="percent.cov2", accept.low=lowerthreshold)
d6z <- SubsetData(d6, subset.name="percent.cov2", accept.low=lowerthreshold)

# Plot with TSNE Dimensionality reduction with proper cells highlight
p1x <- DimPlot(d1, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d1x@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p2x <- DimPlot(d2, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d2x@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p3x <- DimPlot(d3, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d3x@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p5x <- DimPlot(d5, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d5x@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p6x <- DimPlot(d6, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d6x@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)

plot_grid(p1x, p2x, p3x, p5x, p6x, ncol=3)

dev.print(device=pdf, file="lane1/TSNE_virus_cells@0.00001.pdf")

p1y <- DimPlot(d1, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d1y@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p2y <- DimPlot(d2, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d2y@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p3y <- DimPlot(d3, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d3y@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p5y <- DimPlot(d5, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d5y@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p6y <- DimPlot(d6, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d6y@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)

plot_grid(p1y, p2y, p3y, p5y, p6y, ncol=3)

dev.print(device=pdf, file="lane1/TSNE_virus_cells@0.0001.pdf")

p1z <- DimPlot(d1, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d1z@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5) #Error: Insufficient values in manual scale. 5 needed but only 1 provided.
p2z <- DimPlot(d2, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d2z@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p3z <- DimPlot(d3, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d3z@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p5z <- DimPlot(d5, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d5z@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)
p6z <- DimPlot(d6, reduction.use = "tsne", pt.size = 1, cols.use = "lightgrey", cells.highlight = d6z@cell.names, cols.highlight = "cadetblue", sizes.highlight = 0.5)

plot_grid(p1z, p2z, p3z, p5z, p6z, ncol=3)

dev.print(device=pdf, file="lane1/TSNE_virus_cells@0.001.pdf")

# plotting virus side by side with the TSNE Plot from hybrid_pipeline_seurat_R_script.R
plot_grid(t1, p1x, ncol=2)
dev.print(device=pdf, file="lane1/S1_TSNE_virus_cells.pdf")

plot_grid(t2, p2x, p2z, ncol=3)
dev.print(device=pdf, file="lane1/S2_TSNE_virus_cells.pdf")

plot_grid(t3, p3x, p3z, ncol=3)
dev.print(device=pdf, file="lane1/S3_TSNE_virus_cells.pdf")

plot_grid(t5, p5x, p5z, ncol=3)
dev.print(device=pdf, file="lane1/S5_TSNE_virus_cells.pdf")

plot_grid(t6, p6x, p6z, ncol=3)
dev.print(device=pdf, file="lane1/S6_TSNE_virus_cells.pdf")