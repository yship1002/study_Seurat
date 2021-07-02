library(dplyr)
library(Seurat)
library(patchwork)
pdf("./1.pdf", height = 6000, width = 10000)
# Load the PBMC dataset
G6F.data <- Read10X(data.dir = "/lustre/project/wdeng7/jyang10/data/sampleG6F/outs/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
G6F <- CreateSeuratObject(counts = G6F.data, project = "G6F", min.cells = 3, min.features = 200)
G6F
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
G6F[["percent.mt"]] <- PercentageFeatureSet(G6F, pattern = "^mt")
# Visualize QC metrics as a violin plot
VlnPlot(G6F, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(G6F, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(G6F, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
G6F <- subset(G6F, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
G6F <- NormalizeData(G6F, normalization.method = "LogNormalize", scale.factor = 10000)
G6F <- FindVariableFeatures(G6F, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(G6F), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(G6F)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(G6F)
G6F <- ScaleData(G6F, features = all.genes)
G6F <- RunPCA(G6F, features = VariableFeatures(object = G6F))
VizDimLoadings(G6F, dims = 1:2, reduction = "pca")
DimPlot(G6F, reduction = "pca")
DimHeatmap(G6F, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(G6F, dims = 1:15, cells = 500, balanced = TRUE)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time
G6F <- JackStraw(G6F, num.replicate = 100)
G6F <- ScoreJackStraw(G6F, dims = 1:20)
JackStrawPlot(G6F, dims = 1:15)
ElbowPlot(G6F)

G6F <- FindNeighbors(G6F, dims = 1:10)
G6F <- FindClusters(G6F, resolution = 0.5)
# Look at cluster IDs of the first 5 cells
head(Idents(G6F), 5)
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
G6F <- RunUMAP(G6F, dims = 1:10)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(G6F, reduction = "umap")
saveRDS(G6F, file = "/home/jyang10/G6F_tutorial.rds")

dev.off()
