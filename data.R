#load the library

library(Seurat)

#file 1 rename old name to new
list.files("C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/GSE132771_RAW/NML 1")
file.rename(
  "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/GSE132771_RAW/NML 1/genes.tsv.gz",
  "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/GSE132771_RAW/NML 1/features.tsv.gz"
)

#load the data
NML_I <- Read10X(data.dir = "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/GSE132771_RAW/NML 1")


#file 2 rename 

file.rename(
  "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/GSE132771_RAW/NML 2/genes.tsv.gz",
  "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/GSE132771_RAW/NML 2/features.tsv.gz"
)

#load data

NML_II <- Read10X(data.dir = "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/GSE132771_RAW/NML 2")

#file 3 rename no need
#load data 
NML_III <- Read10X(data.dir = "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/GSE132771_RAW/NML 3")

NML_I
dim(NML_I)
dim(NML_II)
dim(NML_III)

# Create Seurat objects
library(tidyverse)
NML_I_obj <- CreateSeuratObject(
  counts = NML_I,
  project = "Normal_I",
  min.cells = 3,
  min.features = 200
)
class(NML_I)
class(NML_I_obj) 

NML_II_obj <- CreateSeuratObject(
  counts = NML_II,
  project = "Normal_II",
  min.cells = 3,
  min.features = 200
)

class(NML_II)
class(NML_II_obj)

NML_III_obj <- CreateSeuratObject(
  counts = NML_III,
  project = "Normal_III",
  min.cells = 3,
  min.features = 200
)

class(NML_III)
class(NML_III_obj)
# View all three
NML_I
colnames(NML_I[])
rownames(NML_I[])
View(NML_I)
View(NML_I_obj@meta.data)

NML_II
colnames(NML_II[])
rownames(NML_II[])
View(NML_II)
View(NML_II_obj@meta.data)

NML_III
colnames(NML_III[])
rownames(NML_III[])
View(NML_III)
View(NML_III_obj@meta.data)

#save

output_dir <- "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/RDS_objects/"
dir.create(output_dir, showWarnings = FALSE)

saveRDS(NML_I_obj,   file = paste0(output_dir, "NML_I_obj.rds"))
saveRDS(NML_II_obj,  file = paste0(output_dir, "NML_II_obj.rds"))
saveRDS(NML_III_obj, file = paste0(output_dir, "NML_III_obj.rds"))

#Read seurat
NML_I_obj <- readRDS("C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/RDS_objects/NML_I_obj.rds")
NML_II_obj <- readRDS("C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/RDS_objects/NML_I_obj.rds")
NML_III_obj <- readRDS("C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/RDS_objects/NML_I_obj.rds")

#repeat

Normalrepeat_I <- readRDS("C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/RDS_objects/NML_I_obj.rds")
rm(Normalrepeat_I)

output_dir <- "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/RDS_objects/"

NML_I_obj   <- readRDS(paste0(output_dir, "NML_I_obj.rds"))
NML_II_obj  <- readRDS(paste0(output_dir, "NML_II_obj.rds"))
NML_III_obj <- readRDS(paste0(output_dir, "NML_III_obj.rds"))

NML_II_obj
NML_III_obj
#Merge seurat objects

MergedNML <- merge(
  NML_I_obj,
  y = c(NML_II_obj, NML_III_obj),
  add.cell.ids = c("NML_I", "NML_II", "NML_III"),
  project = "MergedNML"
)

ls()
MergedNML
view(MergedNML@meta.data)

saveRDS(MergedNML, file = paste0(output_dir, "MergedNML.rds"))

#preprocessing workflow qc quaity control
range(MergedNML$nFeature_RNA)
range(MergedNML$nCount_RNA)
MergedNML <- PercentageFeatureSet(MergedNML, pattern = "^MT-", col.name = "percent.mt")
view(MergedNML@meta.data)
range(MergedNML$percent.mt)

##### Selecting cells for further analysis
# Visualize QC metrics as a violin plot

VlnPlot(MergedNML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

MergedNML

MergedNML <- subset(MergedNML, 
                    subset = nFeature_RNA > 500 & 
                      nFeature_RNA < 4000 & 
                      nCount_RNA < 20000 & 
                      percent.mt < 10
)
MergedNML
VlnPlot(MergedNML, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
MergedNML
saveRDS(MergedNML, file = paste0(output_dir, "MergedNML_filtered.rds"))

#data normalization, feature selection, data scaling
view(MergedNML)

LayerData(MergedNML, layer = "counts")@x
MergedNML <- NormalizeData(MergedNML,
                           normalization.method = "LogNormalize",
                           scale.factor = 10000)
LayerData(MergedNML, layer = "data")@x

MergedNML <- FindVariableFeatures(MergedNML,
                                  selection.method = "vst",
                                  nfeatures = 2000)

VariableFeatures(MergedNML)

# Visualize variable features
plot1 <- VariableFeaturePlot(MergedNML)
# Label the top 10 most variable genes on the plot
top10 <- head(VariableFeatures(MergedNML), 10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

plot2

#scaling
MergedNML <- ScaleData(MergedNML)

LayerData(MergedNML, layer = "scale.data")[1:5, 1:5]
MergedNML@commands

saveRDS(MergedNML, file = paste0(output_dir, "MergedNML_scaled.rds"))
#PCA on scaled data for dimensional reduction
MergedNML <- RunPCA(MergedNML)
DimPlot(MergedNML, reduction = "pca", dims = c(1,2))
DimPlot(MergedNML, reduction = "pca", dims = c(1,10))
DimPlot(MergedNML, reduction = "pca", dims = c(1,50))

MergedNML <- JackStraw(MergedNML, num.replicate = 100)
MergedNML <- ScoreJackStraw(MergedNML, dims = 1:20)
JackStrawPlot(MergedNML, dims = 1:20)
ElbowPlot(MergedNML)
ElbowPlot(MergedNML, ndims = 50, reduction = "pca")

####### Cluster the cells
MergedNML <- FindNeighbors(MergedNML, dims = 1:20)
MergedNML <- FindClusters(MergedNML, resolution = 0.1)
MergedNML <- FindClusters(MergedNML, resolution = 0.3)

####### Run non-linear dimensional reduction (UMAP/tSNE)
MergedNML <- RunUMAP(MergedNML, dims = 1:20)
DimPlot(MergedNML, reduction = "umap", label = TRUE, repel = TRUE)

MergedNML <- RunTSNE(object = MergedNML)
DimPlot(object = MergedNML, reduction = "tsne")

View(MergedNML)
saveRDS(MergedNML, file = paste0(output_dir, "MergedNML_clustered.rds"))

DimPlot(MergedNML, reduction = "umap", label = TRUE)
View(MergedNML@meta.data)

MergedNML.list <- SplitObject(MergedNML, split.by = 'orig.ident')

MergedNML.list
length(MergedNML.list)

names(MergedNML.list)

MergedNML.list <- lapply(X = MergedNML.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = MergedNML.list)
length(features)
head(features, 20)

MergedNML.anchors <- FindIntegrationAnchors(
  object.list = MergedNML.list,
  anchor.features = features
)
MergedNML.integrated <- IntegrateData(anchorset = MergedNML.anchors)
DefaultAssay(MergedNML.integrated) <- "integrated"
rm(MergedNML.anchors)
rm(MergedNML.list)
gc()
ls()

library(Seurat)
library(tidyverse)
library(harmony)

output_dir <- "C:/Users/aliaa/OneDrive/Pictures/Camera Roll/Single_cell transcriptomics/RDS_objects/"
MergedNML <- readRDS(paste0(output_dir, "MergedNML_clustered.rds"))

MergedNML
ls()

DimPlot(MergedNML, reduction = "umap", label = TRUE)

DimPlot(MergedNML, reduction = "umap", group.by = "orig.ident")

DimPlot(MergedNML, reduction = "umap", label = TRUE, group.by = "seurat_clusters")

FeaturePlot(MergedNML,
            features = c("EPCAM", "CLDN5", "COL1A2", "PTPRC"),
            cols = c("lightgrey", "blue"))


plot1 <- DimPlot(MergedNML, 
                 reduction = "umap", 
                 group.by = "orig.ident")
plot2 <- DimPlot(MergedNML, 
                 reduction = "umap", 
                 group.by = "seurat_clusters",
                 label = TRUE)
plot1 + plot2

ggsave(
  filename = paste0(output_dir, "UMAP_comparison.png"),
  width = 14,
  height = 6
)

DefaultAssay(MergedNML) <- "RNA"

MergedNML <- JoinLayers(MergedNML)

Layers(MergedNML)
markers <- FindAllMarkers(MergedNML,
                          only.pos = TRUE,
                          min.pct = 0.25,
                          logfc.threshold = 0.25)
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
View(top5)

markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  print(n = 90)

new.cluster.ids <- c(
  "Monocytes",                    # 0
  "Alveolar Macrophages",         # 1
  "Club Cells",                   # 2
  "Endothelial Cells",            # 3
  "AT2 Cells",                    # 4
  "Fibroblasts",                  # 5
  "Interstitial Macrophages",     # 6
  "T Cells",                      # 7
  "NK Cells",                     # 8
  "Smooth Muscle Cells",          # 9
  "Venous Endothelial Cells",     # 10
  "Transitional Macrophages",     # 11
  "AT1 Cells",                    # 12
  "Plasma B Cells",               # 13
  "Mast Cells",                   # 14
  "Proliferating Cells",          # 15
  "Lymphatic Endothelial Cells",  # 16
  "Ciliated Cells"                # 17
)
names(new.cluster.ids) <- levels(MergedNML)

MergedNML <- RenameIdents(MergedNML, new.cluster.ids)

DimPlot(MergedNML, 
        reduction = "umap", 
        label = TRUE, 
        pt.size = 0.5,
        repel = TRUE) + 
  NoLegend()

saveRDS(MergedNML, file = paste0(output_dir, "MergedNML_annotated.rds"))
DimPlot(MergedNML,
        reduction = "umap",
        label = TRUE,
        label.size = 4,
        pt.size = 0.5,
        repel = TRUE) +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, size = 14)
  ) +
  labs(title = "Normal Human Lung scRNA-seq\n11,612 cells | 18 cell types")
# Save UMAP
ggsave(paste0(output_dir, "UMAP_annotated.pdf"),
       width = 10, height = 8, dpi = 300)


