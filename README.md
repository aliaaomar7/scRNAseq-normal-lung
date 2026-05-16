# scRNAseq-normal-lung
# Normal Human Lung scRNA-seq | Learning Project

I built this as my first single-cell RNA-seq analysis — working through the full Seurat pipeline from raw 10X data to annotated cell types.
The dataset is normal human lung tissue from three donors (GSE132771), NML_I, NML_II, NML_III.

---

## What I did

Started with raw 10X Chromium data, filtered low-quality cells, 
normalized and clustered, then manually annotated 18 cell populations 
using canonical lung markers.

Raw 10X → QC → Normalization → PCA → UMAP → Clustering → Annotation

---

## Dataset

GSE132771 — 3 normal human lung donors, 10X Chromium droplet-seq
~11,612 cells after QC filtering

https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132771

---

## QC filters I used

- Genes per cell: 500–4,000
- UMI counts: < 20,000  
- Mitochondrial %: < 10%

---

## Cell types identified (18 clusters)

Monocytes, Alveolar Macrophages, Club Cells, Endothelial Cells,
AT2 Cells, Fibroblasts, Interstitial Macrophages, T Cells, NK Cells,
Smooth Muscle Cells, Venous Endothelial Cells, Transitional Macrophages,
AT1 Cells, Plasma B Cells, Mast Cells, Proliferating Cells,
Lymphatic Endothelial Cells, Ciliated Cells

---

## Tools

R, Seurat v5, tidyverse, harmony

---
---

Aliaa | PharmD Candidate & Computational Drug Discovery
