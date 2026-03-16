
# load package ------------------------------------------------------------

.libPaths("/path/to/your/R/library")

library(Seurat); options(Seurat.object.assay.version = "v5"); options(future.globals.maxSize = 1e9)
library(harmony)
library(SeuratWrappers)

library(tidyverse)
library(magrittr)
library(data.table)
library(glue)

library(RColorBrewer); library(ggsci); library(scales)
library(ComplexHeatmap); library(circlize)
library(ggrepel)
library(ggrastr)
library(patchwork)

setwd("/path/to/your/project/directory")

# preprocess --------------------------------------------------------------

# HRA003613 seurat object was provided by author
allso <- readRDS("HRA003613.PBMC.seurat_obj.rds")
allso <- readRDS("HRA003613.SG.seurat_obj.rds")

# GSE157278 cellranger output was processed at the next part
allso <- readRDS("pSS_GSE157278.PBMC.allso.rds")
allso@meta.data <- readRDS("pSS_GSE157278.PBMC.metadata.rds")

allso %<>% NormalizeData()

soList <- SplitObject(allso, "sample")
names(soList)

sceList <- purrr::map(soList, ~ SingleCellExperiment::SingleCellExperiment(assays = list(counts = .x[["RNA"]]$counts, logcounts = .x[["RNA"]]$data)))
dblDen <- purrr::map(sceList, ~ scDblFinder::computeDoubletDensity(.x, dims = 20))
dblDen <- purrr::map(dblDen, ~ log10(.x + 1))
soList <- purrr::map2(soList, dblDen, ~ {.x$dbl <- .y; .x})

dblVal <- map(soList, ~ .x$dbl) %>% unname %>% do.call(c, .)
allso$dbl <- dblVal[colnames(allso)]

allso %<>% FindVariableFeatures(nfeatures = 2000)
allso %<>% ScaleData()

allso %<>% RunPCA()
allso %>% ElbowPlot(ndims = 50)
allso %<>% FindNeighbors(dims = 1:20)
allso %<>% FindClusters(resolution = 2)
allso %<>% RunUMAP(dims = 1:20, min.dist = .5, n.neighbors = 50)

allso %<>% RunHarmony("sample")
allso %>% ElbowPlot(ndims = 50, reduction = "harmony")
allso %<>% RunUMAP(dims = 1:30, reduction = "harmony", min.dist = .5, n.neighbors = 50)
allso %<>% FindNeighbors(dims = 1:30, reduction = "harmony")
allso %<>% FindClusters(resolution = 2)

allso %<>% Azimuth::RunAzimuth(reference = "Azimuth_Human-PBMC/")

saveRDS(allso, "pSS.PBMC.allso.rds")
saveRDS(allso@meta.data, "pSS.PBMC.metadata.rds")
allso <- readRDS("pSS.PBMC.allso.rds")
allso@meta.data <- readRDS("pSS.PBMC.metadata.rds")

saveRDS(allso, "pSS.SG.allso.rds")
saveRDS(allso@meta.data, "pSS.SG.metadata.rds")
allso <- readRDS("pSS.SG.allso.rds")
allso@meta.data <- readRDS("pSS.SG.metadata.rds")

# GSE157278 ---------------------------------------------------------------

sampDir <- list.dirs("2_cellranger", recursive = F)
sampName <- str_remove(sampDir, ".*/")

dataList <- map(
  sampDir,
  ~ {
    message("Reading ", .x)
    # Read10X(str_c(.x, "/filtered_feature_bc_matrix"))
    Read10X(str_c(.x, "/raw_feature_bc_matrix"))
  }
)

names(dataList) <- sampName
dataList <- imap(dataList, ~ {set_colnames(.x, str_c(.y, "_", colnames(.x)))})
colnames(dataList[[1]]) %>% head
map_int(dataList, ncol)
map_int(dataList, nrow)

soList <- map(dataList, CreateSeuratObject)
map_int(soList, ncol)

soList[[1]]@meta.data %>% colnames()
# soList <- map(soList, ~ {.x$nFeature_RNA <- colSums(.x[["RNA"]]$counts != 0); .x})
# soList <- map(soList, ~ {.x$nCount_RNA <- colSums(.x[["RNA"]]$counts); .x})

ppGeneSet <- list(
  mt = c(
    "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6", "MT-CO3",
    "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB"),
  malat = "MALAT1",
  xist = "XIST"
)

soList <- purrr::map(soList, ~ calcPercentageGeneSets(.x, geneSets = ppGeneSet))
# soList <- purrr::map(soList, ~ {.x[, .x$nFeature_RNA > 200 & .x$nCount_RNA > 500 & .x$mt < 20]})
soList <- purrr::map(soList, ~ {.x[, .x$nFeature_RNA > 100 & .x$nCount_RNA > 200 & .x$mt < 20]})
soList <- purrr::map(soList, NormalizeData)
soList <- purrr::map(soList, ~ FindVariableFeatures(.x))

sceList <- purrr::map(soList, ~ SingleCellExperiment::SingleCellExperiment(assays = list(counts = .x[["RNA"]]$counts, logcounts = .x[["RNA"]]$data)))
# varGeneList <- purrr::map(soList, VariableFeatures)
# dblDen <- purrr::map2(sceList, varGeneList, ~ scDblFinder::computeDoubletDensity(.x, subset.row = intersect(rownames(.x), .y), dims = 20))
dblDen <- purrr::map(sceList, ~ scDblFinder::computeDoubletDensity(.x, dims = 20))
dblDen <- purrr::map(dblDen, ~ log10(.x + 1))
soList <- purrr::map2(soList, dblDen, ~ {.x$dbl <- .y; .x})

soList <- purrr::map(soList, ScaleData)
soList <- purrr::map(soList, RunPCA)
soList <- purrr::map(soList, ~ FindNeighbors(.x, dims = 1:20))
soList <- purrr::map(soList, ~ FindClusters(.x, resolution = 1))
soList <- purrr::map(soList, ~ RunUMAP(.x, dims = 1:20))

saveRDS(soList, "soList_raw.rds")
soList <- readRDS("soList_raw.rds")

names(soList)

nowSamp <- "HC-1"
nowSamp <- "HC-2"
nowSamp <- "HC-3"
nowSamp <- "HC-4"
nowSamp <- "HC-5"
nowSamp <- "pSS-1"
nowSamp <- "pSS-2"
nowSamp <- "pSS-3"
nowSamp <- "pSS-4"
nowSamp <- "pSS-5"

so <- soList[[nowSamp]]

so$size <- so$nCount_RNA %>% log10()
VlnPlot(so, c("size", "nFeature_RNA", "mt", "dbl", "malat", "xist"), group.by = "seurat_clusters")
DimPlot(so, group.by = "seurat_clusters", label = T)

so$rm <- "no"
so$rm[so$dbl > 0.8] <- "yes"
table(so$rm)
DimPlot(so, group.by = "rm")
so <- so[, so$rm == "no"]

badClu <- c() # Name

so$rm[so$seurat_clusters %in% badClu] <- "yes"
table(so$rm)
DimPlot(so, group.by = "rm")
so <- so[, so$rm == "no"]

cMarker <- FindMarkers(so, ident.1 = c(16), min.pct = .25, only.pos = T)

FeaturePlot(so, "PTPRC", label = T)
FeaturePlot(so, "EPCAM", label = T)
FeaturePlot(so, "KRT19", label = T)
FeaturePlot(so, "WFDC2", label = T)
FeaturePlot(so, "COL1A2", label = T)
FeaturePlot(so, "PECAM1", label = T)
FeaturePlot(so, "RBFOX1", label = T)
FeaturePlot(so, "TOP2A", label = T)
# Ery
FeaturePlot(so, "HBA1", label = T)
# Platelet
FeaturePlot(so, "PPBP", label = T)
FeaturePlot(so, "GATA1", label = T)
# Neu
FeaturePlot(so, "MPO", label = T)
FeaturePlot(so, "DEFA4", label = T)
FeaturePlot(so, "LTF", label = T)
FeaturePlot(so, "AFF3", label = T)
FeaturePlot(so, "MMP9", label = T)
FeaturePlot(so, "FCGR3B", label = T)
FeaturePlot(so, "S100A8", label = T)
FeaturePlot(so, "G0S2", label = T)
FeaturePlot(so, "ISG15", label = T)
FeaturePlot(so, "PI3", label = T)
FeaturePlot(so, "IL1R2", label = T)
# HSPC
FeaturePlot(so, "CD34", label = T)
FeaturePlot(so, "SPINK2", label = T)
FeaturePlot(so, "CLC", label = T)
FeaturePlot(so, "MS4A2", label = T)
FeaturePlot(so, "TOP2A", label = T)
FeaturePlot(so, "MKI67", label = T)
# T
FeaturePlot(so, "CD3D", label = T)
FeaturePlot(so, "IL7R", label = T)
FeaturePlot(so, "CD8A", label = T)
# NK
FeaturePlot(so, "NKG7", label = T)
FeaturePlot(so, "NCAM1", label = T)
FeaturePlot(so, "KLRC2", label = T)
# B
FeaturePlot(so, "CD79A", label = T)
FeaturePlot(so, "MS4A1", label = T)
FeaturePlot(so, "MZB1", label = T)
FeaturePlot(so, "IGHA2", label = T)
FeaturePlot(so, "IGHG1", label = T)
FeaturePlot(so, "JCHAIN", label = T)
# Mono/Mac
FeaturePlot(so, "FCN1", label = T)
FeaturePlot(so, "FABP4", label = T)
FeaturePlot(so, "SPP1", label = T)
FeaturePlot(so, "CD14", label = T)
FeaturePlot(so, "LRMDA", label = T)
FeaturePlot(so, "FCGR3A", label = T)
FeaturePlot(so, "CD68", label = T)
FeaturePlot(so, "C1QA", label = T)
# DC
FeaturePlot(so, "CD74", label = T)
FeaturePlot(so, "CD1C", label = T)
FeaturePlot(so, "LILRA4", label = T)

soList[[nowSamp]] <- so

saveRDS(soList, "soList_qc.rds")
soList <- readRDS("soList_qc.rds")

countList <- map(soList, ~ .x[["RNA"]]$counts)
dim(countList[[1]])
dim(countList[[2]])

countData <- do.call(cbind, countList)
allso <- CreateSeuratObject(countData)
allso$sample <- str_remove(colnames(allso), "_.*")
allso$sample %>% table

allso$sampType <- str_remove(allso$sample, "-.")
allso@meta.data[, c("sample", "sampType")] %>% table

map(soList, ~ .x@meta.data) %>% map_int(ncol)
allMeta <- map(soList, ~ .x@meta.data) %>% do.call(rbind, .)
rownames(allMeta) <- colnames(allso)
allso@meta.data <- cbind(allso@meta.data, allMeta[, c("size", "mt", "dbl", "malat", "xist")])
allso %<>% NormalizeData()

saveRDS(allso, "pSS_GSE157278.PBMC.allso.rds")
saveRDS(allso@meta.data, "pSS_GSE157278.PBMC.metadata.rds")
allso <- readRDS("pSS_GSE157278.PBMC.allso.rds")
allso@meta.data <- readRDS("pSS_GSE157278.PBMC.metadata.rds")

# anno --------------------------------------------------------------------

allso$predicted.celltype.l1 %>% table %>% sort

so <- allso[, allso$predicted.celltype.l1 %in% c("CD8 T", "CD4 T", "NK", "other T")]
so <- allso[, allso$predicted.celltype.l1 %in% c("B")]
so <- allso[, allso$predicted.celltype.l1 %in% c("Mono", "DC")]
so <- allso[, allso$predicted.celltype.l1 %in% c("other")]

so %<>% FindVariableFeatures(nfeatures = 2000)
so %<>% ScaleData()

so %<>% RunPCA()
so %>% ElbowPlot(ndims = 50)
so %<>% FindNeighbors(dims = 1:20)
so %<>% FindClusters(resolution = 2)
so %<>% RunUMAP(dims = 1:20, min.dist = .5, n.neighbors = 50)

so %<>% RunHarmony("sample")
so %>% ElbowPlot(ndims = 50, reduction = "harmony")
so %<>% RunUMAP(dims = 1:20, reduction = "harmony", min.dist = .5, n.neighbors = 50)
so %<>% FindNeighbors(dims = 1:20, reduction = "harmony")
so %<>% FindClusters(resolution = 2)

DimPlot(so, group.by = "predicted.celltype.l1", label = T)
DimPlot(so, group.by = "predicted.celltype.l2", label = T)
DimPlot(so, group.by = "seurat_clusters", label = T)

# so$nFeature_RNA <- colSums(so[["RNA"]]$counts != 0)
# so$nCount_RNA <- colSums(so[["RNA"]]$counts)

so$size <- so$nCount_RNA %>% log10()
VlnPlot(so, c("size", "nFeature_RNA", "mt", "dbl", "malat", "xist"), group.by = "seurat_clusters")

so$celltype <- "CD8 Memory"
so$celltype[so$seurat_clusters %in% c()] <- "dbl"
so$celltype[so$seurat_clusters %in% c()] <- "no"
so$celltype[so$seurat_clusters %in% c()] <- "CD4 Naive"
so$celltype[so$seurat_clusters %in% c()] <- "CD4 Memory"
so$celltype[so$seurat_clusters %in% c()] <- "CD4 Treg"
so$celltype[so$seurat_clusters %in% c()] <- "CD4 CTL"
so$celltype[so$seurat_clusters %in% c()] <- "CD8 Naive"
so$celltype[so$seurat_clusters %in% c()] <- "CD8 Memory"
so$celltype[so$seurat_clusters %in% c()] <- "MAIT"
so$celltype[so$seurat_clusters %in% c()] <- "gdT"
so$celltype[so$seurat_clusters %in% c()] <- "NKT"
so$celltype[so$seurat_clusters %in% c()] <- "NK CD56lo"
so$celltype[so$seurat_clusters %in% c()] <- "NK CD56hi"
so$celltype[so$seurat_clusters %in% c()] <- "Prolif T"
so$celltype[so$seurat_clusters %in% c()] <- "Prolif NK"
so$celltype[so$predicted.celltype.l2 == "ILC"] <- "ILC"

Idents(so) <- so$seurat_clusters
cMarker <- FindMarkers(so, ident.1 = c(), min.pct = .25, only.pos = T)

FeaturePlot(so, "PTPRC", label = T)
FeaturePlot(so, "EPCAM", label = T)
FeaturePlot(so, "KRT19", label = T)
FeaturePlot(so, "WFDC2", label = T)
FeaturePlot(so, "PECAM1", label = T)
FeaturePlot(so, "COL1A1", label = T)
FeaturePlot(so, "COL1A2", label = T)
FeaturePlot(so, "DCN", label = T)
FeaturePlot(so, "ACTA2", label = T)
FeaturePlot(so, "ACTC1", label = T)
FeaturePlot(so, "MYL9", label = T)
FeaturePlot(so, "TOP2A", label = T)
FeaturePlot(so, "PCNA", label = T)
FeaturePlot(so, "MKI67", label = T)
## Red Blood Cell
FeaturePlot(so, "HBA1", label = T)
FeaturePlot(so, "HBB", label = T)
## Platelet
FeaturePlot(so, "PPBP", label = T)
FeaturePlot(so, "PF4", label = T)
## Neutrophil
FeaturePlot(so, "CD24", label = T)
FeaturePlot(so, "LTF", label = T)
FeaturePlot(so, "FCGR3B", label = T)
## T
FeaturePlot(so, "CD2", label = T)
FeaturePlot(so, "CD3D", label = T)
FeaturePlot(so, "CD8A", label = T)
FeaturePlot(so, "CD4", label = T)
# T naive
FeaturePlot(so, "CCR7", label = T)
FeaturePlot(so, "LEF1", label = T)
FeaturePlot(so, "TCF7", label = T)
# CD8 TEM
FeaturePlot(so, "CCL5", label = T)
FeaturePlot(so, "GZMH", label = T)
FeaturePlot(so, "GZMK", label = T)
FeaturePlot(so, "GZMB", label = T)
# CD4 TCM
FeaturePlot(so, "ITGB1", label = T)
FeaturePlot(so, "LTB", label = T)
FeaturePlot(so, "CD27", label = T)
# CD4 TEM
FeaturePlot(so, "IL7R", label = T)
FeaturePlot(so, "KLRG1", label = T)
# Th17 & Th2
FeaturePlot(so, "GATA3", label = T)
FeaturePlot(so, "CCR6", label = T)
FeaturePlot(so, "CTSH", label = T)
# Tph
FeaturePlot(so, "CXCL13", label = T)
FeaturePlot(so, "ADGRG1", label = T)
FeaturePlot(so, "PRDM1", label = T)
FeaturePlot(so, "HLA-DRA", label = T)
FeaturePlot(so, "HLA-DRB1", label = T)
# Tfh
FeaturePlot(so, "CXCR5", label = T)
FeaturePlot(so, "PDCD1", label = T)
# gdT
FeaturePlot(so, "TRDC", label = T)
FeaturePlot(so, "TRDV1", label = T)
FeaturePlot(so, "TRDV2", label = T)
# Treg
FeaturePlot(so, "FOXP3", label = T)
FeaturePlot(so, "RTKN2", label = T)
# MAIT
FeaturePlot(so, "KLRB1", label = T)
FeaturePlot(so, "SLC4A10", label = T)
## NK
FeaturePlot(so, "NKG7", label = T)
FeaturePlot(so, "GNLY", label = T)
FeaturePlot(so, "NCR1", label = T)
FeaturePlot(so, "NCAM1", label = T)
FeaturePlot(so, "KLRC2", label = T)
FeaturePlot(so, "XCL1", label = T)
# ILC
FeaturePlot(so, "KIT", label = T)
FeaturePlot(so, "HDC", label = T)
## B
FeaturePlot(so, "CD79A", label = T)
FeaturePlot(so, "CD79B", label = T)
FeaturePlot(so, "MS4A1", label = T)
# imm B
FeaturePlot(so, "EBF1", label = T)
FeaturePlot(so, "AFF3", label = T)
# Plasmablast
FeaturePlot(so, "IGHG2", label = T)
FeaturePlot(so, "IGHA2", label = T)
FeaturePlot(so, "MZB1", label = T)
FeaturePlot(so, "JCHAIN", label = T)
# B memory
FeaturePlot(so, "COCH", label = T)
FeaturePlot(so, "AIM2", label = T)
FeaturePlot(so, "CD27", label = T)
# B naive
FeaturePlot(so, "IGHM", label = T)
FeaturePlot(so, "IGHD", label = T)
# B intermediate
FeaturePlot(so, "TNFRSF13B", label = T)
# DNB2
FeaturePlot(so, "ZEB2", label = T)
FeaturePlot(so, "TBX21", label = T)
FeaturePlot(so, "ITGAX", label = T)
# DNB1
FeaturePlot(so, "HOPX", label = T)
FeaturePlot(so, "CXCR5", label = T)
## Mono
FeaturePlot(so, "CTSS", label = T)
FeaturePlot(so, "FCN1", label = T)
# CD14 Mono
FeaturePlot(so, "CD14", label = T)
FeaturePlot(so, "S100A9", label = T)
FeaturePlot(so, "CD83", label = T)
FeaturePlot(so, "LRMDA", label = T)
# CD16 mono
FeaturePlot(so, "FCGR3A", label = T)
FeaturePlot(so, "CDKN1C", label = T)
# Macrophage
FeaturePlot(so, "LYZ", label = T)
FeaturePlot(so, "CD68", label = T)
FeaturePlot(so, "CCL2", label = T)
FeaturePlot(so, "CCR2", label = T)
FeaturePlot(so, "FABP4", label = T)
FeaturePlot(so, "SPP1", label = T)
FeaturePlot(so, "MERTK", label = T)
FeaturePlot(so, "FOLR2", label = T)
FeaturePlot(so, "CCL3", label = T)
FeaturePlot(so, "CCL4", label = T)
FeaturePlot(so, "C1QA", label = T)
FeaturePlot(so, "C1QB", label = T)
FeaturePlot(so, "C1QC", label = T)
## DC
FeaturePlot(so, "CD74", label = T)
FeaturePlot(so, "HLA-DPA1", label = T)
# ASDC
FeaturePlot(so, "AXL", label = T)
# pDC
FeaturePlot(so, "ITM2C", label = T)
FeaturePlot(so, "PLD4", label = T)
FeaturePlot(so, "SERPINF1", label = T)
FeaturePlot(so, "CLEC4C", label = T)
FeaturePlot(so, "LILRA4", label = T)
FeaturePlot(so, "IL3RA", label = T)
FeaturePlot(so, "PTPRS", label = T)
# cDC-CD1C+
FeaturePlot(so, "CLEC9A", label = T) # cDC1
FeaturePlot(so, "CLEC10A", label = T) # cDC2
FeaturePlot(so, "FCER1A", label = T)
FeaturePlot(so, "CD1C", label = T)
# Mig-DC
FeaturePlot(so, "FSCN1", label = T)
## HSPC
FeaturePlot(so, "CD34", label = T)
FeaturePlot(so, "SPINK2", label = T)
FeaturePlot(so, "CD69", label = T)
FeaturePlot(so, "CLC", label = T)
FeaturePlot(so, "GATA2", label = T)
FeaturePlot(so, "PBX1", label = T)
# CLP
FeaturePlot(so, "ACY3", label = T)
FeaturePlot(so, "PRSS2", label = T)
# LMPP
FeaturePlot(so, "AVP", label = T)
FeaturePlot(so, "CRHBP", label = T)
# EMP
FeaturePlot(so, "MYCT1", label = T)
FeaturePlot(so, "GATA1", label = T)
FeaturePlot(so, "KLF1", label = T)
FeaturePlot(so, "GATA2", label = T)
# GMP
FeaturePlot(so, "SERPINB10", label = T)
FeaturePlot(so, "RETN", label = T)
# MPO+ Neutrophil
FeaturePlot(so, "TOP2A", label = T)
FeaturePlot(so, "CD24", label = T)
FeaturePlot(so, "MPO", label = T)
FeaturePlot(so, "ELANE", label = T)
FeaturePlot(so, "AFF2", label = T)
# DEFA3+ Neutrophil
FeaturePlot(so, "DEFA4", label = T)
FeaturePlot(so, "DEFA3", label = T)
FeaturePlot(so, "MS4A3", label = T)
# LTF+ Neutrophil
FeaturePlot(so, "LTF", label = T)
FeaturePlot(so, "CEACAM1", label = T)
FeaturePlot(so, "ITGAM", label = T)
# FCGR3B+ Neutrophil
FeaturePlot(so, "CSF3R", label = T)
FeaturePlot(so, "FCGR3B", label = T)
FeaturePlot(so, "S100A8", label = T)
FeaturePlot(so, "G0S2", label = T)
FeaturePlot(so, "MMP9", label = T)
FeaturePlot(so, "MME", label = T)
FeaturePlot(so, "IL1R2", label = T)
# Eos
FeaturePlot(so, "SERPINB2", label = T)
FeaturePlot(so, "GATA2", label = T)
# Mast
FeaturePlot(so, "RGS13", label = T)

allso$celltype <- ""
allso@meta.data[colnames(so), "celltype"] <- so$celltype

# imm share ---------------------------------------------------------------

pbmc <- readRDS("pSS.PBMC.allso_clear.rds")
pbmc@meta.data <- readRDS("pSS.PBMC.metadata_clear.rds")

pbmc_pub <- readRDS("pSS_GSE157278.PBMC.allso.rds")
pbmc_pub@meta.data <- readRDS("pSS_GSE157278.PBMC.metadata.rds")

pbmc_pub <- pbmc_pub[, pbmc_pub$celltype %ni% c("Mega", "Ery")]

sgimm <- allso[, allso$celltype %in% c(
  "CD4 Memory", "CD4 Naive", "CD8 Memory GZMB", "Prolif T",
  "B Naive", "B Memory", "B Inte", "Plasma",
  "MDM SPP1", "MDM FCN1", "IM", "Prolif Mac", "Mast"
)]
pbmc$data <- "PBMC"
pbmc_pub$data <- "PBMC_pub"
sgimm$data <- "SG"

pbmc$sample %>% table
sgimm$sample %>% table
pbmc_pub$sample %>% table

countList <- list(
  pbmc[["RNA"]]$counts,
  pbmc_pub[["RNA"]]$counts,
  sgimm[["RNA"]]$counts
)

countData <- do.call(cbind, countList)

imm <- CreateSeuratObject(countData)
imm$sample <- c(
  pbmc$sample %>% unname,
  pbmc_pub$sample %>% unname,
  sgimm$sample %>% unname
)

imm$sampType <- c(
  pbmc$sampType %>% unname,
  pbmc_pub$sampType %>% unname,
  sgimm$sampType %>% unname
)
imm$sampType %>% table

imm$data <- c(
  pbmc$data %>% unname,
  pbmc_pub$data %>% unname,
  sgimm$data %>% unname
)
imm$data %>% table

imm$celltype <- c(
  pbmc$celltype %>% unname,
  pbmc_pub$celltype %>% unname,
  sgimm$celltype %>% unname
)

imm$dbl <- c(
  pbmc$dbl %>% unname,
  pbmc_pub$dbl %>% unname,
  sgimm$dbl %>% unname
)

imm %<>% calcPercentageGeneSets(geneSets = ppGeneSet)
imm %<>% NormalizeData()

imm %<>% FindVariableFeatures(nfeatures = 3000)
imm %<>% ScaleData()

imm %<>% RunPCA()
imm %>% ElbowPlot(ndims = 50)
imm %<>% FindNeighbors(dims = 1:20)
imm %<>% FindClusters(resolution = 2)
imm %<>% RunUMAP(dims = 1:20, min.dist = .5, n.neighbors = 50)

imm %<>% RunHarmony("sample")
imm %>% ElbowPlot(ndims = 50, reduction = "harmony")
imm %<>% RunUMAP(dims = 1:20, reduction = "harmony", min.dist = .5, n.neighbors = 50)
imm %<>% FindNeighbors(dims = 1:20, reduction = "harmony")
imm %<>% FindClusters(resolution = 2)

imm %<>% Azimuth::RunAzimuth(reference = "Azimuth_Human-PBMC/")

saveRDS(imm, "pSS.imm.allso_clear.rds")
saveRDS(imm@meta.data, "pSS.imm.metadata_clear.rds")
imm <- readRDS("pSS.imm.allso_clear.rds")
imm@meta.data <- readRDS("pSS.imm.metadata_clear.rds")

imm@meta.data[, c("celltype", "data")] %>% table

allso$celltype_adj %>% table
allso@meta.data[, c("celltype_adj", "sampType")] %>% table

# celltype_l1 -------------------------------------------------------------

allso$celltype_l1 <- ""

allso$celltype_l1[allso$celltype_adj %in% c(
  "CD4 Naive", "CD4 Th1", "CD4 Th2", "CD4 preTfh", "CD4 Tfh", "CD4 Treg", "CD4 CTL"
)] <- "CD4 T"
allso$celltype_l1[allso$celltype_adj %in% c(
  "CD8 Naive", "CD8 Memory GZMK", "CD8 Memory GZMB"
)] <- "CD8 T"
allso$celltype_l1[allso$celltype_adj %in% c(
  "dnT", "gdT", "MAIT", "Prolif T"
)] <- "T other"
allso$celltype_l1[allso$celltype_adj %in% c(
  "NK", "NK Memory", "NK CD56hi", "ILC", "Prolif NK"
)] <- "NK"
allso$celltype_l1[allso$celltype_adj %in% c(
  "B Naive", "B Memory", "Plasma"
)] <- "B"
allso$celltype_l1[allso$celltype_adj %in% c(
  "cDC1", "cDC2", "pDC", "ASDC"
)] <- "DC"
allso$celltype_l1[allso$celltype_adj %in% c(
  "CD14 Mono", "CD16 Mono", "Mac FOLR2", "Mac TREM2", "Prolif Mac"
)] <- "Mono/Mac"
allso$celltype_l1[allso$celltype_adj %in% c(
  "HSPC", "BaEoMa", "Mast"
)] <- "HSPC/Gran"

allso$celltype_l1[allso$celltype_adj %in% c(
  "Basal", "Duct", "Serous LTF", "Serous LPO", "Serous PCDH10", "Mucous", "Iono", "Myoepi"
)] <- "Epi"
allso$celltype_l1[allso$celltype_adj %in% c(
  "EC Art", "EC Ven"
)] <- "Endo"
allso$celltype_l1[allso$celltype_adj %in% c(
  "Fib Act", "Fib APOE", "Fib APOE stress", "Fib F2R", "Fib F2R stress", "Fib PCOLCE2", "Peri", "Muscle"
)] <- "Str"
allso$celltype_l1[allso$celltype_adj %in% c(
  "Neuron", "Melano"
)] <- "Other"

allso@meta.data[, c("celltype_adj", "celltype_l1")] %>% table

imm$celltype_l1 <- ""

imm$celltype_l1[imm$celltype %in% c(
  "CD4 Naive", "CD4 Th1", "CD4 Th2", "CD4 preTfh", "CD4 Tfh", "CD4 Treg", "CD4 CTL"
)] <- "CD4 T"
imm$celltype_l1[imm$celltype %in% c(
  "CD8 Naive", "CD8 Memory GZMK", "CD8 Memory GZMB"
)] <- "CD8 T"
imm$celltype_l1[imm$celltype %in% c(
  "dnT", "gdT", "MAIT", "Prolif T"
)] <- "T other"
imm$celltype_l1[imm$celltype %in% c(
  "NK", "NK Memory", "NK CD56hi", "ILC", "Prolif NK"
)] <- "NK"
imm$celltype_l1[imm$celltype %in% c(
  "B Naive", "B Memory", "Plasma"
)] <- "B"
imm$celltype_l1[imm$celltype %in% c(
  "cDC1", "cDC2", "pDC", "ASDC"
)] <- "DC"
imm$celltype_l1[imm$celltype %in% c(
  "CD14 Mono", "CD16 Mono", "Mac FOLR2", "Mac TREM2", "Prolif Mac"
)] <- "Mono/Mac"
imm$celltype_l1[imm$celltype %in% c(
  "HSPC", "BaEoMa", "Mast"
)] <- "HSPC/Gran"

imm$celltype_l1[imm$celltype %in% c(
  "Basal", "Duct", "Serous LTF", "Serous LPO", "Serous PCDH10", "Mucous", "Iono", "Myoepi"
)] <- "Epi"
imm$celltype_l1[imm$celltype %in% c(
  "EC Art", "EC Ven"
)] <- "Endo"
imm$celltype_l1[imm$celltype %in% c(
  "Fib Act", "Fib APOE", "Fib APOE stress", "Fib F2R", "Fib F2R stress", "Fib PCOLCE2", "Peri", "Muscle"
)] <- "Str"
imm$celltype_l1[imm$celltype %in% c(
  "Neuron", "Melano"
)] <- "Other"

imm@meta.data[, c("celltype", "celltype_l1")] %>% table

# color -------------------------------------------------------------------

allso$sampType %>% table

groupCol <- c(
  "#0077B6", "#e54a00"
) %>% set_names(
  c(
    "HC", "pSS"
  )
)

allso$celltype_l1 %>% table

ctl1Col <- c(
  "#D22E2E", "#F47B00", "#FF9F00", "#FABF2C", "#7a1fa1",
  "#303f9f", "#1976d2", "#0097a6",
  "#388d3b", "#679f38", "#aeb32b", "#606060"
) %>% set_names(
  c(
    "CD4 T", "CD8 T", "T other", "NK", "B",
    "Mono/Mac", "DC", "HSPC/Gran",
    "Endo", "Str", "Epi",
    "Other"
  )
)

fibCol <- c(
  "#f9a620", "#c51b7d", "#f1b6da", "#4d9221", "#b8e186",
  "#ffd449", "#762a83", "#00b4d8"
) %>% set_names(
  c(
    "Fib Act", "Fib APOE stress", "Fib F2R stress",
    "Fib APOE", "Fib F2R",
    "Fib PCOLCE2",
    "Peri", "Muscle"
  )
)
fibCol %>% show_col()

pss$tls %>% table
tlsCol <- c(
  "#00b4d8", "#ff5400", "#ff6d00", "#ff8500", "#ff9100", "#ff9e00"
) %>% set_names(
  c(
    "no", "TLS1", "TLS2", "TLS3", "TLS4", "TLS5"
  )
)

tlsclsCol <- c(
  "#20854E", "#0072B5", "#E18727", "#BC3C29", "gray70"
) %>% set_names(
  c(
    "VEC & Stromal (T1)",
    "Plasma enrich (T2)",
    "B & T cell enrich (T3)",
    "Duct niche (T4)",
    "Other"
  )
)

st$type %>% unique %>% sort
clsCol <- c(
  "#224c25", "#4e6441", "#8bb174", "#b5ca8d", "#e3dbc3",
  "#f5bd1f", "#ffd93d", "#ffe863",
  "#caaded","#a266c7", "#784b93", "#4e305f"
) %>% set_names(
  c(
    "Serous acini (E1)", "Serous acini (E2)", "Serous acini (E3)", "Serous acini (E4)", "Serous acini (E5)",
    "Mucous acini (E6)", "Mucous acini (E7)", "Mucous acini (E8)",
    "VEC & Stromal (T1)",
    "Plasma enrich (T2)",
    "B & T cell enrich (T3)",
    "Duct niche (T4)"
  )
)
clsCol %>% show_col()

typeCol <- c(
  "#caaded","gray80", "#784b93", "#4e305f"
) %>% set_names(
  c(
    "v", "p", "b", "d"
  )
)

cls_groupCol <- rep(
  c("#2364aa", "#ea7317"), c(8, 4)
) %>% set_names(
  c(
    "Serous acini (E1)", "Serous acini (E2)", "Serous acini (E3)", "Serous acini (E4)", "Serous acini (E5)",
    "Mucous acini (E6)", "Mucous acini (E7)", "Mucous acini (E8)",
    "VEC & Stromal (T1)",
    "Plasma enrich (T2)",
    "B & T cell enrich (T3)",
    "Duct niche (T4)"
  )
)
cls_groupCol %>% show_col()

# pct plot ----------------------------------------------------------------

imm$data %>% table
imm$paper <- ""
imm$paper[imm$data == "PBMC"] <- "PBMC (Xiang. et al.)"
imm$paper[imm$data == "PBMC_pub"] <- "PBMC (Hong. et al.)"
imm$paper[imm$data == "SG"] <- "SG (Xiang. et al.)"

metaData <- imm@meta.data %>% as.data.table()
Ndata <- metaData[, .N, by = .(paper, celltype_l1, sampType)]

Ndata[, pct := N / sum(N) * 100, by = .(paper, sampType)]
Ndata$celltype_l1 %<>% factor(levels = names(ctl1Col))

ggplot(Ndata, aes(x = paper, y = pct)) +
  geom_col(aes(fill = celltype_l1), width = .8) +
  facet_wrap(~ sampType, nrow = 1) +
  scale_fill_manual(values = ctl1Col) +
  scale_y_continuous(expand = expansion(c(0, 0))) +
  labs(x = "", y = "Celltype Percentage (%)", fill = "") +
  theme(
    aspect.ratio = 2/1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key = element_blank(),
    strip.background = element_blank()
  )
ggsave("imm.pct.col.pdf", width = 8, height = 5)


metaData <- allso@meta.data %>% as.data.table()
Ndata <- metaData[, .N, by = .(celltype_l1, sampType)]

Ndata[, pct := N / sum(N) * 100, by = .(sampType)]
Ndata$celltype_l1 %<>% factor(levels = names(ctl1Col))

ggplot(Ndata, aes(x = sampType, y = pct)) +
  geom_col(aes(fill = celltype_l1), width = .8) +
  scale_fill_manual(values = ctl1Col) +
  scale_y_continuous(expand = expansion(c(0, 0))) +
  labs(x = "", y = "Celltype Percentage (%)", fill = "") +
  theme(
    aspect.ratio = 2/1,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key = element_blank(),
    strip.background = element_blank()
  )
ggsave("SG.pct.col.pdf", width = 8, height = 5)


metaData <- imm@meta.data %>% as.data.table()
Ndata <- metaData[data != "SG", .N, by = .(paper, celltype, celltype_l1, sampType)]

Ndata[, pct := N / sum(N) * 100, by = .(paper, sampType)]
# Ndata$celltype_l1 %<>% factor(levels = names(ctl1Col))

pctData <- dcast(Ndata[, -5], paper + celltype + celltype_l1 ~ sampType, value.var = "pct")
pctData[is.na(pctData)] <- 0
pctData[, lfc := log2(pSS + 0.1) - log2(HC + 0.1)]
pctData$celltype_l1 %<>% factor(levels = names(ctl1Col))

ggplot(pctData, aes(x = fct_reorder(celltype, lfc, .desc = T), y = lfc)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_col(aes(fill = celltype_l1, color = paper), width = .8) +
  scale_fill_manual(values = ctl1Col) +
  scale_color_manual(values = c("gray50", "black")) +
  labs(x = "", y = "Log2 (pSS / HC)", title = "Celltype Enrichment", fill = "", color = "") +
  guides(color = guide_legend(override.aes = list(fill = NA))) +
  theme(
    aspect.ratio = 1/3,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key = element_blank()
  )
ggsave("PBMC.celltype_enrich_l2fc.col.pdf", width = 12, height = 4)


metaData <- allso@meta.data %>% as.data.table()
Ndata <- metaData[, .N, by = .(celltype_adj, celltype_l1, sampType)]

Ndata[, pct := N / sum(N) * 100, by = .(sampType)]

pctData <- dcast(Ndata[, -4], celltype_adj + celltype_l1 ~ sampType, value.var = "pct")
pctData[is.na(pctData)] <- 0
pctData[, lfc := log2(pSS + 0.1) - log2(HC + 0.1)]
pctData$celltype_l1 %<>% factor(levels = names(ctl1Col))

ggplot(pctData, aes(x = fct_reorder(celltype_adj, lfc, .desc = T), y = lfc)) +
  geom_hline(yintercept = 0, color = "black") +
  geom_col(aes(fill = celltype_l1), color = "black", width = .8) +
  scale_fill_manual(values = ctl1Col) +
  # scale_color_manual(values = c("gray50", "black")) +
  labs(x = "", y = "Log2 (pSS / HC)", title = "Celltype Enrichment", fill = "", color = "") +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  theme(
    aspect.ratio = 1/5,
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.grid.major.y = element_line(color = "gray90"),
    panel.grid.major.x = element_line(color = "gray90"),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.key = element_blank()
  )
ggsave("SG.celltype_enrich_l2fc.col.pdf", width = 18, height = 4)

# GO ----------------------------------------------------------------------

library(clusterProfiler)
library(org.Hs.eg.db)

egoList <- purrr::map(geneList, ~ {
  enrichGO(
    gene = na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = .x, columns = "ENTREZID", keytype = "SYMBOL")$ENTREZID),
    OrgDb = "org.Hs.eg.db", ont = "BP", pvalueCutoff = 1, qvalueCutoff = 1, readable = T
  )
})
names(egoList)

iwalk(egoList, ~ write.csv(.x@result, str_c("GO/", .y, ".GO.csv")))

# CCC ---------------------------------------------------------------------

library(CellChat)

sub <- allso[, allso$sampType == "pSS"]
sub <- allso[, allso$sampType == "HC"]
sub$celltype_adj %>% table %>% sort

keepType <- sub$celltype_adj %>% table %>% {names(.)[. > 10]}

cco <- createCellChat(
  object = sub[["RNA"]]$data[, sub$celltype_adj %in% keepType],
  meta = sub@meta.data[sub$celltype_adj %in% keepType, ], group.by = "celltype_adj")

CellChatDB <- CellChatDB.human
cco@DB <- CellChatDB

cco <- subsetData(cco)
# future::plan("multisession", workers = 4)
# future::plan("multisession", workers = 1)
cco <- identifyOverExpressedGenes(cco)
cco <- identifyOverExpressedInteractions(cco)

cco <- smoothData(cco, adj = PPI.human)

# Compute the communication probability and infer cellular communication network
options(future.globals.maxSize = 2 * 1e9)
cco <- computeCommunProb(cco, type = "truncatedMean", trim = .1, raw.use = F)
cco <- filterCommunication(cco, min.cells = 10)

# Infer the cell-cell communication at a signaling pathway level
cco <- computeCommunProbPathway(cco)

# Calculate the aggregated cell-cell communication network
cco <- aggregateNet(cco)

saveRDS(cco, "pSS_SG_scRNA.cco.rds")
saveRDS(cco, "HC_SG_scRNA.cco.rds")

# CytoSig -----------------------------------------------------------------

expMat <- so[["RNA"]]$data %>% as.matrix
expMat %>% dim
rowm <- rowMeans(expMat)
expMat <- expMat - rowm

write.table(expMat, "cytosig/expMat.csv", sep = ",")

"CytoSig_run.py -i expMat.csv -o out -r 1000 -a 10000 -e 0 -s 2"

cytoZ <- fread("cytosig/out.Zscore")

Zmat <- cytoZ %>% {set_rownames(as.matrix(.[, -1]), .[[1]])}
dim(Zmat)
so[["cytoZ"]] <- CreateAssay5Object(data = Zmat[, colnames(so)])

# to scanpy ---------------------------------------------------------------

# perform tangram analysis in python following the official tutorial

library(reticulate)
use_python("~/app/anaconda3/envs/tangram/bin/python", required = T)
py_config()

sc <- import("scanpy")

so <- allso[, allso$sampType == "pSS"]
so <- allso[, allso$sampType == "HC"]

cellMeta <- so@meta.data[, c("celltype_adj", "sample", "sampType")]
geneMeta <- data.frame(gene = rownames(so))
rownames(geneMeta) <- geneMeta$gene

adata <- sc$AnnData(
  X = as(Matrix::t(so[["RNA"]]$data), "RsparseMatrix"),
  obs = cellMeta,
  var = geneMeta
)

adata$obsm$update(X_umap = Embeddings(so, "umap"))
adata$obsm$update(X_pca = Embeddings(so, "harmony"))

adata$write_h5ad("SG_pSS.sc.h5ad")
adata$write_h5ad("SG_HC.sc.h5ad")

# adata st
adata <- sc$read_h5ad("SG_pSS.tangram_dc.adata_st.h5ad")
adata <- sc$read_h5ad("SG_HC.tangram_dc.adata_st.h5ad")

spaceCord <- adata$obsm[["spatial"]]
spaceCord[, 2] <- -spaceCord[, 2]
spaceCord %>% plot

pss$x_cord <- spaceCord[, 1]
pss$y_cord <- spaceCord[, 2]
pss$cell_count <- adata$obs$cell_count

hc$x_cord <- spaceCord[, 1]
hc$y_cord <- spaceCord[, 2]
hc$cell_count <- adata$obs$cell_count

ct_pred <- adata$obsm[["tangram_ct_pred"]]

identical(colnames(pss), rownames(ct_pred))
pss[["ct_pred"]] <- CreateAssay5Object(data = t(ct_pred))
DefaultAssay(pss) <- "ct_pred"
VariableFeatures(pss) <- pss[["ct_pred"]]$data %>% rownames()

identical(colnames(hc), rownames(ct_pred))
hc[["ct_pred"]] <- CreateAssay5Object(data = t(ct_pred))
DefaultAssay(hc) <- "ct_pred"
VariableFeatures(hc) <- hc[["ct_pred"]]$data %>% rownames()

pss %<>% ScaleData()
pss %<>% RunPCA(npcs = 6)
pss %>% ElbowPlot(ndims = 6)
pss %<>% RunUMAP(dims = 1:6)
pss %<>% FindNeighbors(dims = 1:6)
pss %<>% FindClusters(resolution = 1)

pss$CD4_all <- ct_pred[, colnames(ct_pred) %in% c("CD4 Naive", "CD4 Th1", "CD4 Th2", "CD4 preTfh", "CD4 Tfh", "CD4 Treg", "CD4 CTL")] %>% rowSums
pss$CD8_all <- ct_pred[, colnames(ct_pred) %in% c("CD8 Naive", "CD8 Memory GZMK", "CD8 Memory GZMB")] %>% rowSums
pss$NK_all <- ct_pred[, colnames(ct_pred) %in% c("NK", "NK Memory", "NK CD56hi", "ILC", "Prolif NK")] %>% rowSums
pss$B_all <- ct_pred[, colnames(ct_pred) %in% c("B Naive", "B Memory", "Plasma")] %>% rowSums
pss$DC_all <- ct_pred[, colnames(ct_pred) %in% c("cDC1", "cDC2", "pDC", "ASDC")] %>% rowSums
pss$Mac_all <- ct_pred[, colnames(ct_pred) %in% c("CD14 Mono", "CD16 Mono", "Mac FOLR2", "Mac TREM2", "Prolif Mac")] %>% rowSums
pss$Endo_all <- ct_pred[, colnames(ct_pred) %in% c("EC Art", "EC Ven")] %>% rowSums
pss$Str_all <- ct_pred[, colnames(ct_pred) %in% c("Fib Act", "Fib APOE", "Fib APOE stress", "Fib F2R", "Fib F2R stress", "Fib PCOLCE2", "Peri", "Muscle")] %>% rowSums
pss$Epi_all <- ct_pred[, colnames(ct_pred) %in% c("Basal", "Duct", "Serous LTF", "Serous LPO", "Serous PCDH10", "Mucous", "Iono", "Myoepi")] %>% rowSums

hc %<>% ScaleData()
hc %<>% RunPCA(npcs = 6)
hc %>% ElbowPlot(ndims = 6)
hc %<>% RunUMAP(dims = 1:6)
hc %<>% FindNeighbors(dims = 1:6)
hc %<>% FindClusters(resolution = 1)

hc$CD4_all <- ct_pred[, colnames(ct_pred) %in% c("CD4 Naive", "CD4 Th1", "CD4 Th2", "CD4 preTfh", "CD4 Tfh", "CD4 Treg", "CD4 CTL")] %>% rowSums
hc$CD8_all <- ct_pred[, colnames(ct_pred) %in% c("CD8 Naive", "CD8 Memory GZMK", "CD8 Memory GZMB")] %>% rowSums
hc$NK_all <- ct_pred[, colnames(ct_pred) %in% c("NK", "NK Memory", "NK CD56hi", "ILC", "Prolif NK")] %>% rowSums
hc$B_all <- ct_pred[, colnames(ct_pred) %in% c("B Naive", "B Memory", "Plasma")] %>% rowSums
hc$DC_all <- ct_pred[, colnames(ct_pred) %in% c("cDC1", "cDC2", "pDC", "ASDC")] %>% rowSums
hc$Mac_all <- ct_pred[, colnames(ct_pred) %in% c("CD14 Mono", "CD16 Mono", "Mac FOLR2", "Mac TREM2", "Prolif Mac")] %>% rowSums
hc$Endo_all <- ct_pred[, colnames(ct_pred) %in% c("EC Art", "EC Ven")] %>% rowSums
hc$Str_all <- ct_pred[, colnames(ct_pred) %in% c("Fib Act", "Fib APOE", "Fib APOE stress", "Fib F2R", "Fib F2R stress", "Fib PCOLCE2", "Peri", "Muscle")] %>% rowSums
hc$Epi_all <- ct_pred[, colnames(ct_pred) %in% c("Basal", "Duct", "Serous LTF", "Serous LPO", "Serous PCDH10", "Mucous", "Iono", "Myoepi")] %>% rowSums

# ad ge
ad_ge <- sc$read_h5ad("SG_pSS.tangram_dc.ad_ge.h5ad")
ad_ge <- sc$read_h5ad("SG_HC.tangram_dc.ad_ge.h5ad")

x <- ad_ge$X %>% t
rownames(x) <- ad_ge$var_names$values %>% toupper() %>% str_replace("ORF", "orf") %>% str_replace("Morf", "MORF")
colnames(x) <- ad_ge$obs_names$values
table(colnames(x) %in% colnames(pss))
table(colnames(x) %in% colnames(hc))

varGene <- ad_ge$uns$training_genes %>% toupper() %>% str_replace("ORF", "orf") %>% str_replace("Morf", "MORF")

pss[["ge"]] <- CreateAssay5Object(data = x[, colnames(pss)])
DefaultAssay(pss) <- "ge"
VariableFeatures(pss) <- varGene

hc[["ge"]] <- CreateAssay5Object(data = x[, colnames(hc)])
DefaultAssay(hc) <- "ge"
VariableFeatures(hc) <- varGene

# segment
ad_seg <- sc$read_h5ad("SG_pSS.tangram_dc.adata_segment.h5ad")
ad_seg <- sc$read_h5ad("SG_HC.tangram_dc.adata_segment.h5ad")

segData <- ad_seg$obs %>% as.data.table()

segData[, spot := str_remove(centroids, "_.*")]
table(segData$spot %in% colnames(pss))
table(segData$spot %in% colnames(hc))

spot2tls <- pss$tls
segData[, tls := spot2tls[spot]]

segData[, celltype_l1 := type2l1[cluster]]
segData$celltype_l1 %>% table
segData$celltype_l1 %<>% factor(levels = names(ctl1Col))

ggplot(segData[sample(1:.N)], aes(x = y, y = -x)) +
  geom_point_rast(aes(color = celltype_l1), shape = 16, size = .5) +
  facet_wrap(~ celltype_l1) +
  scale_color_manual(values = ctl1Col[-12]) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  labs(color = "") +
  coord_fixed() +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.key = element_blank()
  )
ggsave("SG_pSS.segment.celltype_l1.cord.pdf", width = 10, height = 6)
ggsave("SG_HC.segment.celltype_l1.cord.pdf", width = 10, height = 6)

saveRDS(segData, "SG_pSS.segmentData.rds")
segData <- readRDS("SG_pSS.segmentData.rds")

saveRDS(segData, "SG_HC.segmentData.rds")
segData <- readRDS("SG_HC.segmentData.rds")

# space -------------------------------------------------------------------

pss <- Load10X_Spatial("space/pSS/")
hc <- Load10X_Spatial("space/HC/")

SpatialFeaturePlot(pss, features = "nCount_Spatial", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(hc, features = "nCount_Spatial", image.alpha = 0.3, pt.size.factor = 2.2)

DefaultAssay(pss) <- "Spatial"
pss %<>% NormalizeData()
# pss %<>% FindVariableFeatures()

pss %<>% ScaleData()
pss %<>% RunPCA(npcs = 50)
pss %>% ElbowPlot(ndims = 50)
pss %<>% RunUMAP(dims = 1:20)
pss %<>% FindNeighbors(dims = 1:20)
pss %<>% FindClusters(resolution = 1)

plotMeta(pss, "seurat_clusters")
plotMeta(pss, "tls")

# TLS detect step 1
pss$cell_count %>% hist(breaks = 100)
pss$tls <- "out"
pss$tls[pss$cell_count >= 45] <- "in"

posData <- fread("space/pSS/spatial/tissue_positions.csv")
posData <- posData[match(colnames(pss), barcode)]
pss$array_row <- posData$array_row
pss$array_col <- posData$array_col

posData <- fread("space/HC/spatial/tissue_positions.csv")
posData <- posData[match(colnames(hc), barcode)]
hc$array_row <- posData$array_row
hc$array_col <- posData$array_col

SpatialDimPlot(pss, group.by = "seurat_clusters", cols = set_names(pal_d3()(9), 0:8), image.alpha = 0.3, pt.size.factor = 2.2)
SpatialDimPlot(hc, group.by = "seurat_clusters", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialDimPlot(pss, group.by = "tls", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(pss, features = "x_cord", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(pss, features = "y_cord", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(pss, features = "array_row", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(pss, features = "array_col", image.alpha = 0.3, pt.size.factor = 2.2)

SpatialFeaturePlot(pss, features = "PC_1", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(pss, features = "PC_2", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(pss, features = "PC_3", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(pss, features = "PC_4", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(pss, features = "PC_5", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialFeaturePlot(pss, features = "PC_6", image.alpha = 0.3, pt.size.factor = 2.2)

saveRDS(pss, "ST.pSS.rds")
saveRDS(pss@meta.data, "ST.pSS_metadata.rds")
pss <- readRDS("ST.pSS.rds")
pss@meta.data <- readRDS("ST.pSS_metadata.rds")

saveRDS(hc, "ST.HC.rds")
saveRDS(hc@meta.data, "ST.HC_metadata.rds")
hc <- readRDS("ST.HC.rds")
hc@meta.data <- readRDS("ST.HC_metadata.rds")

# detect TLS --------------------------------------------------------------

# step 2
metaData <- pss@meta.data %>% as.data.table(T)
tlsCord <- metaData[tls == "in", .(rn, x_cord, y_cord)]

library(dbscan)
xy <- tlsCord[, 2:3] %>% as.matrix()
knnRes <- dbscan::kNN(xy, k = 3, sort = T)
nnDist <- min(knnRes$dist[, 1])

dbscanRes <- dbscan::dbscan(xy, eps = nnDist*2*0.9, minPts = 3)
tlsCord$cluster <- dbscanRes$cluster

pss$TLScluster <- "no"
pss@meta.data[tlsCord$rn, "TLScluster"] <- tlsCord$cluster %>% as.character()
SpatialDimPlot(pss, group.by = "TLScluster", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialDimPlot(pss, group.by = "tls", image.alpha = 0.3, pt.size.factor = 2.2)

pss$TLScluster %>% table %>% sort

pss$tls <- "no"
pss$tls[pss$TLScluster == 2] <- "TLS1"
pss$tls[pss$TLScluster == 1] <- "TLS2"
pss$tls[pss$TLScluster == 5] <- "TLS3"
pss$tls[pss$TLScluster == 4] <- "TLS4"
pss$tls[pss$TLScluster == 6] <- "TLS5"

# step 3
pss$array_col[pss$tls == "TLS3"] %>% min
pss$array_row[pss$tls == "TLS3"] %>% min
pss$array_row[pss$tls == "TLS3"] %>% max
pss$array_col[pss$tls == "TLS3" & pss$array_row == 27] %>% max
pss$array_col[pss$tls == "TLS3" & pss$array_row == 39] %>% max
pss$array_col[pss$tls == "TLS3" & pss$array_row == 39] %>% min
pss$array_col[pss$tls == "TLS3" & pss$array_row == 38] %>% min
pss$array_col[pss$tls == "TLS3" & pss$array_row == 37] %>% min
pss$array_col[pss$tls == "TLS3" & pss$array_row == 37] %>% max
pss$array_col[pss$tls == "TLS3" & pss$array_row == 35] %>% max
pss$array_row[pss$tls == "TLS2"] %>% max
pss$array_col[pss$tls == "TLS2" & pss$array_row == 48] %>% min
pss$array_col[pss$tls == "TLS1"] %>% min
pss$array_row[pss$tls == "TLS1"] %>% min
pss$array_col[pss$tls == "TLS1" & pss$array_row == "59"] %>% max
pss$array_col[pss$tls == "TLS4"] %>% max

pss$tls[pss$tls == "TLS3" & pss$array_col <= 80] <- "no"
pss$tls[pss$tls == "TLS3" & pss$array_row <= 26] <- "no"
pss$tls[pss$tls == "TLS3" & pss$array_row == 27 & pss$array_col == 107] <- "no"
pss$tls[pss$tls == "TLS3" & pss$array_row >= 40] <- "no"
pss$tls[pss$tls == "TLS3" & pss$array_row == 39 & pss$array_col %in% c(105, 85)] <- "no"
pss$tls[pss$tls == "TLS3" & pss$array_row == 38 & pss$array_col == 84] <- "no"
pss$tls[pss$tls == "TLS3" & pss$array_row == 37 & pss$array_col == 83] <- "no"
pss$tls[pss$tls == "TLS3" & pss$array_row == 37 & pss$array_col == 109] <- "no"
pss$tls[pss$tls == "TLS3" & pss$array_row == 35 & pss$array_col == 111] <- "no"
pss$tls[pss$tls == "TLS2" & pss$array_row >= 51] <- "no"
pss$tls[pss$tls == "TLS2" & pss$array_row == 48 & pss$array_col == 22] <- "no"
pss$tls[pss$tls == "TLS1" & pss$array_col <= 58] <- "no"
pss$tls[pss$tls == "TLS1" & pss$array_row <= 58] <- "no"
pss$tls[pss$tls == "TLS1" & pss$array_row == 59 & pss$array_col == 89] <- "no"
pss$tls[pss$tls == "TLS4" & pss$array_col >= 65] <- "no"

# step4: expand
metaData <- pss@meta.data %>% as.data.table(T)
metaData[, id := str_c(array_row, "_", array_col)]

tlsData <- metaData[tls == "TLS1"]
tlsData <- metaData[tls == "TLS2"]
tlsData <- metaData[tls == "TLS3"]
tlsData <- metaData[tls == "TLS4"]
tlsData <- metaData[tls == "TLS5"]
arrayData <- tlsData[, .(array_row, array_col)]

arrayData1 <- copy(arrayData)
arrayData2 <- copy(arrayData)
arrayData3 <- copy(arrayData)
arrayData4 <- copy(arrayData)
arrayData5 <- copy(arrayData)
arrayData6 <- copy(arrayData)

# top right
arrayData1[, array_col := array_col + 1]
arrayData1[, array_row := array_row + 1]
# top
arrayData2[, array_col := array_col + 2]
# top left
arrayData3[, array_col := array_col + 1]
arrayData3[, array_row := array_row - 1]
# bottom left
arrayData4[, array_col := array_col - 1]
arrayData4[, array_row := array_row - 1]
# bottom
arrayData5[, array_col := array_col - 2]
# bottom right
arrayData6[, array_col := array_col - 1]
arrayData6[, array_row := array_row + 1]

arrayData <- rbind(
  arrayData1,
  arrayData2,
  arrayData3,
  arrayData4,
  arrayData5,
  arrayData6
)
arrayData[, id := str_c(array_row, "_", array_col)]

metaData[id %in% arrayData$id, tls := "TLS1"]
metaData[id %in% arrayData$id, tls := "TLS2"]
metaData[id %in% arrayData$id, tls := "TLS3"]
metaData[id %in% arrayData$id, tls := "TLS4"]
metaData[id %in% arrayData$id, tls := "TLS5"]

ggplot(metaData, aes(x = x_cord, y = y_cord)) +
  geom_point(aes(color = tls))

pss$tls <- metaData$tls

# space combine -----------------------------------------------------------

countData <- cbind(
  pss[["Spatial"]]$counts %>% {set_colnames(., str_c("pSS_", colnames(.)))},
  hc[["Spatial"]]$counts %>% {set_colnames(., str_c("HC_", colnames(.)))}
)

st <- CreateSeuratObject(countData)
st$sampType <- rep(c("pSS", "HC"), c(ncol(pss), ncol(hc)))
st$sampType %>% table
st$tls <- c(pss$tls %>% unname(), rep("no"))

allType <- union(rownames(pss[["ct_pred"]]$data), rownames(hc[["ct_pred"]]$data))
pssMat <- matrix(0, nrow = length(allType), ncol = ncol(pss)) %>% set_rownames(allType) %>% set_colnames(colnames(st)[st$sampType == "pSS"])
pssMat[rownames(pss[["ct_pred"]]$data), ] <- pss[["ct_pred"]]$data %>% as.matrix
hcMat <- matrix(0, nrow = length(allType), ncol = ncol(hc)) %>% set_rownames(allType) %>% set_colnames(colnames(st)[st$sampType == "HC"])
hcMat[rownames(hc[["ct_pred"]]$data), ] <- hc[["ct_pred"]]$data %>% as.matrix

st[["ct_pred"]] <- CreateAssay5Object(data = cbind(pssMat, hcMat))

DefaultAssay(st) <- "ct_pred"
VariableFeatures(st) <- st[["ct_pred"]]$data %>% rownames()

st %<>% ScaleData()
st %<>% RunPCA(npcs = 6)
st %>% ElbowPlot(ndims = 6)
st %<>% RunUMAP(dims = 1:5)
st %<>% FindNeighbors(dims = 1:5)
st %<>% FindClusters(resolution = 1)

totGene <- union(rownames(pss[["ge"]]$data), rownames(hc[["ge"]]$data))
zGene <- totGene %>% setdiff(rownames(pss[["ge"]]$data))
zMat <- matrix(0, nrow = length(zGene), ncol = ncol(pss)) %>% as("dgCMatrix") %>% set_rownames(zGene) %>% set_colnames(colnames(pss))
pssMat <- rbind(pss[["ge"]]$data, zMat)[totGene, ] %>% set_colnames(colnames(st)[st$sampType == "pSS"])
zGene <- totGene %>% setdiff(rownames(hc[["ge"]]$data))
zMat <- matrix(0, nrow = length(zGene), ncol = ncol(hc)) %>% as("dgCMatrix") %>% set_rownames(zGene) %>% set_colnames(colnames(hc))
hcMat <- rbind(hc[["ge"]]$data, zMat)[totGene, ] %>% set_colnames(colnames(st)[st$sampType == "HC"])

st[["ge"]] <- CreateAssay5Object(data = cbind(pssMat, hcMat))

st@meta.data[, c("seurat_clusters", "sampType")] %>% table

pss$cluster <- ""
pss$cluster <- st$seurat_clusters[st$sampType == "pSS"] %>% unname %>% as.vector()
hc$cluster <- st$seurat_clusters[st$sampType == "HC"] %>% unname %>% as.vector()

pss$type <- clus2type[pss$cluster] %>% unname
hc$type <- clus2type[hc$cluster] %>% unname

SpatialDimPlot(pss, group.by = "cluster", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialDimPlot(hc, group.by = "cluster", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialDimPlot(pss, group.by = "type", image.alpha = 0.3, pt.size.factor = 2.2)
SpatialDimPlot(hc, group.by = "type", image.alpha = 0.3, pt.size.factor = 2.2)

st$type <- str_c("C", st$seurat_clusters)
st$type %>% table
st$type[st$seurat_clusters %in% c(10, 4, 1, 7, 0)] %<>% str_c("Serous acini (", ., ")")
st$type[st$seurat_clusters %in% c(2, 5, 9)] %<>% str_c("Mucous acini (", ., ")")
st$type[st$seurat_clusters %in% c(6)] %<>% str_c("VEC & Stromal (", ., ")")
st$type[st$seurat_clusters %in% c(8)] %<>% str_c("Plasma enrich (", ., ")")
st$type[st$seurat_clusters %in% c(3)] %<>% str_c("B & T cell enrich (", ., ")")
st$type[st$seurat_clusters %in% c(11)] %<>% str_c("Duct niche (", ., ")")

# cluster enrich
metaData <- st@meta.data %>% as.data.table(T)
Ndata <- metaData[, .N, by = .(seurat_clusters, type, sampType)]
Ndata[, pct := N / sum(N) * 100, by = sampType]
enrichMat <- Ndata[, .(sampType, seurat_clusters, type, pct)] %>% dcast(seurat_clusters + type ~ sampType, value.var = "pct")
enrichMat[is.na(enrichMat)] <- 0
enrichMat[, lfc := log2(pSS + 0.1) - log2(HC + 0.1)]

# cluster - celltype
Idents(st) <- st$seurat_clusters
allMark <- FindAllMarkers(st, test.use = "t", only.pos = T) %>% as.data.table()

usedType <- allMark[, gene[1:5], by = cluster][["V1"]] %>% na.omit %>% unique
usedType %<>% setdiff(c("NK", "dnT", "Prolif NK", "CD4 preTfh", "pDC", "CD4 Th2", "NK Memory", "Prolif Mac"))

usedMat <- st[["ct_pred"]]$data %>% as.matrix()
usedMat <- apply(usedMat, 2, rescale)
clusMat <- apply(usedMat, 1, function(x) tapply(x, st$seurat_clusters, mean)) %>% t

colFun <- colorRamp2(seq(0, 1, len = 11), brewer.pal(11, "RdYlBu") %>% rev)

pdf("SD_pSS+HC.cluster_celltype.enrich.heatmap.pdf", width = 10, height = 9)
Heatmap(
  clusMat[usedType, ] %>% set_colnames(enrichMat$type), col = colFun,
  rect_gp = gpar(color = "black"),
  cluster_rows = T, cluster_columns = T,
  cluster_column_slices = F,
  column_split = ifelse(enrichMat$lfc > 0, "pSS Enriched", "HC Enriched"),
  column_names_rot = 45,
  top_annotation = HeatmapAnnotation(
    enrich = enrichMat$lfc,
    col = list(
      enrich = colorRamp2(seq(-5, 5, len = 5), c("#2364aa","#3da5d9","#73bfb8","#fec601","#ea7317"))
    )
  ),
  bottom_annotation = HeatmapAnnotation(
    type = enrichMat$type,
    col = list(
      type = clsCol
    )
  ),
  heatmap_legend_param = list(
    title = "Celltype Prediction",
    title_position = "lefttop-rot"
  ),
  width = unit(5, "in"),
  height = unit(5, "in")
)
dev.off()

saveRDS(st, "SG_pSS+HC.st.rds")
saveRDS(st@meta.data, "SG_pSS+HC.metadata.rds")
st <- readRDS("SG_pSS+HC.st.rds")
st@meta.data <- readRDS("SG_pSS+HC.metadata.rds")
