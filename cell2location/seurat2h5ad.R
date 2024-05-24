library(Seurat)
library(tidyverse)
library(Matrix)
library(biomaRt)
library(org.Pt.eg.db)
library(SeuratDisk)

setwd("/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/ma_etal2022_chimp_dlPFC_snRNA/processed/PFC-snRNA-seq_Chimpanzee/")
exp_matrix = ReadMtx(mtx = "snRNA-seq_Chimpanzee_counts.mtx.gz",
                     cells = "snRNA-seq_Chimpanzee_cell_meta.txt.gz",
                     features = "snRNA-seq_Chimpanzee_genes.txt.gz", feature.column = 1, skip.cell = 1)

seurat_object = CreateSeuratObject(counts = exp_matrix)
seurat_object_S5 = UpdateSeuratObject(seurat_object)
seurat_object_S5[["RNA"]] <- as(object = seurat_object_S5[["RNA"]], Class = "Assay")

# Convert SYMBOLs to ENSEMBL IDs
seurat_object_S5$SYMBOL = rownames(seurat_object_S5)
require(org.Pt.eg.db)
ensIDs = mapIds(org.Pt.eg.db,
                keys = seurat_object_S5$SYMBOL,
                column = 'ENSEMBL',
                keytype = 'SYMBOL')
all(unname(seurat_object_S5$SYMBOL) == names(ensIDs))

keep <- !is.na(ensIDs)
ensIDs <- ensIDs[keep]
seurat_object_S5 <- seurat_object_S5[keep,]
seurat_object_S5$ensID = unname(ensIDs[keep])

SaveH5Seurat(seurat_object_S5, filename = "snRNA-seq_Chimpanzee.h5Seurat", overwrite = T)
Convert("snRNA-seq_Chimpanzee.h5Seurat", dest = "h5ad")
