# This script will convert Seurat objects to 
# AnnData object, for subsequent analyses in 
# scanpy

# Gokberk Alagoz - 13.03.2024
library(Seurat)
library(SeuratObject)
library(anndata)
library(sceasy)
library(data.table)
library(SeuratDisk)
library(reticulate)
use_python("/usr/shared/apps/anaconda/3.2021.05/bin/python")

# Read Ma et al. (2022) RDS/Seurat object
setwd("/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/ma_etal2022_chimp_dlPFC_snRNA/processed")
data = readRDS("PFC_snRNAseq_liftover.rds")
data = UpdateSeuratObject(data)

# Convert Ma et al. (2022) to AnnData and save
SaveH5Seurat(data,
             filename = "PFC_snRNAseq_liftover.h5Seurat", 
             overwrite = T,
             assay = "counts")
Convert("PFC_snRNAseq_liftover.h5Seurat", dest = "h5ad")

# Read the anndata object back into R, remove the seurat
# object from your environment so your R session is faster
anndata_chimp = read_h5ad("PFC_snRNAseq_liftover.h5ad")
#rm(data)

# Subset anndata to chimp data
# rm anndata from your environment
anndata_chimp = anndata[anndata$obs["species"]=="Chimpanzee",]
#rm(anndata)

# Read the chimpanzee genome annotation file
# You need this file to fetch ensembl gene IDs into
# your anndata object, to be able to match genes with your
# visium results later on in cell2location analysis.

gtf_chimp = fread("/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/resources/ref_genomes_from10x/ref_freshf/chimp_ref/genes/protein_coding_genes_wGeneNames.gtf",
                  sep = "\t", header = FALSE, data.table = FALSE)
# clean up
gtf_chimp_clean = data.frame(
  gene_ID = sub('.*gene_id "(.*?)".*', '\\1', gtf_chimp$V9),
  gene_name = sub('.*gene_name "(.*?)".*', '\\1', gtf_chimp$V9)
)

nrow(gtf_chimp_clean) # 17,617 protein coding genes in chimp gtf file

# Intersect gene_name in the gtf file with gene names
# in  you reference anndata file. For the matching ones,
# write the ensembl gene ID to anndata$gene_ID

anndata_chimp$var$gene_ID = gtf_chimp_clean$gene_ID[match(anndata_chimp$var$features, gtf_chimp_clean$gene_name)]

#########################################################################

# Read Caglayan et al. (2023) excitatory neurons RDS object
setwd("/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/caglayan_etal2023_chimp_PCC_snRNA/")
data = readRDS("GSE192773_Seurat_Excitatory_RNA.RDS")
data = UpdateSeuratObject(data)

# Convert SYMBOLs to ENSEMBL IDs
data$SYMBOL = rownames(data)
require(org.Pt.eg.db)
ensIDs = mapIds(org.Pt.eg.db,
                keys = data$SYMBOL,
                column = 'ENSEMBL',
                keytype = 'SYMBOL')
all(data$SYMBOL == names(ensIDs))

keep <- !is.na(ensIDs)
ensIDs <- ensIDs[keep]
data <- data[keep,]
data$ensID = unname(ensIDs[keep])

# Convert RDS to AnnData and save
sceasy::convertFormat(data, from="seurat", to="anndata",
                      outFile='PCC_excitatory_snRNAseq.h5ad')

#########################################################################

# Read Khrameeva et al. (2020) ACC RDS object
setwd("/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/khrameeva_etal2020_chimp_ACC_snRNAseq/")
data = readRDS("GSE127774_ACC_seurat.rds")
data = UpdateSeuratObject(data)

names(data@meta.data) = c("orig_ident",
                          "nCount_RNA",
                          "nFeature_RNA")

# Convert RDS to AnnData and save
sceasy::convertFormat(data, from="seurat", to="anndata",
                      outFile='PCC_excitatory_snRNAseq.h5ad')