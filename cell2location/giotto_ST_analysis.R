# This script will ...

# Gokberk Alagoz - 15.03.2024
# written following the tutorial:
# https://rubd.github.io/Giotto_site/articles/mouse_visium_brain_201226.html

library(Giotto)
library(Seurat)
library(reticulate)
library(anndata)
library(scRNAseq)
library(zellkonverter)

# 1. set working directory
#results_folder = '/path/to/directory/'
results_folder = "/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/snRNA_integration/giotto_results/"
setwd(results_folder)

# 2. set giotto python path
# set python path to your preferred python version path
# set python path to NULL if you want to automatically install (only the 1st time) and use the giotto miniconda environment
python_path = NULL 
if(is.null(python_path)) {
  installGiottoEnvironment()
}

#-Part1-#######################################################

# Make a cell-type=specific signatures matrix using the external
# snRNAseq data

chimpData = subset(data, species == 'Chimpanzee')

norm <- GetAssayData(object = chimpData, slot = "data", assay = "RNA") # this gives you a sparse matrix - run the next line to convert it to a matrix
norm <- as.matrix(norm) 

signMatrix = makeSignMatrixRank(norm,
                                chimpData$subclass,
                                ties_method = c("random", "max"),
                                gobject = NULL)

#-Part2-#######################################################

raw_counts = assay(sce, "counts")
raw_counts = assay(sce, "counts")[which(!duplicated(row.names(assay(sce, "counts")))),]
raw_counts = raw_counts[,which(!duplicated(colnames(assay(sce, "counts"))))]
cell_locs = data.frame(sce[row.names(raw_counts)]$imagerow, sce[row.names(raw_counts)]$imagecol)

# 2. use an existing matrix and data.table
my_giotto_object = createGiottoObject(raw_exprs = raw_counts,
                                      spatial_locs = cell_locs)

#-Part3-#######################################################

test_object = createGiottoVisiumObject(
  visium_dir = "/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_count/Takki_leFP_003_A1/Takki_leFP_003_A1/outs",
  h5_visium_path = "/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_count/Takki_leFP_003_A1/Takki_leFP_003_A1/outs/filtered_feature_bc_matrix.h5",
  #h5_gene_ids = c("symbols", "ensembl"),
  h5_tissue_positions_path = "/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_count/Takki_leFP_003_A1/Takki_leFP_003_A1/outs/spatial/tissue_positions_list.csv",
  h5_image_png_path = "/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_count/Takki_leFP_003_A1/Takki_leFP_003_A1/outs/spatial/tissue_lowres_image.png"
  )

spatPlot(gobject = test_object, cell_color = 'in_tissue', show_image = T, point_alpha = 0.7,
         save_param = list(save_name = 'spatplot_image'))

spatPlot(gobject = test_object, cell_color = 'in_tissue', point_size = 2,
         cell_color_code = c('0' = 'lightgrey', '1' = 'blue'),
         save_param = list(save_name = '2_c_in_tissue'))

metadata = pDataDT(test_object)
in_tissue_barcodes = metadata[in_tissue == 1]$cell_ID
visium_brain = subsetGiotto(test_object, cell_ids = in_tissue_barcodes)

## filter
visium_brain <- filterGiotto(gobject = visium_brain,
                             expression_threshold = 1,
                             gene_det_in_min_cells = 50,
                             min_det_genes_per_cell = 1000,
                             expression_values = c('raw'),
                             verbose = T)

## normalize
visium_brain <- normalizeGiotto(gobject = visium_brain, scalefactor = 6000, verbose = T)

## add gene & cell statistics
visium_brain <- addStatistics(gobject = visium_brain)

## visualize
spatPlot2D(gobject = visium_brain, show_image = T, point_alpha = 0.7,
           save_param = list(save_name = '2_d_spatial_locations'))

spatPlot2D(gobject = visium_brain, show_image = T, point_alpha = 0.7,
           cell_color = 'nr_genes', color_as_factor = F,
           save_param = list(save_name = '2_e_nr_genes'))

## highly variable genes (HVG)
visium_brain <- calculateHVG(gobject = visium_brain,
                             save_param = list(save_name = '3_a_HVGplot'))

## run PCA on expression values (default)
gene_metadata = fDataDT(visium_brain)
featgenes = gene_metadata[hvg == 'yes' & perc_cells > 3 & mean_expr_det > 0.4]$gene_ID

visium_brain <- runPCA(gobject = visium_brain, 
                       genes_to_use = featgenes, 
                       scale_unit = F, center = T, 
                       method="factominer")

screePlot(visium_brain, ncp = 30, save_param = list(save_name = '3_b_screeplot'))

rank_enrich = runRankEnrich(
  visium_brain,
  signMatrix,
  expression_values = c("raw"),
  reverse_log_scale = F,
  output_enrichment = c("original"),
  p_value = T
)

# 1.5 visualizations
cell_types_subset = colnames(signMatrix)[1:10]
spatCellPlot(gobject = visium_brain, 
             spat_enr_names = 'RANK',
             cell_annotation_values = cell_types_subset,
             cow_n_col = 4,coord_fix_ratio = NULL, point_size = 0.75,
             save_param = list(save_name="7_b_spatcellplot_1"))

####################################################

adata <- zellkonverter::SCE2AnnData(sce)

write_h5ad(adata, paste0(results_folder, "Takki_leFP_ST.h5ad"))

ST_gobject <- anndataToGiotto(adata)