## module load conda_R/3.6.x # devel

## ----Libraries ------------------
BiocManager::install("BiocNeighbors",version="3.10")
library(tidyverse)
library(ggplot2)
library(Matrix)
library(Rmisc)
library(ggforce)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(hdf5r) # must "module load hdf5/1.10.5" to use the package
library(readbitmap)
library(Seurat)
library(tibble)
library(SummarizedExperiment)
library(rtracklayer)  # must "module unload anaconda" to install the package

## Function for plotting
geom_spatial <-  function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = FALSE,
                         ...) {
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}
save(geom_spatial,file='/spatial_trans/geom_spatial.Rdata')

## get sample names
sample_names <- read.delim("/spatial_trans/lenas.txt", as.is=TRUE, header=FALSE)$V1
sample_names
# in this lenas.txt, you need to list your sample names, just one column, such as 
# v10j27_028_a1
# v10j27_028_b1
# v10j27_028_c1
# v10j27_028_d1

##  10x output path
path = "/spatial_trans/data/"

## output
# before load your samples, you need to copy your data to a folder, with each sample in each subfold. In the each subfolder, it contains the following files.
# metrics_summary_csv.csv
# scalefactors_json.json
# tissue_hires_image.png
# tissue_lowres_image.png
# tissue_positions_list.txt
# v10j27_028_a1_analysis__clustering_graphclust_clusters.csv ("v10j27_028_a1" is just the name of your sample)
# v10j27_028_a1_filtered_feature_bc_matrix.h5
# v10j27_028_a1_filtered_feature_bc_matrix__barcodes.tsv.gz
# v10j27_028_a1_filtered_feature_bc_matrix__features.tsv.gz
# v10j27_028_a1_filtered_feature_bc_matrix__matrix.mtx.gz
# v10j27_028_a1_raw_feature_bc_matrix__features.tsv.gz
# v10j27_028_a1_web_summary.html
# all of these files can be obtained from preprocessing results for each sample, Meggie did the preprocessing analysis.
image_paths <- paste0(path, sample_names, "/tissue_lowres_image.png")
scalefactor_paths <- paste0(path, sample_names, "/scalefactors_json.json")
tissue_paths <- paste0(path, sample_names, "/tissue_positions_list.txt")
cluster_paths <- paste0(path, sample_names, "/", sample_names, "_analysis__clustering_graphclust_clusters.csv")
matrix_paths <- paste0(path, sample_names, "/", sample_names, "_filtered_feature_bc_matrix.h5")

all(file.exists(c(image_paths, scalefactor_paths, tissue_paths, cluster_paths, matrix_paths)))
# TRUE

## get annotation (you need to revise this part based on the naming rules of your samples)
map = read.delim("/spatial_trans/data/v19s18_052_a1/v19s18_052_a1_raw_feature_bc_matrix__features.tsv.gz",
	as.is=TRUE, header=FALSE,col.names=c("EnsemblID", "Symbol", "Type")) # just use this v19s18_052_a1 as an example, you can also load one of your samples, and different samples have the same set of gene_id and gene_name
## get GTF, this seems like what they used
gtf = import("/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/primary_data/reference_genome/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
gtf = gtf[gtf$type == "gene"] # this gtf is S4 class in R, we can use [] to extract the data, we only extract "genes", we don't want "transcript", "exon" or others 
# 'map' contains our genes, 'gtf' contains reference genes. I find map$EnsemblID does not match gtf$gene_id, so I need to find overlapping genes
names(gtf) = gtf$gene_id # extract gtf$gene_id information to give it to @NAMES of gtf, so each row now has a name
# I find the same gene_id corresponds different gene_name in 'gtf' and 'map', and some gene_name has different gene_id, but I finally only use gene_id
inter_gene_id=intersect(map$EnsemblID,gtf$gene_id)
inter_gene_symbol_map=map$Symbol[ match(inter_gene_id,map$EnsemblID)]
inter_gene_symbol_gtf=gtf$gene_name[ match(inter_gene_id,gtf$gene_id)]
ind=which(inter_gene_symbol_map==inter_gene_symbol_gtf)
inter_gene_id_final=inter_gene_id[ind]
gtf = gtf[inter_gene_id_final] # only want overlapping genes between ref and our dataset; NOTE: intersect(map$Symbol,gtf$gene_name) is not intersect(map$EnsemblID,gtf$gene_id)
seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25]) # add "chr" in front of chormosome to the seqnames of gtf
mcols(gtf) = mcols(gtf)[,c(5:9)] # extract the 5th to 9th columns

## ------------------------------------------------------------------------
images_cl <- lapply(image_paths, read.bitmap) # we can see there are three dimensions, including length, width and channel (RGB)
dims = t(sapply(images_cl, dim)) # calculate the dimensions of the images
colnames(dims) = c("height", "width", "channel") # replace the colnames
dims = as.data.frame(dims) # convert the matrix to data.frame


## ------------------------------------------------------------------------ 
# gather infomation from low-resolusion images of capture areas and save them in the "rseList"
grobs <- lapply(images_cl, rasterGrob, width=unit(1,"npc"), height=unit(1,"npc")) # Render a bitmap image at the given location, size, and orientation
images_tibble <- tibble(sample=sample_names, grob=grobs) # construct a new data.frame, the same as base::data.frame()
images_tibble$height = dims$height
images_tibble$width = dims$width
images_tibble


## ------------------------------------------------------------------------
scales <- lapply(scalefactor_paths, function(x) fromJSON(file=x))


## ------------------------------------------------------------------------
# gather the clustering results (only have spots with transcriptome data) and save it in the "rseList"
clusters = lapply(cluster_paths, read.csv)
head(clusters[[1]])


## ------------------------------------------------------------------------
bcs <- list()
for (i in 1:length(sample_names)) {
   bcs[[i]] <- read.csv(tissue_paths[i],col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
   bcs[[i]]$sample_name <- sample_names[i]
   bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_lowres_scalef # scale tissue coordinates for lowres image
   bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_lowres_scalef # scale tissue coordinates for lowres image
   bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue) # convert it to factor ("0" means no-tissue, "1" means yes-tissue)
   bcs[[i]] <- merge(bcs[[i]], clusters[[i]], by.x = "barcode", by.y = "Barcode", all = TRUE) # merge data
   bcs[[i]]$height <- images_tibble$height[i]
   bcs[[i]]$width <- images_tibble$width[i]
}
names(bcs) <- sample_names

head(bcs[[1]])


## ------------------------------------------------------------------------

## keep in regular genomics formatting
# we need to know that the row name of umiList is gene_name, rather than gene_id
umiList <- lapply(matrix_paths, Read10X_h5) # collect UMI information for each sample (capture area); umiList contains large sparse matrix
names(umiList) = sample_names
sapply(umiList, dim)
# because we only include a subset of genes that overlap with genes in the reference, so we need to extract
for (i in 1:length(bcs)) {
   umiList[[i]]=umiList[[i]][gtf$gene_name,]
}
# add the image information and clustering results into the "rseList" and save it
# the full name of rse is RangedSummarizedExperiment
rseList = mapply(function(u, bc) {
	rownames(bc) = bc$barcode
	bc = bc[colnames(u),c(1,7,2:6,8:ncol(bc))]
	rse = SummarizedExperiment(
		assays = list('umis' = u),
		rowRanges = gtf, colData = bc) # gtf and umiList should have the same subset of genes
	rse$sum_umi = colSums(u)
	rse$sum_gene = colSums(u > 0)
	return(rse)
}, umiList, bcs)

## add images
for(i in seq(along=rseList)) {
	metadata(rseList[[i]])$image = images_tibble[i,]
}
## save out
save(rseList, geom_spatial,
	file = "/spatial_trans/generate_sce/Human_DLPFC_Visium_processedData_rseList.rda")

