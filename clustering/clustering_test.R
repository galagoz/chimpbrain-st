library(BayesSpace)
library(harmony)
library(ggplot2)
library(patchwork)
library(scater)

# use Harmony to normalize the batch effect and run graph-based clustering
# first load the sce object that you generated from the scripts that I send to you.
outDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/results/"
load(paste0(outDir, "generate_sce/Chimp_leFP_Visium_processedData_sce_scran.Rdata"))

# remove genes
count_p1=rowSums(as.matrix(counts(sce)[1:10000,]>0)+0)
count_p2=rowSums(as.matrix(counts(sce)[10001:20000,]>0)+0)
count_p3=rowSums(as.matrix(counts(sce)[20001:length(rowData(sce)$gene_name),]>0)+0)
count_spot=c(as.vector(count_p1), as.vector(count_p2), as.vector(count_p3))
keep_gene_ind=which(count_spot > round(length(sce$sample_name)*0.0001))
sce=sce[keep_gene_ind,]

## Drop mitochondrial genes
ix_mito <- grep("^MT-", rowData(sce)$gene_name)
sce <- sce[-ix_mito, ]

## lognormalize the data, get PCs

# BayesSpace is a very good toolbox, which has helped us to do the data integration, dimentionality reduction, normalzing the gene expression data, and clustering, you need to learn about this "spatialPreprocess" function, as it integrate all of these processings from Scater and Scran. To better understand this "spatialPreprocess" function, I recommend reading this manual to learn more (https://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html).
set.seed(102)
sce = spatialPreprocess(sce, n.PCs = 50, n.HVGs=2000) # top 10% genes by setting "prop = 0.1", or you can also just include top 2000 genes (suggested by the Harmony toolbox)

## Batch correction

# How does the data look prior to batch normalization?
sce = runUMAP(sce, dimred = "PCA")
colnames(reducedDim(sce, "UMAP")) = c("UMAP1", "UMAP2")

UMAP_PCA = ggplot(data.frame(reducedDim(sce, "UMAP")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce$sample_name))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()
#pdf(paste0(outDir, "clustering/UMAPpca.pdf"), height = 10, width = 10)
UMAP_PCA
#dev.off()

sce = RunHarmony(sce, "sample_name", verbose = F) # run on the PCA of sce

sce = runUMAP(sce, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(sce, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

UMAP_harmony = ggplot(data.frame(reducedDim(sce, "UMAP.HARMONY")), 
       aes(x = UMAP1, y = UMAP2, color = factor(sce$sample_name))) +
  geom_point() +
  labs(color = "Sample") +
  theme_bw()
#pdf(paste0(outDir, "clustering/UMAPharmony.pdf"), height = 10, width = 10)
UMAP_harmony
#dev.off()

# set the spot position
# BayesSpace use the "row" and "col" to do the clustering, we must let these figures separate
sample_names=unique(sce$sample_name)
for (i in 1:length(sample_names) - 1) {
  sce$col[sce$sample_name == sample_names[i + 1]] =
    150 * i + sce$col[sce$sample_name == sample_names[i + 1]]
}

# clustering analysis
set.seed(149)
# selecting the number of clusters
chimp_leFP <- qTune(sce, qs=seq(2, 20), burn.in=10, nrep=100)
qPlot(chimp_leFP)

#pdf(paste0(outDir, "clustering/qPlot_seq2-20_burnin10_nrep100.pdf"), height = 10, width = 10)
qPlot(chimp_leFP)
#dev.off()

sce = spatialCluster(chimp_leFP, use.dimred = "HARMONY", q = 12, nrep = 20000, burn.in = 1000) #use HARMONY and set the specific number of clusters you wanted to get

# save clustering results
#save(sce, file = paste0(outDir, "generate_sce/Chimp_leFP_Visium_processedData_sce_harmony_BayesSpace.Rdata"))
load(paste0(outDir, "generate_sce/Chimp_leFP_Visium_processedData_sce_harmony_BayesSpace.Rdata"))

sce = runUMAP(sce, dimred = "HARMONY", name = "UMAP.HARMONY")
colnames(reducedDim(sce, "UMAP.HARMONY")) = c("UMAP1", "UMAP2")

UMAP_harmony_clusters = ggplot(data.frame(reducedDim(sce, "UMAP.HARMONY")),
                               aes(x = UMAP1, y = UMAP2, color = factor(sce$spatial.cluster))) +
  geom_point() +
  labs(color = "Clusters") +
  theme_bw()

#pdf(paste0(outDir, "clustering/UMAPspatialCluster.pdf"), height = 10, width = 10)
UMAP_harmony_clusters
#dev.off()

# plot clusters, so you can plot these clusters based on the index from BayesSpace.
#pdf(paste0(outDir, "clustering/clusterPlot_q12_burnin1000_nrep20000.pdf"), height = 10, width = 40)
clusterPlot(sce)
#dev.off()

#pdf(paste0(outDir, "clustering/blackBorders_clusterPlot_q12_burnin1000_nrep20000.pdf"), height = 10, width = 40)
clusterPlot(sce, color="black") +
  theme_bw() +
  xlab("Column") +
  ylab("Row") +
  labs(fill="BayesSpace\ncluster", title="Spatial clustering of ST_chimp_leftFP")
#dev.off()

pdf(paste0(outDir, "clustering/featurePlot_TBR1.pdf"), height = 10, width = 40)
featurePlot(sce, "TBR1")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_MBP.pdf"), height = 10, width = 40)
featurePlot(sce, "MBP")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_RORB.pdf"), height = 10, width = 40)
featurePlot(sce, "RORB")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_MGP.pdf"), height = 10, width = 40)
featurePlot(sce, "MGP")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_CCDC68.pdf"), height = 10, width = 40)
featurePlot(sce, "CCDC68")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_SETBP1.pdf"), height = 10, width = 40)
featurePlot(sce, "SETBP1")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_FOXP1.pdf"), height = 10, width = 40)
featurePlot(sce, "FOXP1")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_FOXP3.pdf"), height = 10, width = 40)
featurePlot(sce, "FOXP3")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_FOXP4.pdf"), height = 10, width = 40)
featurePlot(sce, "FOXP4")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_CHD3.pdf"), height = 10, width = 40)
featurePlot(sce, "CHD3")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_ZIC4.pdf"), height = 10, width = 40)
featurePlot(sce, "ZIC4")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_ZIC1.pdf"), height = 10, width = 40)
featurePlot(sce, "ZIC1")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_KRT17.pdf"), height = 10, width = 40)
featurePlot(sce, "KRT17")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_AQP4.pdf"), height = 10, width = 40)
featurePlot(sce, "AQP4")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_HPCAL1.pdf"), height = 10, width = 40)
featurePlot(sce, "HPCAL1")
dev.off()

pdf(paste0(outDir, "clustering/featurePlot_TRABD2A.pdf"), height = 10, width = 40)
featurePlot(sce, "TRABD2A")
dev.off()

# spatialEnhance did not work on lux14, because it's too memory intensive
# find another way to run this
sce.enhanced <- spatialEnhance(sce, q=12, platform="ST", gamma=2,
                               nrep=100000, save.chain=TRUE)
