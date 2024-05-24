# screen -S scran
# qrsh -l mem_free=60G,h_vmem=60G,h_fsize=100G -pe local 4
# module load conda_R/3.6.x
library('SingleCellExperiment')
library('scran')
library('scater')
library('BiocParallel')
library('PCAtools')
library('igraph')
library('ggplot2')
library('cowplot')
library('jaffelab') ## for ss(), splitit(), myplclust(); actually they are from rafalib
library('pheatmap')
library("spatialLIBD") # this toolbox also includes functions, we must load it

# library('zinbwave')
# library('clusterExperiment')
#
# library('RColorBrewer')

dir.create('/spatial_trans/pdf_scran', showWarnings = FALSE)
dir.create('/spatial_trans/rda_scran', showWarnings = FALSE)

## From convert_sce.R
load('/spatial_trans/geom_spatial.Rdata', verbose = TRUE)

## For resuming and other analyses
load('/spatial_trans/generate_sce/Human_DLPFC_Visium_processedData_sce.Rdata', verbose = TRUE)
## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#2_setting_up_the_data
qcstats <- perCellQCMetrics(sce) # doing QC for each spot ('perCellQCMetrics' is the function of scuttle package)
qcfilter <- quickPerCellQC(qcstats)
colSums(as.matrix(qcfilter))

with(qcfilter, table(low_lib_size, low_n_features, discard))

table(sce$sample_name[qcfilter$discard])

## Plot discarded umis
source('/spatial_trans/test/spatialLIBD_global_plot_code.R')
rm(sce_image_clus_p)
sce$discard <- qcfilter$discard
plots_discard <-
  lapply(unique(sce$sample_name), function(sampleid) {
    sce_image_clus(sce, sampleid, 'discard', colors = c('light blue', 'red'))
  })
pdf('/spatial_trans/pdf_scran/discarded_cells_grid.pdf',
    height = 50,
    width = 55) # if there are more figures, we need to scale the height and width, otherwise the spots do not fit
plot_grid(plotlist = plots_discard)
dev.off()

## We have decided not to filter umis since they seem to be layer specific, but if we want to exclude, please use the following code
sce <- sce[,!qcfilter$discard]

## if you want to save spots in a given sample, please use the following code
# aa=qcfilter$discard==T & sce@colData$sample_name==c("v10j29_114_b1")
# qcfilter$discard[aa==T]=F
# sce$discard <- qcfilter$discard

## remove the folded regions
sampleid=unique(sce$sample_name)
for (i in 1:length(sampleid)) {
  sce_sub=sce[,sce$sample_name==sampleid[i]]
  sce_sub_info=data.frame(barcode=sce_sub$barcode, imagecol=sce_sub$imagecol, imagerow=sce_sub$imagerow)
  write.csv(sce_sub_info,file=paste('/spatial_trans/folded_region/data/',sampleid[i],'.csv',sep=""),row.names=F,quote=F)
}

# manually label the folded region
# add the labels back to sce 
fold=c()
for (i in 1:length(sampleid)) {
  fold_label=read.csv(paste('/spatial_trans/folded_region/data/',sampleid[i],'_fold.csv',sep=""))
  fold=c(fold,as.logical(fold_label$fold))
}
sce$fold <- fold
sce <- sce[,!sce$fold]

# take a look at the gene names 
aa=data.frame(gene_id=rowData(sce)$gene_id,gene_name=rowData(sce)$gene_name)
write.table(aa,file="/spatial_trans/geneid_genename.txt",quote=F,sep="\t",row.names=F,col.names=F)


# plots umis after rmoving these spots
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
plots_umis=list()
for (i in 1:length(sampleid)){
  sce_sub=sce[, sce$sample_name == sampleid[i]]
  d = as.data.frame(colData(sce_sub))
  plots_umis[[i]]=ggplot(d, aes(x=imagecol,y=imagerow,fill=sum_umi)) + # draw the frame
    geom_spatial(data=metadata(sce_sub)$image[i,],
                 aes(grob=grob), x=0.5, y=0.5) + # add the background image
    geom_point(shape = 21,  size = 1.25, stroke = 0.25) + # draw the circles
    coord_cartesian(expand=FALSE)+ # expand to the boundary of the image
    scale_fill_gradientn(colours = myPalette(100), limits=c(0,30000))+ # adjust the color bar
    xlim(0,max(sce_sub$width)) + # scale the x axis to fit the background
    ylim(max(sce_sub$height),0) + # scale the y axis to fit the background
    xlab("") + ylab("") + # remove the label of x and y axis
    labs(fill = "Total UMI")+
    ggtitle(unique(sce_sub$sample_name)) +
    theme_set(theme_bw(base_size = 10))+
    theme(panel.grid.major = element_blank(), # set the background as white
          panel.grid.minor = element_blank(),
          panel.background = element_blank(), 
          axis.line = element_line(colour = "black"), # set the outline as black
          axis.text = element_blank(),
          axis.ticks = element_blank())
}
pdf('/data/clusterfs/lag/users/zhisha/spatial_trans/pdf_scran/umi_by_sample_exc.pdf',
    height = 50,
    width = 55) # if there are more figures, we need to scale the height and width, otherwise the spots do not fit
plot_grid(plotlist = plots_umis)
dev.off()


## Read in the number of cells per spot (if you do not need this, you can also remove this section, it does not affect follow-up analyses)
cells <- do.call(rbind, lapply(dir('/spatial_trans/histology/cell_spot'), function(sampleid) {
  x <- read.csv(file.path('/spatial_trans/histology/cell_spot', sampleid, 'tissue_spot_counts2.csv'))
  x$key <- paste0(sampleid, '_', x$barcode)
  return(x[, c('key', 'count')])
}))


## Used in plotly code in spatialLIBD
sce$key <- paste0(sce$sample_name, '_', colnames(sce))
m <- match(sce$key, cells$key)
stopifnot(!all(is.na(m)))
sce$cell_count <- cells$count[m]

tapply(sce$cell_count, sce$sample_name, summary)


# save this sce, which is the most important part of all data.
# if you have done this part, congratulations! Processing is done. Next step is to run the analyses that you are interested in based on "sca" structure.
save(sce, file = '/spatial_trans/generate_sce/Human_DLPFC_Visium_processedData_sce_scran.Rdata')

