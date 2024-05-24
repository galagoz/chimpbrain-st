# module load conda_R/3.6.x
library('SingleCellExperiment')
library('ggplot2')

## load rse list
load("/spatial_trans/generate_sce/Human_DLPFC_Visium_processedData_rseList.rda", verbose = TRUE)

sceList <- lapply(rseList, function(rse) {
    SingleCellExperiment(
        assays = list(counts = assays(rse)$umis),# save this large UMI sparse matrix
        rowData = rowData(rse), # save the gene info
        colData = colData(rse), # save the info of sample image
        metadata = metadata(rse) # save
    )
})
save(sceList, file = '/spatial_trans/generate_sce/Human_DLPFC_Visium_processedData_sceList.Rdata')
# combine all spots of different samples (only spots with transcriptome data) into one SCE
sce <- do.call(cbind, sceList) # "cbind" is to combine elements by rows
metadata(sce) <- list('image' = do.call(rbind, metadata(sce)))
table(colData(sce)$sample_name)

## Add design info (you need to design by yourself)
# in this csv file, you need to have some information describing the position of your samples, as we collected several spatial replicates within each tissue block
# 10xID	Description	Rep
# v10j27_028_a1	anw945_gfi1_rep1 - section acquired at position 0 um	replicate 01
# v10j27_028_b1	anw945_gfi1_rep2 - section acquired at position 0 um	replicate 02
# v10j27_028_c1	anw945_gfi1_rep3 - section acquired at position 300 um	replicate 01
# v10j27_028_d1	anw945_gfi1_rep4 - section acquired at position 300 um	replicate 02
study <-
    read.csv(
        '/spatial_trans/image_index_10xID.csv'
    )

## same order
stopifnot(identical(match(names(sceList), study$X10xID),
    seq_len(length(sceList))))


gsub('anw_|_rep.*', '', study$Description)
gsub('.*position | um', '', study$Description)
sce$subject <- rep(gsub('anw_|_rep.*', '', study$Description),
    table(colData(sce)$sample_name))
with(colData(sce), table(sample_name, subject))


sce$position <- rep(gsub('.*position | um', '', study$Description),
    table(colData(sce)$sample_name))
with(colData(sce), table(sample_name, position))


sce$replicate <- rep(gsub('.*0', '', study$Rep),
    table(colData(sce)$sample_name))
with(colData(sce), table(sample_name, replicate))


## For blocking later with scran
sce$subject_position <- paste0(sce$subject, '_pos', sce$position)

save(sce, file = '/spatial_trans/generate_sce/Human_DLPFC_Visium_processedData_sce.Rdata')

