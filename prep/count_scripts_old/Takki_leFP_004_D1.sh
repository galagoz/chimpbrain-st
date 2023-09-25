#$ -N spaceranger_count
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

module load spaceranger/1.2.2

cd /data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_count/Takki_leFP_004_D1

spaceranger count --id=Takki_leFP_004_D1 --transcriptome=/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/resources/ref_genomes_from10x/ref_freshf --fastqs=/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_rawdata/Takki_leFP_004_D1 --sample=Takki_leFP_004_D1 --image=/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/slide_images/processed/V13J17-280_Takki_leFP_004_D1_processed.tif --slide=V13J17-280 --area=D1 --localcores=16 --localmem=128
