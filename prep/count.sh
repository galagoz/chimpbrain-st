#!/bin/bash
#
# This script will run Space Ranger count
# on fastq files.
# Gokberk Alagoz - 7.7.23

##########################################
# PATHS
inDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_rawdata"
refDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/resources/ref_genomes_from10x/ref_freshf/chimp_ref"
codeDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/scripts/chimpbrain-st/prep"
imageDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/slide_images/processed"
outDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_count"

##########################################

# create a directory for grid scripts
mkdir "${codeDir}/count_scripts"

echo "Starting to run count"

for sample in ${inDir}/*; do
   
   # set current sample name
   tmp_sample=$(basename "$sample")
   echo $tmp_sample

   # get slide area
   tmp_areaID="$(cut -d'_' -f4- <<<"$tmp_sample")"

   # make output and log directories for each sample
   mkdir "${outDir}/${tmp_sample}"
   mkdir "${outDir}/${tmp_sample}/logs"
   shellFile="${codeDir}/count_scripts/${tmp_sample}.sh"
   logFile="${outDir}/${tmp_sample}/logs/${tmp_sample}.log"

   # set remaining inputs for spaceranger count
   tmp_img=$(find "$imageDir" -type f -name "*$tmp_areaID*")
   tmp_img_fname=$(basename "$tmp_img")
   slideID="$(cut -d'_' -f1 <<<"$tmp_img_fname")"

   echo '#$ -N spaceranger_count
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

module load spaceranger/1.2.2

cd '$outDir'/'${tmp_sample}'

spaceranger count --id='$tmp_sample' --transcriptome='$refDir' --fastqs='$sample' --sample='$tmp_sample' --image='$tmp_img' --slide='$slideID' --area='$tmp_areaID' --localcores=16 --localmem=128' > $shellFile
   	
   chmod a+x ${shellFile}
   echo "Created the script for cluster -> submitting ${tmp_sample} to the Grid"
   qsub -o ${logFile} -j y ${shellFile}
	   
done

echo "Done!"

##########################################
