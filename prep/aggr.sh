#!/bin/bash
#
# This script will run the aggr function of space ranger
# to aggregate data from consecutive tissue sections into 
# a single feature-barcode matrix.
#
# To run this, you need to make the "aggregation csv".
# This csv should contain library_id, molecule_h5, cloupe_file, spatial_folder columns.
#
# IMPORTANT NOTE: We did not use this script (i.e. aggregate data from consecutive sections).
#
# For details, see https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/using/aggregate
# Gokberk Alagoz - 10.07.23

##########################################
# PATHS
inDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data"
codeDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/scripts/chimpbrain-st/prep"
outDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_aggr"

##########################################

# create a directory for grid scripts
mkdir "${codeDir}/aggr_scripts"

echo "Starting to run aggr"

for block in ${inDir}/aggregation_*.csv; do

   # set current tissue block name
   tmp_block=$(basename "$block") # Get rid of the path name
   tmp_block="${tmp_block#*_}" # Use Parameter Expansion (PE) to strip off the part before '-'
   tmp_block="${tmp_block%%.*}" # Use PE again to strip after the first '.'
   echo $tmp_block

   # make output and log directories for each block
   mkdir "${outDir}/${tmp_block}"
   mkdir "${outDir}/${tmp_block}/logs"
   shellFile="${codeDir}/aggr_scripts/${tmp_block}.sh"
   logFile="${outDir}/${tmp_block}/logs/${tmp_block}.log"

   echo '#$ -N spaceranger_aggr
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

module load spaceranger/1.2.2

cd '$outDir'/'${tmp_block}'

spaceranger aggr --id='$tmp_block' --csv='$block' --normalize=mapped' > $shellFile

   chmod a+x ${shellFile}
   echo "Created the script for cluster -> submitting ${tmp_block} to the Grid"
   qsub -o ${logFile} -j y ${shellFile}

done

echo "Done!"

##########################################
