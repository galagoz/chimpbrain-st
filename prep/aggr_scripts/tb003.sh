#$ -N spaceranger_aggr
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash

module load spaceranger/1.2.2

cd /data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/takki_st_aggr/tb003

spaceranger aggr --id=tb003 --csv=/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/aggregation_tb003.csv --normalize=mapped
