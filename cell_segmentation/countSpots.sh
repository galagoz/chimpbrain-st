#!/bin/bash
#$ -l mem_free=50G
#$ -l h_vmem=50G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-12
matlab -nodisplay -nodesktop -r "addpath(genpath('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/scripts/chimpbrain-st')); tic; try countSpots('/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/H&EbinaryImages/', '/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/data/input/', 'Takki_leFP_003_A1'); catch e, disp(e.message); end; toc; addpath /home/gokala; leave"
