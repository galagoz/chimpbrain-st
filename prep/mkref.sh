#!/bin/bash
#
# This script will make a reference genome
# assembly for Space Ranger.
# Gokberk Alagoz - 7.7.23

module load spaceranger/1.2.2

##########################################
# PATHS
rawgenomeDir="/data/workspaces/lag/workspaces/lg-spatial-transcriptomics/working_data/chimp_brain/gokberk/resources/GCF_028858775.1/ncbi_dataset/data/GCF_028858775.1"

spaceranger mkref --genome=mPanTro3-v1.1 --fasta=${rawgenomeDir}/GCF_028858775.1_NHGRI_mPanTro3-v1.1-hic.freeze_pri_genomic.fna --genes=${rawgenomeDir}/genomic.gtf
