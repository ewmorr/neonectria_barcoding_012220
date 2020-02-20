#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="RDP_class"
#SBATCH --output=dada2.out
#SBATCH --partition=shared
##SBATCH --cpus-per-task=8
##SBATCH --mem=128000

module purge
module load anaconda/colsa

###
#Rerun taxonomic classification in dada2

conda activate dada2-check

# first run dada2 for the concatenated files
cd ~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/sep_run_pool #this was originally run with last year's UNITE db

srun Rscript ~/repo/neonectria_barcoding_012220/UNITE_taxonomic_classification-slurm.r
