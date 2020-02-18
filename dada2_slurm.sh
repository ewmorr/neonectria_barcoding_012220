#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name="r-dada2"
#SBATCH --output=dada2.out
#SBATCH --partition=shared

module purge
module load anaconda/colsa

###
#This script runs the primary dada2 algorithm, i.e., assuming that quality filtering and ITS extraction has already been performed

conda activate r-dada2_env

# first run dada2 for the concatenated files
cd ~/GARNAS_neonectria_barcoding_runOneAndTwo_020320/run1_run2_files_cat
mkdir dada2_processing_tables_figs
mkdir intermediate_RDS

Rscript ~/repo/neonectria_barcoding_012220/dada2-slurm.r
Rscript ~/repo/neonectria_barcoding_012220/UNITE_taxonomic_classification-slurm.r

mkdir dada2_out
Rscript ~/repo/neonectria_barcoding_012220/dada2_tables_to_file.r
