#!/bin/bash
#SBATCH --partition=pi_hall
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --output=/gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/logs/run_FUSION_analysis-%J.log

# Activate conda environment
module reset
module load miniconda
conda activate /home/amm422/.conda/envs/R

# Run R script
Rscript /gpfs/gibbs/pi/ycgh/amm422/project/centromere/prediction/FUSION_analysis.R