#!/bin/bash

#SBATCH --partition=short
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --time=0-01:00:00
#SBATCH --output=%j_%x.log.out
#SBATCH --error=%j_%x.log.err

module load R-base/4.3.0
module load R-cbrg/current

# Check script for number of cores needed
Rscript --vanilla update_bioc.R
