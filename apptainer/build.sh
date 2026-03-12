#!/bin/bash
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --output=build.log

apptainer build --force Apptainer.sif Apptainer.def