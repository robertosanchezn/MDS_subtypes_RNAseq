#!/bin/bash
#SBATCH -p short
#SBATCH -o apptainer/build.stdout.log
#SBATCH --mem 64GB

apptainer build apptainer/Apptainer.sif apptainer/Apptainer.def 
