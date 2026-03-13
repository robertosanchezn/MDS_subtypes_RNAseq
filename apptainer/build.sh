#!/bin/bash
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --output=build.log

# If Apptainer_builder.sif doesn't exist, build it

if [ ! -f "Apptainer_builder.sif" ]; then
    echo "Building builder container"
    apptainer build Apptainer_builder.sif Apptainer_builder.def
fi

# Export locked environments
echo "Exporting locked environment"

apptainer exec Apptainer_builder.sif /opt/conda/bin/conda env export -n base > environment.yml

echo "Building definitive container"
# Build the Apptainer container
apptainer build --force Apptainer.sif Apptainer.def


