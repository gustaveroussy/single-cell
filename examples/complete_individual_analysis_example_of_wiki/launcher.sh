#!/bin/bash

########################################################################
## Single-cell script to launch single-cell pipeline
##
## using: sbatch /mnt/beegfs02/pipelines/bigr_single-cell/2.0.0/examples/complete_individual_analysis_example_of_wiki/launcher.sh
########################################################################
#SBATCH --job-name=pipeline_sc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=mediumq

source /mnt/beegfs02/software/recherche/miniconda/25.1.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs02/pipelines/bigr_single-cell/2.0.0/envs/compiled_conda/snakemake
module load singularity/3.10.5

#print environment tools versions
python --version
snakemake --version
singularity --version

#parameters
path_to_configfile="/mnt/beegfs02/pipelines/bigr_single-cell/2.0.0/examples/complete_individual_analysis_example_of_wiki/Params.yaml"
path_to_pipeline="/mnt/beegfs02/pipelines/bigr_single-cell/2.0.0"

#launch
snakemake --profile ${path_to_pipeline}/profiles/slurm -s ${path_to_pipeline}/Snakefile --configfile ${path_to_configfile}


conda deactivate
