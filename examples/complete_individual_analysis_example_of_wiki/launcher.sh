#!/bin/bash

########################################################################
## Single-cell script to launch single-cell pipeline
##
## using: sbatch /mnt/beegfs/userdata/m_aglave/pipeline/single-cell/examples/complete_individual_analysis_example_of_wiki/launcher.sh
########################################################################
#SBATCH --job-name=pipeline_sc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=mediumq

conda activate /mnt/beegfs/userdata/m_aglave/.environnement_conda/single_cell_user
module load singularity


path_to_configfile="/mnt/beegfs/userdata/m_aglave/pipeline/single-cell/examples/complete_individual_analysis_example_of_wiki/Params.yaml"
path_to_pipeline="/mnt/beegfs/pipelines/single-cell"

snakemake --profile ${path_to_pipeline}/profiles/slurm -s ${path_to_pipeline}/Snakefile --configfile ${path_to_configfile}

conda deactivate
