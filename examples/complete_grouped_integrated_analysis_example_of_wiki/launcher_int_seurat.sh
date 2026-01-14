#!/bin/bash

########################################################################
## Single-cell RNA-seq script to launch single-cell RNA-seq pipeline
##
## using: sbatch /mnt/beegfs02/userdata/m_aglave/pipeline/bigr_single-cell/2.0.0/examples/complete_grouped_integrated_analysis_example_of_wiki/launcher_int_seurat.sh
##
########################################################################

## JOB PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#SBATCH --job-name=dev_seurat
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=longq

source /mnt/beegfs02/software/recherche/miniconda/25.1.1/etc/profile.d/conda.sh
conda activate /mnt/beegfs02/pipelines/bigr_single-cell/2.0.0/envs/compiled_conda/snakemake
module load singularity/3.10.5

python --version
snakemake --version
singularity --version

#parameters
path_to_configfile="/mnt/beegfs02/pipelines/bigr_single-cell/2.0.0/examples/complete_grouped_integrated_analysis_example_of_wiki/Params_int_seurat.yaml"
path_to_pipeline="/mnt/beegfs02/pipelines/bigr_single-cell/2.0.0"

#launch
snakemake --profile ${path_to_pipeline}/profiles/slurm -s ${path_to_pipeline}/Snakefile --configfile ${path_to_configfile}

conda deactivate
