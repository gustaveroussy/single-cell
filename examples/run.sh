#!/bin/bash

########################################################################
## Single-cell RNA-seq script to launch single-cell RNA-seq pipeline
##
## using: sbatch /mnt/beegfs/userdata/m_aglave/pipeline/run.sh
##
########################################################################

## JOB PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#SBATCH --job-name=pipeline_sc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=mediumq

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/m_aglave/.environnement_conda/scRNAseq_10X_user
module load singularity

#print environment tools versions
python --version
snakemake --version
singularity --version

#parameters
path_to_configfile="/mnt/beegfs/userdata/m_aglave/pipeline/scRNAseq_10X/Param_snakfile_alignment.yaml"
path_to_pipeline="/mnt/beegfs/pipelines/single-cell"

#launch
snakemake --profile ${path_to_pipeline}/profiles/slurm -s ${path_to_pipeline}/Snakefile --configfile ${path_to_configfile}


conda deactivate

