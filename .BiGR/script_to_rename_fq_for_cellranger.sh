#!/bin/bash
########################################################################
## Script to rename fastq file for CellRanger.
##
## using:
## FQ_PATH="/mnt/beegfs/scratch/bioinfo_core/B21061_NADR_11/data_input/"
## OUTPUT_PATH="/mnt/beegfs/scratch/bioinfo_core/B21061_NADR_11/CellRanger/0_rename_fastq/"
## bash script_to_rename_fq_for_cellranger.sh
##
## examples fo results:
## srun rsync -cvr HIP_1_S13_I1_001.fastq.gz > HIP_S13_L001_I1_001.fastq.gz
## srun rsync -cvr HIP_2_S13_I1_001.fastq.gz > HIP_S13_L002_I1_001.fastq.gz
## srun rsync -cvr HIP_3_S13_I1_001.fastq.gz > HIP_S13_L003_I1_001.fastq.gz
## srun rsync -cvr HIP_4_S13_I1_001.fastq.gz > HIP_S13_L004_I1_001.fastq.gz
## srun rsync -cvr HIP_1_S13_R1_001.fastq.gz > HIP_S13_L001_R1_001.fastq.gz
########################################################################

#FQ_PATH="/mnt/beegfs/scratch/bioinfo_core/B21061_NADR_11/data_input/"
#OUTPUT_PATH="/mnt/beegfs/scratch/bioinfo_core/B21061_NADR_11/CellRanger/0_rename_fastq/"

mkdir -p ${OUTPUT_PATH}

for FILE in $(find ${FQ_PATH}'/' -name "*.f*q.gz" ! -name "*L00*.f*q.gz" | sort) ; do
FILE=$(basename $FILE)
echo ${FILE};
[[ ${FILE} =~ _[0-4]_ ]];
BEGIN_FILE=${FILE%${BASH_REMATCH}*};
[[ ${BASH_REMATCH} =~ [0-4]_ ]];
NB=${BASH_REMATCH};
[[ ${FILE} =~ _S[0-9]+_ ]];
S=${BASH_REMATCH};
END_FILE=${FILE#*${BASH_REMATCH}};
NEW_NAME=${BEGIN_FILE}${S}'L00'${NB}${END_FILE};
echo ${NEW_NAME};
rsync -cvr ${FQ_PATH}'/'${FILE} ${OUTPUT_PATH}'/'${NEW_NAME};
done
