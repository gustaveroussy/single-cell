#!/bin/bash
########################################################################
## Script to concat fastq file from 2 lanes into one file for CellRanger.
##
## using:
## FQ_PATH="/mnt/beegfs/scratch/bioinfo_core/B21061_NADR_11/data_input/"
## OUTPUT_PATH="/mnt/beegfs/scratch/bioinfo_core/B21061_NADR_11/CellRanger/0_rename_fastq/"
## bash script_to_concat_rename_fq_for_cellranger.sh
##
## examples fo results:
## srun cat HIP_1_S13_L001_I1_001.fastq.gz HIP_1_S13_L002_I1_001.fastq.gz > HIP_S13_L001_I1_001.fastq.gz
## srun cat HIP_2_S13_L001_I1_001.fastq.gz HIP_2_S13_L002_I1_001.fastq.gz > HIP_S13_L002_I1_001.fastq.gz
## srun cat HIP_3_S13_L001_I1_001.fastq.gz HIP_3_S13_L002_I1_001.fastq.gz > HIP_S13_L003_I1_001.fastq.gz
## srun cat HIP_4_S13_L001_I1_001.fastq.gz HIP_4_S13_L002_I1_001.fastq.gz > HIP_S13_L004_I1_001.fastq.gz
########################################################################

#FQ_PATH="/mnt/beegfs/scratch/bioinfo_core/B21061_NADR_11/data_input/"
#OUTPUT_PATH="/mnt/beegfs/scratch/bioinfo_core/B21061_NADR_11/CellRanger/0_rename_fastq/"

mkdir -p ${OUTPUT_PATH}
cd ${OUTPUT_PATH}

#concat files
nb_tot_files=$(ls ${FQ_PATH}*L001*.f*q.gz | wc -l)
nb_files=0
echo "nb_tot_files: " $nb_tot_files
echo "nb_files: " $nb_files
while ((${nb_files} != ${nb_tot_files}))
do
  nb_jobs=$(squeue -u ${USER} | grep "cat"| wc -l) #check your job in squeue
  echo "number of active jobs: " $nb_jobs
  if [[ ${nb_jobs} -ge 10 ]] #if there are more than 10 jobs
  then
    echo "sleep... "
    sleep 30s
  else
    nb_files=$((nb_files + 30))
    echo "nb_files: " $nb_files
    if [[ ${nb_files} -le ${nb_tot_files} ]] #si ça dépasse pas nb_tot_files: alors on est pas sur les derniers fichiers
    then
      list_files=$(ls ${FQ_PATH}*L001*.f*q.gz | head -n${nb_files} | tail -n30)
    else
      echo "Last files..."
      nb_files=$((nb_files - 30))
      nb_remaining_files=$((${nb_tot_files} - ${nb_files}))
      echo "nb_tot_files: " $nb_tot_files
      echo "nb_files: " $nb_files
      echo "nb_remaining_files: " $nb_remaining_files
      nb_files=${nb_tot_files}
      list_files=$(ls ${FQ_PATH}*L001*.f*q.gz | tail -n${nb_remaining_files})
    fi
    if [ ! -z "$list_files" ] # if remaining files
    then
      echo "list_files: " $list_files
      for L001_FILE in ${list_files}
      do
        L001_FILE=$(basename ${L001_FILE})
        L002_FILE=${L001_FILE/L001/L002}
        L003_FILE=${L001_FILE/L001/L003}
        L004_FILE=${L001_FILE/L001/L004}
        if [ ! -e ${FQ_PATH}${L003_FILE} ];then
            L003_FILE=''
        fi
        if [ ! -e ${FQ_PATH}${L004_FILE} ];then
            L004_FILE=''
        fi
        FILE=${L001_FILE/L001_/}
        echo "L001 File: ${L001_FILE}" >> ${OUTPUT_PATH}"log.txt"
        echo "L002 File: ${L002_FILE}" >> ${OUTPUT_PATH}"log.txt"
        echo "L003 File: ${L003_FILE}" >> ${OUTPUT_PATH}"log.txt"
        echo "L004 File: ${L004_FILE}" >> ${OUTPUT_PATH}"log.txt"

        #rename files for CellRanger
        [[ ${FILE} =~ _[0-4]_ ]];
        BEGIN_FILE=${FILE%${BASH_REMATCH}*};
        [[ ${BASH_REMATCH} =~ [0-4]_ ]];
        NB=${BASH_REMATCH};
        [[ ${FILE} =~ _S[0-9]+_ ]];
        S=${BASH_REMATCH};
        END_FILE=${FILE#*${BASH_REMATCH}};
        NEW_NAME=${BEGIN_FILE}${S}'L00'${NB}${END_FILE};
        echo "New File name: ${NEW_NAME}" >> ${OUTPUT_PATH}"log.txt"
        echo "#############" >> ${OUTPUT_PATH}"log.txt"
        #cat
        if [ -z ${L003_FILE} ] && [ -z ${L004_FILE} ]; then
          srun cat ${FQ_PATH}${L001_FILE} ${FQ_PATH}${L002_FILE} > ${OUTPUT_PATH}${NEW_NAME} &
        else
          srun cat ${FQ_PATH}${L001_FILE} ${FQ_PATH}${L002_FILE} ${FQ_PATH}${L003_FILE} ${FQ_PATH}${L004_FILE}> ${OUTPUT_PATH}${NEW_NAME} &
        fi

      done
    fi
  fi
done
echo "Finish to concat files"
