#!/bin/bash
set -e
#usage:
# export  script_folder="/mnt/beegfs/scratch/m_aglave/P33_AMJO_1/script/"
# export  CR_input="/mnt/beegfs/scratch/m_aglave/P33_AMJO_1/data_input/"
# export  CR_output="/mnt/beegfs/scratch/m_aglave/P33_AMJO_1/data_output/CellRanger/"
# bash /mnt/beegfs/scratch/m_aglave/P33_AMJO_1/script/create_scripts_for_cellranger_multi_by_sample.sh --ge --ge_1_2_3_4 "sample_GE" --adt --adt_1_2_3_4 "sample_ADT" --tcr --tcr_1_2_3_4 "sample_TCR" --bcr --bcr_1_2_3_4 "sample_BCR" --species human/mouse

## Usage
usage()
{
cat << EOF
Usage: $0 options
This is a script to run CellRanger multi.
OPTIONS:
  --ge            Follow by the sample name of Gene Expression if there is one library.
  --ge_1_2_3_4    Follow by the sample name of Gene Expression if there are 4 libraries.
  --adt           Follow by the sample name of Antibody if there is one library.
  --adt_1_2_3_4   Follow by the sample name of Antibody if there are 4 libraries.
  --tcr           Follow by the sample name of TCR if there is one library.
  --tcr_1_2_3_4   Follow by the sample name of TCR if there is are 4 libraries.
  --bcr           Follow by the sample name of BCR if there is one library.
  --bcr_1_2_3_4   Follow by the sample name of BCR if there is are 4 libraries.
  --species       Follow by the species of the sample (human or mouse).
  --help          Print this message.
example:
#fatsq files are:
sample_GE_1_S13_L001_I1_001.fastq.gz
sample_GE_1_S13_L001_R1_001.fastq.gz
sample_GE_1_S13_L001_R2_001.fastq.gz
sample_GE_2_S13_L001_I1_001.fastq.gz
sample_GE_2_S13_L001_R1_001.fastq.gz
sample_GE_2_S13_L001_R2_001.fastq.gz
sample_GE_3_S13_L001_I1_001.fastq.gz
sample_GE_3_S13_L001_R1_001.fastq.gz
sample_GE_3_S13_L001_R2_001.fastq.gz
sample_GE_4_S13_L001_I1_001.fastq.gz
sample_GE_4_S13_L001_R1_001.fastq.gz
sample_GE_4_S13_L001_R2_001.fastq.gz
sample_ADT_S14_L001_I1_001.fastq.gz
sample_ADT_S14_L001_R1_001.fastq.gz
sample_ADT_S14_L001_R2_001.fastq.gz
sample_TCR_S14_L001_I1_001.fastq.gz
sample_TCR_S14_L001_R1_001.fastq.gz
sample_TCR_S14_L001_R2_001.fastq.gz
sample_BCR_S14_L001_I1_001.fastq.gz
sample_BCR_S14_L001_R1_001.fastq.gz
sample_BCR_S14_L001_R2_001.fastq.gz
#run
bash create_scripts_for_cellranger_multi_by_sample.sh --ge_1_2_3_4 "sample_GE" --adt "sample_ADT" --tcr "sample_TCR" --bcr "sample_BCR" --species human
#Then change the feature_ref.csv file if you have ADT data.
EOF
}

## Get arguments
all_arg=($( for i in ${*} ; do echo $i ; done | egrep -e "--ge|--ge_1_2_3_4|--adt|--adt_1_2_3_4|--tcr|--tcr_1_2_3_4|--bcr|--bcr_1_2_3_4" ))

POSSIBLE_ARG=ge::,ge_1_2_3_4::,adt::,adt_1_2_3_4::,tcr::,tcr_1_2_3_4::,bcr::,bcr_1_2_3_4::,species:,help
OPTS=$(getopt -n --longoptions $POSSIBLE_ARG -- "$@")
eval set -- "$OPTS"
echo $OPTS
while :
do
  case "$1" in
    "--ge" | "--ge_1_2_3_4" )
      sample_GE="$2"
      shift 2
      ;;
    "--adt" | "--adt_1_2_3_4" )
      sample_ADT="$2"
      shift 2
      ;;
    "--tcr" | "--tcr_1_2_3_4" )
      sample_TCR="$2"
      shift 2
      ;;
    "--bcr" | "--bcr_1_2_3_4" )
      sample_BCR="$2"
      shift 2
      ;;
    "--species" )
      species="$2"
      shift 2
      ;;
    "--help" )
      usage
      exit 0
      ;;
    --)
      shift 1
      ;;
    * )
    #  echo "Unexpected option: $1"
      break
      ;;
  esac
done

## Create output folder
mkdir -p ${CR_output}

## Set params
#references
if [[ ${species} == "human" ]];then
  ref_ge='/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-gex-GRCh38-2020-A'
  ref_vdj='/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0'
fi
if [[ ${species} == "mouse" ]];then
  ref_ge='/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-gex-mm10-2020-A'
  ref_vdj='/mnt/beegfs/database/bioinfo/cellranger/6.0.2/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0'
fi

## Make launcher script
echo "sample_GE: ${sample_GE}"
cat << EOT > ${script_folder}'CR_launcher_multi_'${sample_GE}'.sh'
#!/bin/bash
#SBATCH --job-name=CRm_${sample_GE}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=24G
#SBATCH --partition=longq

cd ${CR_output}

/mnt/beegfs/software/cellranger/6.0.2/cellranger-6.0.2/bin/cellranger multi \
--id=${sample_GE} \
--csv=${script_folder}CR_params_multi_${sample_GE}.csv \
--localcores=3 \
--localmem=24

EOT

## Make params script
#GE
cat << EOT > ${script_folder}'CR_params_multi_'${sample_GE}'.csv'
[gene-expression]
reference,${ref_ge}
expect-cells,10000
no-secondary,false
no-bam,true

EOT

#ADT
if [[ "${all_arg[@]}" =~ "--adt " || "${all_arg[@]}" =~ --adt$ || "${all_arg[@]}" =~ "--adt_1_2_3_4 " || "${all_arg[@]}" =~ --adt_1_2_3_4$ ]]; then
touch ${script_folder}feature_ref.csv
cat << EOT >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'
[feature]
reference,${script_folder}feature_ref.csv

EOT
fi

#VDJ
if [[ "${all_arg[@]}" =~ --tcr$ || "${all_arg[@]}" =~ "--tcr " || "${all_arg[@]}" =~ --bcr$ || "${all_arg[@]}" =~ "--bcr " || "${all_arg[@]}" =~ --tcr_1_2_3_4$ || "${all_arg[@]}" =~ "--tcr_1_2_3_4 " || "${all_arg[@]}" =~ --bcr_1_2_3_4$ || "${all_arg[@]}" =~ "--bcr_1_2_3_4 " ]]; then
cat << EOT >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'
[vdj]
reference,${ref_vdj}

EOT
fi

# Add libraries section
#header
cat << EOT >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'
[libraries]
fastq_id,fastqs,feature_types
EOT
#ge
if [[ "${all_arg[@]}" =~ --ge$ || "${all_arg[@]}" =~ "--ge " ]]; then echo "${sample_GE},${CR_input},gene expression" >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'; fi
if [[ "${all_arg[@]}" =~ --ge_1_2_3_4$ || "${all_arg[@]}" =~ "--ge_1_2_3_4 " ]]; then
cat << EOT >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'
${sample_GE}_1,${CR_input},gene expression
${sample_GE}_2,${CR_input},gene expression
${sample_GE}_3,${CR_input},gene expression
${sample_GE}_4,${CR_input},gene expression
EOT
fi
#adt
if [[ "${all_arg[@]}" =~ --adt$ || "${all_arg[@]}" =~ "--adt " ]]; then echo "${sample_ADT},${CR_input},antibody capture" >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'; fi
if [[ "${all_arg[@]}" =~ --adt_1_2_3_4$ || "${all_arg[@]}" =~ "--adt_1_2_3_4 " ]]; then
cat << EOT >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'
${sample_ADT}_1,${CR_input},antibody capture
${sample_ADT}_2,${CR_input},antibody capture
${sample_ADT}_3,${CR_input},antibody capture
${sample_ADT}_4,${CR_input},antibody capture
EOT
fi
#tcr
if [[ "${all_arg[@]}" =~ --tcr$ || "${all_arg[@]}" =~ "--tcr " ]]; then echo "${sample_TCR},${CR_input},vdj-t" >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'; fi
if [[ "${all_arg[@]}" =~ --tcr_1_2_3_4$ || "${all_arg[@]}" =~ "--tcr_1_2_3_4 " ]]; then
cat << EOT >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'
${sample_TCR}_1,${CR_input},vdj-t
${sample_TCR}_2,${CR_input},vdj-t
${sample_TCR}_3,${CR_input},vdj-t
${sample_TCR}_4,${CR_input},vdj-t
EOT
fi
#bcr
if [[ "${all_arg[@]}" =~ --bcr$ || "${all_arg[@]}" =~ "--bcr " ]]; then echo "${sample_BCR},${CR_input},vdj-b" >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'; fi
if [[ "${all_arg[@]}" =~ --bcr_1_2_3_4$ || "${all_arg[@]}" =~ "--bcr_1_2_3_4 " ]]; then
cat << EOT >> ${script_folder}'CR_params_multi_'${sample_GE}'.csv'
${sample_BCR}_1,${CR_input},vdj-b
${sample_BCR}_2,${CR_input},vdj-b
${sample_BCR}_3,${CR_input},vdj-b
${sample_BCR}_4,${CR_input},vdj-b
EOT
fi
