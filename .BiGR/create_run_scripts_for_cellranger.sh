#!/bin/bash
# script to create cellranger scripts (one by sample), and launch it on flamingo cluster.
# usage: 
# export CR_input_fq="path_to_fastq_files"
# export CR_output="path_to_output"
# export script_folder="path_to_script_folder"
# export transcriptome="/mnt/beegfs/database/bioinfo/cellranger/3.0.0/refdata-cellranger-mm10-3.0.0" #human by default
# bash create_run_scripts_for_cellranger.sh "name_sample1" "name_sample2" "name_sample3"

all_samples=${*}

if [ -z "$transcriptome" ]; then transcriptome="/mnt/beegfs/database/bioinfo/cellranger/3.0.0/refdata-cellranger-GRCh38-3.0.0"; fi

for sample in ${all_samples}
do
echo ${sample}

cat << EOT > ${script_folder}'/launcher_'${sample}'.sh'
#!/bin/bash
#SBATCH --job-name=CR_${sample}
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=mediumq

mkdir -p ${CR_output}
cd ${CR_output}

/mnt/beegfs/software/cellranger/3.1.0/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger count \
--id=${sample} \
--transcriptome=${transcriptome} \
--fastqs=${CR_input_fq} \
--sample=${sample} \
--expect-cells=10000 \
--localcores=4 \
--localmem=32

EOT

sbatch ${script_folder}'/launcher_'${sample}'.sh'

done

