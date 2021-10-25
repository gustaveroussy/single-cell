#!/bin/bash
# script to create cellranger scripts (one by sample), and launch it on flamingo cluster.
# usage: 
# CR_input_fq="path_to_fastq_files"
# CR_output="path_to_output"
# script_folder="path_to_script_folder"
# bash create_run_scripts_for_cellranger.sh "name_sample1" "name_sample2" "name_sample3"

all_samples=${*}

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
--transcriptome=/mnt/beegfs/database/bioinfo/cellranger/3.0.0/refdata-cellranger-GRCh38-3.0.0 \
--fastqs=${CR_input_fq} \
--sample=${sample} \
--expect-cells=10000 \
--localcores=4 \
--localmem=32

EOT

done

sbatch ${script_folder}'/launcher_'${sample}'.sh'
