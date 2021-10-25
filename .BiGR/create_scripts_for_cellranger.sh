#!/bin/bash
#usage: 
#script_folder="path_to_script_folder"
#bash create_scripts_for_cellranger.sh

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
#SBATCH --partition=longq

cd /mnt/beegfs/scratch/bioinfo_core/B21060_NADR_10/CellRanger

/mnt/beegfs/software/cellranger/3.1.0/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger count \
--id=${sample} \
--transcriptome=/mnt/beegfs/database/bioinfo/cellranger/3.0.0/refdata-cellranger-GRCh38-3.0.0 \
--fastqs=/mnt/beegfs/scratch/bioinfo_core/B21060_NADR_10/CellRanger/0_rename_fastq \
--sample=${sample} \
--expect-cells=10000 \
--localcores=4 \
--localmem=32

EOT

done