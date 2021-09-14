#!/usr/bin/env bash
set -e

########################################################################
## Script to automatically launch single-cell pipeline
##
## using: bash /mnt/beegfs/pipelines/single-cell/.BiGR/run.sh -adt -tcr -bcr -human -mouse -qc_only
########################################################################

#### Conda environments creation/loading ####
source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
#env for pipeline
echo "Loading 'single-cell-pipeline' environment"
conda activate /mnt/beegfs/userdata/${USER}/.environnement_conda/single-cell-pipeline || {
  echo "'single-cell-pipeline' environment not found! Creating the environment..."
  mkdir -p /mnt/beegfs/userdata/${USER}/.environnement_conda
  conda env create -f /mnt/beegfs/pipelines/single-cell/envs/conda/single-cell_user.yaml  --prefix="/mnt/beegfs/userdata/${USER}/.environnement_conda/single-cell-pipeline"
  conda activate /mnt/beegfs/userdata/${USER}/.environnement_conda/single-cell-pipeline || {
    echo "'single-cell-pipeline' environment can't be launched..."
    exit 1
  }
}
conda deactivate
#env for ADT
if [[ $* == *-adt* ]]
then
  echo "Loading 'kallisto' environment"
  conda activate /mnt/beegfs/userdata/${USER}/.environnement_conda/kallisto || {
    echo "'kallisto' environment not found! Creating the environment..."
    conda create -y --prefix="/mnt/beegfs/userdata/${USER}/.environnement_conda/kallisto" -c bioconda kallisto
    conda activate /mnt/beegfs/userdata/${USER}/.environnement_conda/kallisto || {
      echo "'kallisto' environment can't be launched..."
      exit 1
    }
  }
  conda deactivate
fi

#### Get path and samples names ####
project_path=$(dirname $(pwd))
data_input_folder=${project_path}'/data_input/'
data_output_folder=${project_path}'/data_output/'
script_folder=${project_path}'/script/'
ADT_file=${project_path}'/script/ADT.csv' #a csv file with protein_name;tag;gene_name, with header
tmp_folder=${project_path}'/tmp/'
yaml_file=${script_folder}'Params.yaml'
#remplacer "_S[1-4]" par "_[1-4]_S[1-4]"
samples_names=$(ls ${data_input_folder} | sed 's/\_S[1-4]/%/'|cut -d'%' -f1 | uniq) #from pf
#create folders if not exist
mkdir -p ${data_output_folder} ${script_folder} ${tmp_folder}

#### Creating ADT index ####
if [[ $* == *-adt* ]]
then
  conda activate /mnt/beegfs/userdata/${USER}/.environnement_conda/kallisto || {
    echo "'kallisto' environment can't be launched..."
    exit 1
  }
  echo "Creating the ADT index thanks to the ${ADT_file} file..."
  #creating folder
  ADT_INDEX_folfer=${data_output_folder}'ADT_INDEX/'
  mkdir -p ${ADT_INDEX_folfer}
  cd ${ADT_INDEX_folfer}
  #creating fasta file
  > "ADT.fa"
  protein_name_tag=$(cat ${ADT_file} | tail -n +2 | cut -d";" -f1,2)
  for i in ${protein_name_tag};do
    protein_name=$(echo ${i} | cut -d";" -f1)
    tag=$(echo ${i} | cut -d";" -f2)
    echo ">"${protein_name} >> "ADT.fa"
    echo ${tag} >> "ADT.fa"
  done
  #creating tr2gs file
  protein_name=$(cat ${ADT_file} | tail -n +2 | cut -d";" -f1)
  echo "transcript gene gene_name" > "ADT_tr2gs.txt"
  for i in ${protein_name};do
    echo $i" "$i" "$i >> "ADT_tr2gs.txt"
  done
  #creating index
  kallisto index 'ADT.fa' -i 'ADT.kidx' --kmer-size=15
  cd ${project_path}'/logs'
  # si pas -qc_only, get gene name from ADT.csv file
  if [[ ! $* == *-qc_only* ]]
  then
    ADT_gene_name=$(cat ${ADT_file} | tail -n +2 | cut -d";" -f3)
    ADT_gene_name=$(echo $ADT_gene_name | sed 's/ /,/g')
  fi
fi

#### Make parameters configfile (Params.yaml) ####
echo "Creating parameters configfile..."
steps='Steps: ['
echo 'Tmp: "'${tmp_folder}'"'> ${yaml_file}
echo ''>> ${yaml_file}

# For Gene Expression
if [[ $* == *-adt* ]] || [[ $* == *-tcr* ]] || [[ $* == *-bcr* ]]
then
  samples_names_ge=$(for i in $samples_names; do echo $i; done | grep "_GE$")
else
  samples_names_ge=${samples_names}
fi
if [[ -z "${samples_names_ge//}" ]] #s'il reste des noms d'échantillons on poursuit sinon on arrête avec une erreur.
then
  echo "No sample name found for Gene Expession!"
  exit 1
else
  steps=${steps}'"Alignment_countTable_GE","Droplets_QC_GE"'
  echo 'Alignment_countTable_GE:'>> ${yaml_file}
  echo '  sample.name.ge : ["'$(echo ${samples_names_ge} | sed 's/ /\",\"/g')'"]'>> ${yaml_file}
  echo '  input.dir.ge : "'${data_input_folder}'"'>> ${yaml_file}
  echo '  output.dir.ge : "'${data_output_folder}'"'>> ${yaml_file}
  if [[ $* == *-tcr* ]] || [[ $* == *-bcr* ]]
  then
    echo '  sctech : "10xv2"'>> ${yaml_file}
  else
    echo '  sctech : "10xv3"'>> ${yaml_file}
  fi
  if [[ $* == *-human* ]] && [[ $* == *-mouse* ]]
  then
    echo "Only one species must be set! (-human or -mouse)"
    exit 1
  elif [[ $* == *-human* ]]
  then
    echo '  kindex.ge : "/mnt/beegfs/database/bioinfo/single-cell/INDEX/KB-python_KALLISTO/0.24.4_0.46.2/homo_sapiens/GRCh38/Ensembl/r99/cDNA_LINCs_MIRs/GRCH38_r99_cDNA_linc_mir.kidx"'>> ${yaml_file}
    echo '  tr2g.file.ge : "/mnt/beegfs/database/bioinfo/single-cell/INDEX/KB-python_KALLISTO/0.24.4_0.46.2/homo_sapiens/GRCh38/Ensembl/r99/cDNA_LINCs_MIRs/GRCH38_r99_cDNA_linc_mir_tr2gs.txt"'>> ${yaml_file}
    echo '  reference.txt: "Ensembl reference transcriptome v99 corresponding to the homo sapiens GRCh38 build"'>> ${yaml_file}
    echo ''>> ${yaml_file}
    echo 'Droplets_QC_GE:'>> ${yaml_file}
    echo '  species: "homo_sapiens"'>> ${yaml_file}
    echo '  author.name: '${USER}>> ${yaml_file}
  elif [[ $* == *-mouse* ]]
  then
    echo '  kindex.ge : "/mnt/beegfs/database/bioinfo/single-cell/INDEX/KB-python_KALLISTO/0.24.4_0.46.2/mus_musculus/GRCm38/Ensembl/r99/cDNA_LINCs_MIRs/GRCm38_r99_cDNA_linc_mir.kidx"'>> ${yaml_file}
    echo '  tr2g.file.ge : "/mnt/beegfs/database/bioinfo/single-cell/INDEX/KB-python_KALLISTO/0.24.4_0.46.2/mus_musculus/GRCm38/Ensembl/r99/cDNA_LINCs_MIRs/GRCm38_r99_cDNA_linc_mir_tr2gs.txt"'>> ${yaml_file}
    echo '  reference.txt: "Ensembl reference transcriptome v99 corresponding to the homo sapiens GRCm38 build"'>> ${yaml_file}
    echo ''>> ${yaml_file}
    echo 'Droplets_QC_GE:'>> ${yaml_file}
    echo '  species: "mus_musculus"'>> ${yaml_file}
    echo '  author.name: "'${USER}'"'>> ${yaml_file}
  else
    echo "Species option mandatory! (-human or -mouse)"
    exit 1
  fi
  echo ''>> ${yaml_file}
  if [[ ! $* == *-qc_only* ]]
  then
    steps=${steps}',"Filtering_GE","Norm_DimRed_Eval_GE","Clust_Markers_Annot_GE"'
    echo 'Clust_Markers_Annot_GE:'>> ${yaml_file}
    echo '  keep.dims: 21'>> ${yaml_file}
    echo '  keep.res: 0.8'>> ${yaml_file}
    echo ''>> ${yaml_file}
  fi
fi

# For ADT -adt
if [[ $* == *-adt* ]]
then
  samples_names_adt=$(for i in $samples_names; do echo $i; done | grep "_ADT$")
  #s'il reste des noms d'échantillons on poursuit sinon on arrête avec une erreur.
  if [[ -z "${samples_names_adt//}" ]]
  then
    echo "No sample name found for Antibody-Derived Tags! Don't use the -adt option if there are no ADT fastq files."
    exit 1
  else
    steps=${steps}',"Alignment_countTable_ADT"'
    echo 'Alignment_countTable_ADT:'>> ${yaml_file}
    echo '  sample.name.adt : ["'$(echo ${samples_names_adt} | sed 's/ /\",\"/g')'"]'>> ${yaml_file}
    echo '  input.dir.adt : "'${data_input_folder}'"'>> ${yaml_file}
    echo '  output.dir.adt : "'${data_output_folder}'"'>> ${yaml_file}
    if [[ $* == *-tcr* ]] || [[ $* == *-bcr* ]]
    then
      echo '  sctech : "10xv2"'>> ${yaml_file}
    else
      echo '  sctech : "10xv3"'>> ${yaml_file}
    fi
    ### TO DO index ()faire l'index automatiquement à partir du .fa? + Adding_ADT
    echo '  kindex.adt: "'${ADT_INDEX_folfer}'ADT.kidx"'>> ${yaml_file}
    echo '  tr2g.file.adt: "'${ADT_INDEX_folfer}'ADT_tr2gs.txt"'>> ${yaml_file}
    echo ''>> ${yaml_file}
    if [[ ! $* == *-qc_only* ]]
    then
      steps=${steps}',"Adding_ADT"'
      echo 'Adding_ADT:'>> ${yaml_file}
      echo '  gene.names: "'${ADT_gene_name}'"'>> ${yaml_file}
      echo ''>> ${yaml_file}
    fi

  fi

fi

# For TCR/BCR -tcr -bcr
if [[ $* == *-tcr* ]] || [[ $* == *-bcr* ]]
then
  samples_names_tcr=$(for i in $samples_names; do echo $i; done | grep "_TCR$")
  samples_names_bcr=$(for i in $samples_names; do echo $i; done | grep "_BCR$")
  #s'il reste des noms d'échantillons on poursuit sinon on arrête avec une erreur.
  if [[ -z "${samples_names_tcr//}" ]] || [[ -z "${samples_names_bcr//}" ]]
  then
    echo "No sample name found for Immune profiling! Don't use the -tcr/-bcr option if there are no TCR/BCR fastq files."
    exit 1
  else
    steps=${steps}',"Alignment_annotations_TCR_BCR"'
    echo 'Alignment_annotations_TCR_BCR:'>> ${yaml_file}

    if [[ $* == *-tcr* ]]
    then
      #s'il reste des noms d'échantillons on poursuit sinon on arrête avec une erreur.
      if [[ -z "${samples_names_tcr//}" ]]
      then
        echo "No sample name found for TCR! Don't use the -tcr option if there are no TCR fastq files."
        exit 1
      else
        echo '  sample.name.tcr : ["'$(echo ${samples_names_tcr} | sed 's/ /\",\"/g')'"]'>> ${yaml_file}
        echo '  input.dir.tcr : "'${data_input_folder}'"'>> ${yaml_file}
        if [[ ! $* == *-qc_only* ]]
        then
            steps=${steps}',"Adding_TCR"'
        fi
      fi
    fi
    if [[ $* == *-bcr* ]]
    then
      #s'il reste des noms d'échantillons on poursuit sinon on arrête avec une erreur.
      if [[ -z "${samples_names_bcr//}" ]]
      then
        echo "No sample name found for BCR! Don't use the -bcr option if there are no BCR fastq files."
        exit 1
      else
        echo '  sample.name.bcr : ["'$(echo ${samples_names_bcr} | sed 's/ /\",\"/g')'"]'>> ${yaml_file}
        echo '  input.dir.bcr : "'${data_input_folder}'"'>> ${yaml_file}
        if [[ ! $* == *-qc_only* ]]
        then
            steps=${steps}',"Adding_BCR"'
        fi
      fi
    fi
    echo '  output.dir.tcr_bcr : "'${data_output_folder}'"'>> ${yaml_file}
    if [[ $* == *-human* ]] && [[ $* == *-mouse* ]]
    then
      echo "Only one species must be set! (-human or -mouse)"
      exit 1
    elif [[ $* == *-human* ]]
    then
      echo '  crindex.tcr_bcr : "/mnt/beegfs/database/bioinfo/single-cell/TCR_REFERENCES/refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0"'>> ${yaml_file}
    elif [[ $* == *-mouse* ]]
    then
      echo '  crindex.tcr_bcr : "/mnt/beegfs/database/bioinfo/single-cell/TCR_REFERENCES/refdata-cellranger-vdj-GRCm38-alts-ensembl-3.1.0"'>> ${yaml_file}
    else
      echo "Species option mandatory! (-human or -mouse)"
      exit 1
    fi

  fi

fi

# Add Cerebro step
if [[ ! $* == *-qc_only* ]]
then
    steps=${steps}',"Cerebro"]'
else
  steps=${steps}']'
fi

# Write Steps at the beginning
sed -i -e '1 i'"${steps}" ${yaml_file}


#### Run single-cell pipeline ####
echo "Run single-cell pipeline"
cat << EOT > ${script_folder}'launcher.sh'
#!/bin/bash

#SBATCH --job-name=pipeline_sc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --partition=mediumq

source /mnt/beegfs/software/conda/etc/profile.d/conda.sh
conda activate /mnt/beegfs/userdata/${USER}/.environnement_conda/single-cell-pipeline
module load singularity

#print environment tools versions
python --version
echo "snakemake version:"
snakemake --version
singularity --version

#launch
snakemake --profile /mnt/beegfs/pipelines/single-cell/profiles/slurm -s /mnt/beegfs/pipelines/single-cell/Snakefile --configfile ${yaml_file}

conda deactivate

EOT

sbatch ${script_folder}'launcher.sh'

