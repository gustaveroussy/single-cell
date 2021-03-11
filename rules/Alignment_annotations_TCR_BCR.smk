"""
##########################################################################
These rules make the alignment of TCR and BCR in single-cell RNA-seq.
##########################################################################
"""
wildcard_constraints:
    sample_name_tcr_bcr_R=".+(_TCR|_BCR)",
    sample_name_tcr_bcr=".+(_TCR|_BCR)"

"""
This rule makes the fastqc control-quality.
"""
rule fastqc_tcr_bcr:
    input:
        fq = os.path.join(ALIGN_INPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr_R}{lane_R_complement}.fastq.gz")
    output:
        html_file = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr_R}/QC/fastqc/{sample_name_tcr_bcr_R}{lane_R_complement}_fastqc.html"),
        zip_file = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr_R}/QC/fastqc/{sample_name_tcr_bcr_R}{lane_R_complement}_fastqc.zip")
    threads:
        8
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "mkdir -p {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr_R}/QC/fastqc && fastqc --quiet -o {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr_R}/QC/fastqc -t {threads} {input}"


"""
This rule makes the fastq-screen control-quality on R2 files.
"""
rule fastqscreen_tcr_bcr:
    input:
        R2_fq = os.path.join(ALIGN_INPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr_R}{lane_R_complement}.fastq.gz")
    output:
        html_file = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr_R}/QC/fastqscreen/{sample_name_tcr_bcr_R}{lane_R_complement}_screen.html"),
        #png_file = os.path.join(OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr_R}/QC/fastqscreen/{sample_name_tcr_bcr_R}{lane_R_complement}_screen.png"),
        txt_file = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr_R}/QC/fastqscreen/{sample_name_tcr_bcr_R}{lane_R_complement}_screen.txt")
    threads:
        8
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "mkdir -p {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr_R}/QC/fastqscreen && fastq_screen --quiet --threads {threads} --force --outdir {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr_R}/QC/fastqscreen --subset 100000 --conf {FASTQSCREEN_INDEX} {input.R2_fq}"


"""
This rule makes the multiqc from the fastqc and the fastq-screen results.
The function allows to get all QC input files for one specific sample (wildcards).
"""
def multiqc_inputs_tcr_bcr(wildcards):
    name_R1_R2=[elem for elem in ALL_FILES_TCR_BCR if re.search(wildcards.sample_name_tcr_bcr, elem)]
    name_R2=[elem for elem in name_R1_R2 if re.search("R2", elem)]
    files=[]
    for name in name_R1_R2:
        #fastqc
        files.append(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,wildcards.sample_name_tcr_bcr,"QC/fastqc",name) + "_fastqc.html")
        files.append(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,wildcards.sample_name_tcr_bcr,"QC/fastqc", name) + "_fastqc.zip")
    for name in name_R2:
        files.append(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,wildcards.sample_name_tcr_bcr,"QC/fastqscreen", name) + "_screen.html")
        #files.append(os.path.join(OUTPUT_DIR_TCR_BCR,wildcards.sample_name_tcr_bcr,"QC/fastqscreen", name) + "_screen.png")
        files.append(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,wildcards.sample_name_tcr_bcr,"QC/fastqscreen", name) + "_screen.txt")
    return files

rule multiqc_tcr_bcr:
    input:
        qc_files = multiqc_inputs_tcr_bcr
    output:
        html_file = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/QC/multiqc/{sample_name_tcr_bcr}_RAW.html"),
        zip_file = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/QC/multiqc/{sample_name_tcr_bcr}_RAW_data.zip")
    threads:
        1
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "mkdir -p {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC/multiqc && multiqc -n {wildcards.sample_name_tcr_bcr}'_RAW' -i {wildcards.sample_name_tcr_bcr}' RAW FASTQ' -p -z -f -o {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC/multiqc {input}"

"""
This rule makes the alignment and annotation by cellranger.
The function alignment_annotations_inputs_tcr_bcr allows to get all fastq input files for one specific sample (wildcards).
"""
def alignment_annotations_inputs_tcr_bcr(wildcards):
    files=[]
    files=[elem for elem in PATH_ALL_FILES_TCR_BCR_FQ_GZ if re.search(wildcards.sample_name_tcr_bcr, elem)]
    return sorted(files)

rule alignment_annotations_tcr_bcr:
    input:
        fq = alignment_annotations_inputs_tcr_bcr,
        html_file = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/QC/multiqc/{sample_name_tcr_bcr}_RAW.html")
    output:
        csv_file = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/filtered_contig_annotations.csv"),
        html_file = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/web_summary.html"),
	    MandM = os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/Materials_and_Methods.txt")
    params:
        sample_folder = directory(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}"))
    threads:
        8
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        """
        # /mnt/beegfs/software/cellranger/3.1.0/cellranger-3.1.0/sourceme.bash
        # /home/m_aglave/Softwares/cellranger-3.1.0/sourceme.bash
        cd {params}
        rm -r {wildcards.sample_name_tcr_bcr}_CellRanger
        /home/m_aglave/Softwares/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger vdj \
                 --id={wildcards.sample_name_tcr_bcr}_CellRanger \
                 --reference={CRINDEX_TCR_BCR} \
                 --fastqs={ALIGN_INPUT_DIR_TCR_BCR} \
                 --sample={wildcards.sample_name_tcr_bcr} \
                 --localmem=10 \
                 --localcores={threads} #&& \
        conda list
        FASTQC_V=$(conda list "fastqc" | grep "^fastqc " | sed -e "s/fastqc *//g" | sed -e "s/ .*//g")
        FASTQSCREEN_V=$(conda list "fastq-screen" | grep "^fastq-screen " | sed -e "s/fastq-screen *//g" | sed -e "s/ .*//g")
        CELLRANGER_V=`/home/m_aglave/Softwares/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger vdj --version | grep "cellranger vdj (" | sed -e "s/cellranger vdj (//g" | sed -e "s/)//g"`
        echo "Raw BCL-files were demultiplexed and converted to Fastq format using bcl2fastq (version 2.20.0.422 from Illumina). 
Reads quality control was performed using fastqc (version $FASTQC_V) and assignment to the expected genome species evaluated with fastq-screen (version $FASTQSCREEN_V).
CellRanger (version $CELLRANGER_V from 10X Genomics) was used to generate single-cell V(D)J sequences and annotations." > {output.MandM}
        """
