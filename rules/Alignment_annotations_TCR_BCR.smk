"""
##########################################################################
These rules make the alignment of TCR and BCR in single-cell RNA-seq.
##########################################################################
"""

"""
This rule makes the copy of fastq files with the good sample name.
"""
def copy_rename_inputs_tcr_bcr(wildcards):
    res = "no result find!"
    to_test = os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + wildcards.lane_R_complement + ".fastq.gz")
    for i in range(0,len(ALIGN_SYMLINK_FILES_TCR_BCR),1):
        if os.path.normpath(ALIGN_SYMLINK_FILES_TCR_BCR[i]) == to_test :
            res = str(ALIGN_ORIG_FILES_TCR_BCR[i])
    return res

rule copy_rename_fq_tcr_bcr:
    input:
        fq = copy_rename_inputs_tcr_bcr
    output:
        fq_copied = temp(os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}{lane_R_complement}.fastq.gz"))
    log:
        "logs/copy_rename_fq_tcr_bcr/{sample_name_tcr_bcr}{lane_R_complement}.log"
    benchmark:
        "benchmark/copy_rename_fq_tcr_bcr/{sample_name_tcr_bcr}{lane_R_complement}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 256, 2048)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    shell:
        """
        rsync -az -c --quiet {input.fq} {output.fq_copied} &> {log}
        """

"""
This rule makes the fastqc control-quality.
"""
rule fastqc_tcr_bcr:
    input:
        fq_copied = os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}{lane_R_complement}.fastq.gz")
    output:
        html_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/QC_reads/{sample_name_tcr_bcr}{lane_R_complement}_fastqc.html")),
        zip_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/QC_reads/{sample_name_tcr_bcr}{lane_R_complement}_fastqc.zip"))
    log:
        "logs/fastqc_tcr_bcr/{sample_name_tcr_bcr}{lane_R_complement}.log"
    benchmark:
        "benchmark/fastqc_tcr_bcr/{sample_name_tcr_bcr}{lane_R_complement}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_FASTQC
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        fastqc --outdir {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads --dir $TMPDIR --threads {threads} {input.fq_copied} &> {log}
        """

"""
This rule makes the fastq-screen control-quality on R2 files.
"""
rule fastqscreen_tcr_bcr:
    input:
        R2_fq = os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}{lane_R_complement}.fastq.gz")
    output:
        html_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR +"/{sample_name_tcr_bcr}/QC_reads/{sample_name_tcr_bcr}{lane_R_complement}_screen.html")),
        txt_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR +"/{sample_name_tcr_bcr}/QC_reads/{sample_name_tcr_bcr}{lane_R_complement}_screen.txt")),
        png_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR +"/{sample_name_tcr_bcr}/QC_reads/{sample_name_tcr_bcr}{lane_R_complement}_screen.png"))
    log:
        "logs/fastqscreen_tcr_bcr/{sample_name_tcr_bcr}{lane_R_complement}.log"
    benchmark:
        "benchmark/fastqscreen_tcr_bcr/{sample_name_tcr_bcr}{lane_R_complement}.tsv"
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: min(2048 + attempt * 2048, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_FASTQ_SCREEN
    shell:
        """
        fastq_screen --quiet --threads {threads} --force --outdir {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads --subset 100000 --conf {FASTQSCREEN_INDEX} {input.R2_fq} &> {log}
        """

"""
This rule makes the multiqc from the fastqc and the fastq-screen results.
The function allows to get all QC input files for one specific sample (wildcards).
"""
def multiqc_inputs_tcr_bcr(wildcards):
    name_R1_R2=[elem for elem in ALIGN_SYMLINK_FILES_NAME_TCR_BCR if re.search(wildcards.sample_name_tcr_bcr, elem)]
    name_R2=[elem for elem in name_R1_R2 if re.search("R2", elem)]
    files=[]
    for name in name_R1_R2:
        #fastqc
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + "/QC_reads/" + name) + "_fastqc.html")
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + "/QC_reads/" + name) + "_fastqc.zip")
    for name in name_R2:
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + "/QC_reads/" + name) + "_screen.html")
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + "/QC_reads/" + name) + "_screen.txt")
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + "/QC_reads/" + name) + "_screen.png")
    return files

rule multiqc_tcr_bcr:
    input:
        qc_files = multiqc_inputs_tcr_bcr
    output:
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/QC_reads/{sample_name_tcr_bcr}_RAW.html"),
    log:
        "logs/multiqc_tcr_bcr/{sample_name_tcr_bcr}.log"
    benchmark:
        "benchmark/multiqc_tcr_bcr/{sample_name_tcr_bcr}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_MULTIQC
    shell:
        """
        export TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        multiqc \
            -n {wildcards.sample_name_tcr_bcr}'_RAW' \
            -i {wildcards.sample_name_tcr_bcr}' RAW FASTQ' \
            -f --no-megaqc-upload --no-data-dir \
            -o {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads/ {input.qc_files} &> {log}
        """

"""
This function allows to determine the singularity binding parameters.
"""
def alignment_annotations_tcr_bcr_params_sing(wildcards):
    input_folder = os.path.dirname(ALIGN_SYMLINK_FILES_TCR_BCR[0])
    output_folder = os.path.dirname(ALIGN_OUTPUT_DIR_TCR_BCR + "/wildcards.sample_name_tcr_bcr")
    concat = " -B " + PIPELINE_FOLDER + "," + input_folder + "," + output_folder
    return concat

"""
This rule makes the alignment and annotation for vdj by cellranger.
The function alignment_annotations_inputs_tcr_bcr allows to get all fastq input files for one specific sample (wildcards).
"""
def alignment_annotations_inputs_tcr_bcr(wildcards):
    files=[]
    files=[elem for elem in ALIGN_SYMLINK_FILES_TCR_BCR if re.search(wildcards.sample_name_tcr_bcr, elem)]
    return sorted(files)

rule alignment_annotations_tcr_bcr:
    input:
        fq = alignment_annotations_inputs_tcr_bcr
    output:
        csv_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/filtered_contig_annotations.csv"),
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/web_summary.html"),
        MandM = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/Materials_and_Methods.txt")
    log:
        "logs/alignment_annotations_tcr_bcr/{sample_name_tcr_bcr}.log"
    benchmark:
        "benchmark/alignment_annotations_tcr_bcr/{sample_name_tcr_bcr}.tsv"
    params:
        sing_bind = alignment_annotations_tcr_bcr_params_sing,
        sample_folder = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}")
    threads:
        3
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 5120, 25600)),
        time_min = (lambda wildcards, attempt: min(5*60 + attempt * 12*60, 7*24*60))
    shell:
        """
        rm -r {params.sample_folder}/{wildcards.sample_name_tcr_bcr}_CellRanger
        echo "ressources: " &> {log}
        echo {resources.mem_mb} &>> {log}

        echo 'cd {params.sample_folder} && \
        /softwares/cellranger-9.0.1/cellranger vdj \
                 --id={wildcards.sample_name_tcr_bcr}_CellRanger \
                 --reference={CRINDEX_TCR_BCR} \
                 --fastqs={ALIGN_INPUT_DIR_TCR_BCR} \
                 --sample={wildcards.sample_name_tcr_bcr} \
                 --localmem=$res \
                 --localcores={threads} \
                 --disable-ui' | singularity exec --env res=$(({resources.mem_mb}/1000)) --contain {params.sing_bind} {SINGULARITY_ENV_CR} bash &>> {log}
        rm -r {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/{wildcards.sample_name_tcr_bcr}_CellRanger/SC_VDJ_ASSEMBLER_CS
        rm -r {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/{wildcards.sample_name_tcr_bcr}_CellRanger/extras
        rm {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/{wildcards.sample_name_tcr_bcr}_CellRanger/_*
        rm {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/{wildcards.sample_name_tcr_bcr}_CellRanger/*_CellRanger.mri.tgz

        FASTQC_V=$(grep fastqc {PIPELINE_FOLDER}/envs/conda/fastqc.yaml | cut -d= -f2)
        FASTQSCREEN_V=$(grep fastq-screen {PIPELINE_FOLDER}/envs/conda/fastq-screen.yaml | cut -d= -f2)
        CELLRANGER_V="9.0.1"
        echo "Reads quality control was performed using fastqc (version $FASTQC_V) and assignment to the expected genome species evaluated with fastq-screen (version $FASTQSCREEN_V).
CellRanger (version $CELLRANGER_V from 10X Genomics) was used to generate single-cell V(D)J sequences and annotations." > {output.MandM}
        """
