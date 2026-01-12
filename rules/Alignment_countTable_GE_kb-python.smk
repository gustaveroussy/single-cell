"""
##########################################################################
These rules make the alignment of genes expression in single-cell RNA-seq.
##########################################################################
"""

"""
This rule makes the copy of fastq files with the good sample name.
"""
def copy_rename_inputs_ge(wildcards):
    for i in range(0,len(ALIGN_SAMPLE_NAME_GE),1):
        if ALIGN_SAMPLE_NAME_GE[i] == wildcards.sample_name_ge :
            return os.path.normpath(ALIGN_INPUT_DIR_GE_RAW + "/" + ALIGN_SAMPLE_NAME_GE_RAW[i] + str("{lane_R_complement}.fastq.gz"))

rule copy_rename_fq_ge:
    input:
        fq = copy_rename_inputs_ge
    output:
        fq_copied = temp(os.path.normpath(ALIGN_INPUT_DIR_GE + "/{sample_name_ge}{lane_R_complement}.fastq.gz"))
    log:
        "logs/copy_rename_fq_ge/{sample_name_ge}{lane_R_complement}.log",
    benchmark:
        "benchmark/copy_rename_fq_ge/{sample_name_ge}{lane_R_complement}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 256, 2048)),
        time_min = (lambda wildcards, attempt: min(attempt * 5, 50))
    shell:
        """
        rsync -az -c --quiet {input.fq} {output.fq_copied} &> {log}
        """

"""
This rule makes the fastqc control-quality.
"""
rule fastqc_ge:
    input:
        fq_copied = os.path.normpath(ALIGN_INPUT_DIR_GE + "/{sample_name_ge}{lane_R_complement}.fastq.gz")
    output:
        html_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}{lane_R_complement}_fastqc.html")),
        zip_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}{lane_R_complement}_fastqc.zip"))
    log:
        "logs/fastqc_ge/{sample_name_ge}{lane_R_complement}.log",
    benchmark:
        "benchmark/fastqc_ge/{sample_name_ge}{lane_R_complement}.tsv"
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
        fastqc --outdir {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads --dir $(mktemp -d $TMPDIR --threads {threads} {input.fq_copied} &> {log}
        """

"""
This rule makes the fastq-screen control-quality on R2 files.
"""
rule fastqscreen_ge:
    input:
        R2_fq = os.path.normpath(ALIGN_INPUT_DIR_GE + "/{sample_name_ge}{lane_R_complement}.fastq.gz")
    output:
        html_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}{lane_R_complement}_screen.html")),
        txt_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}{lane_R_complement}_screen.txt")),
        png_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}{lane_R_complement}_screen.png"))
    log:
        "logs/fastqscreen_ge/{sample_name_ge}{lane_R_complement}.log",
    benchmark:
        "benchmark/fastqscreen_ge/{sample_name_ge}{lane_R_complement}.tsv"
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: min(2048 + attempt * 2048, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_FASTQ_SCREEN
    shell:
        """
        fastq_screen --quiet --threads {threads} --force --outdir {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads --subset 100000 --conf {FASTQSCREEN_INDEX} {input.R2_fq} &> {log}
        """

"""
This rule makes the multiqc from the fastqc and the fastq-screen results.
The function allows to get all QC input files for one specific sample (wildcards).
"""
def multiqc_inputs_ge(wildcards):
    name_R1_R2=[elem for elem in ALIGN_SYMLINK_FILES_NAME_GE if re.search(str( "^" + wildcards.sample_name_ge), elem)]
    name_R2=[elem for elem in name_R1_R2 if re.search("R2", elem)]
    print("name_R1_R2: ")
    print(name_R1_R2)
    print("name_R2 :")
    print(name_R2)
    files=[]
    for name in name_R1_R2:
        #fastqc
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/" + wildcards.sample_name_ge + "/QC_reads/" + name + "_fastqc.html"))
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/" + wildcards.sample_name_ge + "/QC_reads/" + name + "_fastqc.zip"))
    for name in name_R2:
        #fastqscreen
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/" + wildcards.sample_name_ge + "/QC_reads/" + name + "_screen.html"))
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/" + wildcards.sample_name_ge + "/QC_reads/" + name + "_screen.txt"))
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/" + wildcards.sample_name_ge + "/QC_reads/" + name + "_screen.png"))
    print(files)
    return files

rule multiqc_ge:
    input:
        qc_files = multiqc_inputs_ge
    output:
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}_RAW.html")
    log:
        "logs/multiqc_ge/{sample_name_ge}.log",
    benchmark:
        "benchmark/multiqc_ge/{sample_name_ge}.tsv"
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
            -n {wildcards.sample_name_ge}'_RAW' \
            -i {wildcards.sample_name_ge}' RAW FASTQ' \
            -f --no-megaqc-upload --no-data-dir \
            -o {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads {input.qc_files} &> {log}
        """

"""
This rule makes the alignment, sort, correct, and count UMI by kb-python (wrapper of kallisto and bustools)
The function alignment_inputs_ge allows to get all fastq input files for one specific sample (wildcards).
"""
def alignment_inputs_ge(wildcards):
    return sorted([elem for elem in ALIGN_SYMLINK_FILES_GE if re.search(wildcards.sample_name_ge, elem)])

rule alignment_build_count_matrix_ge:
    input:
        fq_copied = alignment_inputs_ge,
        tr2g = TR2GFILE_GE,
        index = KINDEX_GE,
        whitelist = WHITELISTNAME
    output:
        mtx_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.mtx"),
        barcodes_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.barcodes.txt"),
        genes_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.genes.txt"),
        MandM = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/Materials_and_Methods.txt")
    params:
        kallistobus_path = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS")
    log:
        "logs/alignment_build_count_matrix_ge/{sample_name_ge}.log",
    benchmark:
        "benchmark/alignment_build_count_matrix_ge/{sample_name_ge}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 40960, 204800)),
        time_min = (lambda wildcards, attempt: min(attempt * 120, 600))
    conda:
        CONDA_ENV_KB_PYTHON
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        res=$(({resources.mem_mb}-500))
        mkdir -p {params.kallistobus_path} && \
        kb count \
        -i {input.index} \
        -g {input.tr2g} \
        -x {SCTECH} \
        -w {input.whitelist} \
        -o {params.kallistobus_path}/{wildcards.sample_name_ge} \
        --tmp $TMPDIR \
        -t {threads} \
        -m $res"M" \
        --h5ad --overwrite \
        {KB_PYTHON_EXTRA_GE} \
        {input.fq_copied} &> {log}

        FASTQC_V=$(grep fastqc {PIPELINE_FOLDER}/envs/conda/fastqc.yaml | cut -d= -f2)
        FASTQSCREEN_V=$(grep fastq-screen {PIPELINE_FOLDER}/envs/conda/fastq-screen.yaml | cut -d= -f2)
        KBPYTHON_V=$(grep kb-python {PIPELINE_FOLDER}/envs/conda/kb-python.yaml | cut -d= -f2)
        if [[ {SCTECH} = '10xv3' ]];then
            CR="10X Chromium 3' scRNA-Seq v3 chemistry"
        elif [[ {SCTECH} = '10xv2' ]];then
            CR="10X Chromium 5' scRNA-Seq v2 chemistry"
        fi
        echo "Raw BCL-files were demultiplexed and converted to Fastq format using bcl2fastq (version 2.20.0.422 from Illumina).
Reads quality control was performed using fastqc (version $FASTQC_V) and assignment to the expected genome species evaluated with fastq-screen (version $FASTQSCREEN_V).
An index with the {REF_TXT_GE} was made with the kb-python (version $KBPYTHON_V). Index built, reads pseudo-mapping with parameter corresponding to the $CR, barcode correction using whitelist provided by the manufacturer (10X Genomics), and gene-based reads quantification, were performed by kb-python." > {output.MandM}
        """
