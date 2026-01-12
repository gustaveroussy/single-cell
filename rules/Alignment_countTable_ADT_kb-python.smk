"""
##########################################################################
These rules make the alignment of cell surface protein in single-cell RNA-seq.
##########################################################################
"""

"""
This rule makes the copy of fastq files with the good sample name.
"""
def copy_rename_inputs_adt(wildcards):
    for i in range(0,len(ALIGN_SAMPLE_NAME_ADT),1):
        if ALIGN_SAMPLE_NAME_ADT[i] == wildcards.sample_name_adt :
            return os.path.normpath(ALIGN_INPUT_DIR_ADT_RAW + "/" + ALIGN_SAMPLE_NAME_ADT_RAW[i] + str("{lane_R_complement}.fastq.gz"))

rule copy_rename_fq_adt:
    input:
        fq = copy_rename_inputs_adt
    output:
        fq_copied = temp(os.path.normpath(ALIGN_INPUT_DIR_ADT + "/{sample_name_adt}{lane_R_complement}.fastq.gz"))
    log:
        "logs/copy_rename_fq_adt/{sample_name_adt}{lane_R_complement}.log"
    benchmark:
        "benchmark/copy_rename_fq_adt/{sample_name_adt}{lane_R_complement}.tsv"
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
rule fastqc_adt:
    input:
        fq_copied = os.path.normpath(ALIGN_INPUT_DIR_ADT + "/{sample_name_adt}{lane_R_complement}.fastq.gz")
    output:
        html_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/QC_reads/{sample_name_adt}{lane_R_complement}_fastqc.html")),
        zip_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/QC_reads/{sample_name_adt}{lane_R_complement}_fastqc.zip"))
    log:
        "logs/fastqc_adt/{sample_name_adt}{lane_R_complement}.log"
    benchmark:
        "benchmark/fastqc_adt/{sample_name_adt}{lane_R_complement}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_FASTQC
    shell:
        """
        fastqc --outdir {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads --dir $(mktemp -d {resources.tmpdir}/XXXXXX) --threads {threads} {input.fq_copied} &> {log}
        """

"""
This rule makes the multiqc from the fastqc and the fastq-screen results.
The function allows to get all QC input files for one specific sample (wildcards).
"""
def multiqc_inputs_adt(wildcards):
    name_R1_R2=[elem for elem in ALIGN_SYMLINK_FILES_NAME_ADT if re.search(wildcards.sample_name_adt, elem)]
    files=[]
    for name in name_R1_R2:
        #fastqc
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/" + wildcards.sample_name_adt + "/QC_reads/" + name + "_fastqc.html"))
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/" + wildcards.sample_name_adt + "/QC_reads/" + name + "_fastqc.zip"))
    return files

rule multiqc_adt:
    input:
        qc_files = multiqc_inputs_adt
    output:
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/QC_reads/{sample_name_adt}_RAW.html")
    log:
        "logs/multiqc_adt/{sample_name_adt}.log"
    benchmark:
        "benchmark/multiqc_adt/{sample_name_adt}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_MULTIQC
    shell:
        """
        export TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX) && \
        multiqc \
            -n {wildcards.sample_name_adt}'_RAW' \
            -i {wildcards.sample_name_adt}' RAW FASTQ' \
            -f --no-megaqc-upload --no-data-dir \
            -o {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads/ {input.qc_files} &> {log}
        """

"""
This rule makes the ADT index by kb-python (wrapper of kallisto and bustools)
"""

rule kb_index_adt:
    input:
        features_file = FEATURES_ADT
    output:
        tr2g = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT_tr2g.txt"),
        index = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT.kidx"),
        fasta = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT.fa")
    log:
        "logs/kb_index_adt/kb_index_adt.log"
    benchmark:
        "benchmark/kb_index_adt/kb_index_adt.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_KB_PYTHON
    shell:
        """
        mkdir -p {ALIGN_OUTPUT_DIR_ADT}/INDEX_ADT && \
        kb ref \
            -i {output.index} \
            -f1 {output.fasta} \
            -g {output.tr2g} \
            --tmp $(mktemp -d {resources.tmpdir}/XXXXXX)"/INDEX_ADT" \
            --workflow kite \
            {input.features_file} &> {log}
        """

"""
This rule makes the alignment, sort, correct, and count UMI by kb-python (wrapper of kallisto and bustools)
The function alignment_inputs_adt allows to get all fastq input files for one specific sample (wildcards).
"""
def alignment_inputs_adt(wildcards):
    return sorted([elem for elem in ALIGN_SYMLINK_FILES_ADT if re.search(wildcards.sample_name_adt, elem)])

rule alignment_build_count_matrix_adt:
    input:
        fq_copied = alignment_inputs_adt,
        tr2g = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT_tr2g.txt"),
        index = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT.kidx"),
        whitelist = WHITELISTNAME
    output:
        mtx_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.mtx"),
        barcodes_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.barcodes.txt"),
        genes_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.genes.txt"),
        MandM = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/Materials_and_Methods.txt")
    params:
        kallistobus_path = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS")
    log:
        "logs/alignment_build_count_matrix_adt/{sample_name_adt}.log"
    benchmark:
        "benchmark/alignment_build_count_matrix_adt/{sample_name_adt}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 40960, 204800)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_KB_PYTHON
    shell:
        """
        res=$(({resources.mem_mb}-100))
        mkdir -p {params.kallistobus_path} && \
        kb count \
            -i {input.index} \
            -g {input.tr2g} \
            -x {SCTECH} \
            -w {input.whitelist} \
            -o {params.kallistobus_path}/{wildcards.sample_name_adt} \
            --tmp $(mktemp -d {resources.tmpdir}/XXXXXX)/{wildcards.sample_name_adt} \
            -t {threads} \
            -m $res"M" \
            --h5ad --overwrite --workflow kite \
            {KB_PYTHON_EXTRA_ADT} \
            {input.fq_copied} &> {log}

        FASTQC_V=$(grep fastqc {PIPELINE_FOLDER}/envs/conda/fastqc.yaml | cut -d= -f2)
        KBPYTHON_V=$(grep kb-python {PIPELINE_FOLDER}/envs/conda/kb-python.yaml | cut -d= -f2)
        if [[ {SCTECH} = '10xv3' ]];then
            CR="10X Chromium 3' scRNA-Seq v3 chemistry"
        elif [[ {SCTECH} = '10xv2' ]];then
            CR="10X Chromium 5' scRNA-Seq v2 chemistry"
        fi
        echo "Raw BCL-files were demultiplexed and converted to Fastq format using bcl2fastq (version 2.20.0.422 from Illumina).
Reads quality control was performed using fastqc (version $FASTQC_V).
A customized index, with the correspondance between synthetic DNA-transcripts and tagged protein names, was made with the kb-python (version $KBPYTHON_V). Reads were pseudo-mapped to that customized index. Index built, reads pseudo-mapping with parameter corresponding to the $CR, barcode correction using whitelist provided by the manufacturer (10X Genomics), and gene-based reads quantification, were performed by kb-python." > {output.MandM}
        """
