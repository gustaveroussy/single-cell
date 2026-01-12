"""
##########################################################################
These rules make the alignment of genes expression in single-cell RNA-seq.
##########################################################################
"""

"""
This rule makes the symbolic links of fastq files with the good sample name.
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
        "logs/copy_rename_fq_ge/{sample_name_ge}{lane_R_complement}.log"
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
        "logs/fastqc_ge/{sample_name_ge}{lane_R_complement}.log"
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
        fastqc --outdir {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads --dir $TMPDIR --threads {threads} {input.fq_copied} &> {log}
        """

"""
This rule makes the fastq-screen control-quality on R2 files.
"""
rule fastqscreen_ge:
    input:
        fq_copied = os.path.normpath(ALIGN_INPUT_DIR_GE + "/{sample_name_ge}{lane_R_complement}.fastq.gz")
    output:
        html_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}{lane_R_complement}_screen.html")),
        txt_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}{lane_R_complement}_screen.txt")),
        png_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}{lane_R_complement}_screen.png"))
    log:
        "logs/fastqscreen_ge/{sample_name_ge}{lane_R_complement}.log"
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
        fastq_screen --quiet --threads {threads} --force --outdir {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads --subset 100000 --conf {FASTQSCREEN_INDEX} {input.fq_copied} &> {log}
        """

"""
This rule makes the multiqc from the fastqc and the fastq-screen results.
The function allows to get all QC input files for one specific sample (wildcards).
"""
def multiqc_inputs_ge(wildcards):
    name_R1_R2=[elem for elem in ALIGN_SYMLINK_FILES_NAME_GE if re.search(str( "^" + wildcards.sample_name_ge), elem)]
    name_R2=[elem for elem in name_R1_R2 if re.search("R2", elem)]
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
    return files

rule multiqc_ge:
    input:
        qc_files = multiqc_inputs_ge
    output:
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/QC_reads/{sample_name_ge}_RAW.html")
    log:
        "logs/multiqc_ge/{sample_name_ge}.log"
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
            -o {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads/ {input.qc_files} &> {log}
        """

"""
This rule makes the alignment by kallisto.
The function alignment_inputs_ge allows to get all fastq input files for one specific sample (wildcards).
"""
def alignment_inputs_ge(wildcards):
    return sorted([elem for elem in ALIGN_SYMLINK_FILES_GE if re.search(wildcards.sample_name_ge, elem)])

rule alignment_ge:
    input:
        fq_copied = alignment_inputs_ge,
        index = KINDEX_GE,
    output:
        output_bus_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/output.bus")),
        transcripts_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/transcripts.txt")),
        matrix_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/matrix.ec")),
        run_info_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/run_info.json")
    params:
        kbusdir = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS"),
        sctech = "10xv2" if (SCTECH.split("_")[0] == "10xv2") else "10xv3" if (SCTECH.split("_")[0] == "10xv3" or SCTECH.split("_")[0] == "10xv4") else "",
        stranding = "--fr-stranded" if (SCTECH.split("_")[1] == "3p") else "--rf-stranded" if (SCTECH.split("_")[1] == "5p") else ""
    log:
        "logs/alignment_inputs_ge/{sample_name_ge}.log"
    benchmark:
        "benchmark/alignment_inputs_ge/{sample_name_ge}.tsv"
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: min(5120 + attempt * 10240, 61440)),
        time_min = (lambda wildcards, attempt: min(attempt * 180, 200))
    conda:
        CONDA_ENV_KALLISTO
    shell:
        """
        kallisto bus {params.stranding} -i {input.index} -o {params.kbusdir} -x {params.sctech} -t {threads} {input.fq_copied} &> {log}
        """

"""
This rule correct UMI from the sorted results of alignment, by bustools.
"""
rule correct_UMIs_ge:
    input:
        output_bus_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/output.bus"),
        whitelist = WHITELISTNAME
    output:
        corrected_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/output_corrected.bus"))
    log:
        "logs/correct_UMIs_ge/{sample_name_ge}.log"
    benchmark:
        "benchmark/correct_UMIs_ge/{sample_name_ge}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_BUSTOOLS
    shell:
        """
        bustools correct -w {input.whitelist} -o {output} {input.output_bus_file} &> {log}
        """

"""
This rule sort the results of alignment, by bustools.
"""
rule sort_file_ge:
    input:
        corrected_bus_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/output_corrected.bus")
    output:
        sorted_bus_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/output_corrected_sorted.bus"))
    log:
        "logs/sort_file_ge/{sample_name_ge}.log"
    benchmark:
        "benchmark/sort_file_ge/{sample_name_ge}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(10240 + attempt * 5120, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_BUSTOOLS
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        res=$(({resources.mem_mb}-512)) && \
        bustools sort -T $TMPDIR -t {threads} -m $res"M" -o {output} {input} &> {log}
        """

"""
This rule count UMI from the corrected sorted results of alignment, by bustools.
"""
rule build_count_matrix_ge:
    input:
        sorted_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/output_corrected_sorted.bus"),
        transcripts_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/transcripts.txt"),
        matrix_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/matrix.ec"),
        tr2g = TR2GFILE_GE
    output:
        mtx_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.mtx"),
        barcodes_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.barcodes.txt"),
        genes_file = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.genes.txt"),
        MandM = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS/Materials_and_Methods.txt")
    params:
        kallistobus_path = os.path.normpath(ALIGN_OUTPUT_DIR_GE + "/{sample_name_ge}/KALLISTOBUS")
    log:
        "logs/build_count_matrix_ge/{sample_name_ge}.log"
    benchmark:
        "benchmark/build_count_matrix_ge/{sample_name_ge}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_BUSTOOLS
    shell:
        """
        bustools count --genecounts -o {params.kallistobus_path}/{wildcards.sample_name_ge} -g {input.tr2g} -e {input.matrix_file} -t {input.transcripts_file} {input.sorted_file} &> {log}

        if [[ {SCTECH} = "10xv4_3p" ]];then
            CR="10X Chromium 3' scRNA-Seq v4 chemistry"
        elif [[ {SCTECH} = "10xv3_5p" ]];then
            CR="10X Chromium 5' scRNA-Seq v3 chemistry"
        elif [[ {SCTECH} = "10xv3_3p" ]];then
            CR="10X Chromium 3' scRNA-Seq v3 chemistry"    
        elif [[ {SCTECH} = "10xv2_5p" ]];then
            CR="10X Chromium 5' scRNA-Seq v2 chemistry"
        fi

        FASTQC_V=$(grep fastqc {PIPELINE_FOLDER}/envs/conda/fastqc.yaml | cut -d= -f2)
        FASTQSCREEN_V=$(grep fastq-screen {PIPELINE_FOLDER}/envs/conda/fastq-screen.yaml | cut -d= -f2)
        KALLISTO_V=$(grep kallisto {PIPELINE_FOLDER}/envs/conda/kallisto.yaml | cut -d= -f2)
        BUSTOOLS_V=$(grep bustools {PIPELINE_FOLDER}/envs/conda/bustools.yaml | cut -d= -f2)
        KBPYTHON_V=$(grep kb-python {PIPELINE_FOLDER}/envs/conda/kb-python.yaml | cut -d= -f2)

        echo "Reads quality control was performed using fastqc (version $FASTQC_V) and assignment to the expected genome species evaluated with fastq-screen (version $FASTQSCREEN_V).
Reads were pseudo-mapped to the {REF_TXT_GE} with kb-python (version $KBPYTHON_V) using its 'bus' subcommand and parameters corresponding to the $CR. The index was made with kb-python (version $KBPYTHON_V). Barcode correction using whitelist provided by the manufacturer (10X Genomics), sorting, and gene-based reads quantification were performed with BUStools (version $BUSTOOLS_V)." > {output.MandM}
        """
