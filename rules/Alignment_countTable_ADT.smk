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
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        fastqc --outdir {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads --dir $TMPDIR --threads {threads} {input.fq_copied} &> {log}
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
        export TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        multiqc \
            -n {wildcards.sample_name_adt}'_RAW' \
            -i {wildcards.sample_name_adt}' RAW FASTQ' \
            -f --no-megaqc-upload --no-data-dir \
            -o {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads/ {input.qc_files} &> {log}
        """

"""
This rule makes the ADT index by kallisto
"""

rule kallisto_index_adt:
    input:
        features_fa_file = FEATURES_FA_ADT
    output:
        tr2g = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT_tr2g.txt"),
        index = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT.kidx"),
        #fasta = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT.fa")
    log:
        "logs/kallisto_index_adt/kallisto_index_adt.log"
    benchmark:
        "benchmark/kallisto_index_adt/kallisto_index_adt.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_KALLISTO
    shell:
        """
        #index (kallisto remove temporary folder by itself)
        kallisto index \
            -k 15 \
            -i {output.index} \
            --tmp $(mktemp -d {resources.tmpdir}/XXXXXX) \
            {input.features_fa_file} &>> {log}
        #tr2g
        grep ">" {input.features_fa_file} | sed 's/>//' | awk '{{print $1 "\t" $1}}' >> {output.tr2g}
        """


"""
This rule makes the alignment by kallisto.
The function alignment_inputs_adt allows to get all fastq input files for one specific sample (wildcards).
"""
def alignment_inputs_adt(wildcards):
    return sorted([elem for elem in ALIGN_SYMLINK_FILES_ADT if re.search(wildcards.sample_name_adt, elem)])

rule alignment_adt:
    input:
        fq_copied = alignment_inputs_adt,
        index = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT.kidx")
    output:
        output_bus_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/output.bus")),
        transcripts_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/transcripts.txt")),
        matrix_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/matrix.ec")),
        run_info_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/run_info.json")
    params:
        kbusdir = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS"),
        sctech = "10xv2" if (SCTECH.split("_")[0] == "10xv2") else "10xv3" if (SCTECH.split("_")[0] == "10xv3" or SCTECH.split("_")[0] == "10xv4") else ""
    log:
        "logs/alignment_inputs_adt/{sample_name_adt}.log"
    benchmark:
        "benchmark/alignment_inputs_adt/{sample_name_adt}.tsv"
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: min(6144 + attempt * 2048, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_KALLISTO
    shell:
        """
        kallisto bus --unstranded -i {input.index} -o {params.kbusdir} -x {params.sctech} -t {threads} {input.fq_copied} &> {log}
        """

"""
This rule correct UMI from the sorted results of alignment, by bustools.
"""
rule correct_UMIs_adt:
    input:
        output_bus_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/output.bus"),
        whitelist = WHITELISTNAME
    output:
        corrected_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/output_corrected.bus"))
    log:
        "logs/correct_UMIs_adt/{sample_name_adt}.log"
    benchmark:
        "benchmark/correct_UMIs_adt/{sample_name_adt}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 256, 10240)),
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
rule sort_file_adt:
    input:
        corrected_bus_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/output_corrected.bus")
    output:
        sorted_bus_file = temp(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/output_corrected_sorted.bus"))
    log:
        "logs/sort_file_adt/{sample_name_adt}.log"
    benchmark:
        "benchmark/sort_file_adt/{sample_name_adt}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(12288 + attempt * 2048, 20480)),
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
rule build_count_matrix_adt:
    input:
        sorted_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/output_corrected_sorted.bus"),
        transcripts_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/transcripts.txt"),
        matrix_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/matrix.ec"),
        tr2g = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/INDEX_ADT/ADT_tr2g.txt")
    output:
        mtx_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.mtx"),
        barcodes_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.barcodes.txt"),
        genes_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.genes.txt"),
        MandM = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/Materials_and_Methods.txt")
    params:
        kallistobus_path = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS")
    log:
        "logs/build_count_matrix_adt/{sample_name_adt}.log"
    benchmark:
        "benchmark/build_count_matrix_adt/{sample_name_adt}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_BUSTOOLS
    shell:
        """
        bustools count --genecounts -o {params.kallistobus_path}/{wildcards.sample_name_adt} -g {input.tr2g} -e {input.matrix_file} -t {input.transcripts_file} {input.sorted_file} &> {log}

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

        echo "Reads quality control was performed using fastqc (version $FASTQC_V).
Reads were pseudo-mapped to the antibodies' tags with kallisto (version $BUSTOOLS_V) using its 'bus' subcommand and parameters corresponding to the $CR, and the unstranded parameter. The index was made with kallisto (version $BUSTOOLS_V) with the k-mer parameter set to 15. Barcode correction using whitelist provided by the manufacturer (10X Genomics), sorting, and gene-based reads quantification were performed with BUStools (version $BUSTOOLS_V)." > {output.MandM}
        """
