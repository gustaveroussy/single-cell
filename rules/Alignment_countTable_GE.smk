"""
##########################################################################
These rules make the alignment of genes expression in single-cell RNA-seq.
##########################################################################
"""
wildcard_constraints:
    sample_name_ge_R=".+_GE",
    sample_name_ge=".+_GE"

"""
This rule makes the fastqc control-quality.
"""
rule fastqc_ge:
    input:
        fq = os.path.join(ALIGN_INPUT_DIR_GE,"{sample_name_ge_R}{lane_R_complement}.fastq.gz")
    output:
        html_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge_R}/QC_reads/fastqc/{sample_name_ge_R}{lane_R_complement}_fastqc.html"),
        zip_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge_R}/QC_reads/fastqc/{sample_name_ge_R}{lane_R_complement}_fastqc.zip")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "mkdir -p {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge_R}/QC_reads/fastqc && fastqc --quiet -o {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge_R}/QC_reads/fastqc -t {threads} {input}"


"""
This rule makes the fastq-screen control-quality on R2 files.
"""
rule fastqscreen_ge:
    input:
        R2_fq = os.path.join(ALIGN_INPUT_DIR_GE,"{sample_name_ge_R}{lane_R_complement}.fastq.gz")
    output:
        html_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge_R}/QC_reads/fastqscreen/{sample_name_ge_R}{lane_R_complement}_screen.html"),
        txt_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge_R}/QC_reads/fastqscreen/{sample_name_ge_R}{lane_R_complement}_screen.txt")
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: min(2048 + attempt * 2048, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "mkdir -p {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge_R}/QC_reads/fastqscreen && fastq_screen --quiet --threads {threads} --force --outdir {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge_R}/QC_reads/fastqscreen --subset 100000 --conf {FASTQSCREEN_INDEX} {input.R2_fq}"


"""
This rule makes the multiqc from the fastqc and the fastq-screen results.
The function allows to get all QC input files for one specific sample (wildcards).
"""
def multiqc_inputs_ge(wildcards):
    name_R1_R2=[elem for elem in ALL_FILES_GE if re.search(wildcards.sample_name_ge, elem)]
    name_R2=[elem for elem in name_R1_R2 if re.search("R2", elem)]
    files=[]
    for name in name_R1_R2:
        #fastqc
        files.append(os.path.join(ALIGN_OUTPUT_DIR_GE,wildcards.sample_name_ge,"QC_reads/fastqc",name) + "_fastqc.html")
        files.append(os.path.join(ALIGN_OUTPUT_DIR_GE,wildcards.sample_name_ge,"QC_reads/fastqc", name) + "_fastqc.zip")
    for name in name_R2:
        files.append(os.path.join(ALIGN_OUTPUT_DIR_GE,wildcards.sample_name_ge,"QC_reads/fastqscreen", name) + "_screen.html")
        #files.append(os.path.join(ALIGN_OUTPUT_DIR_GE,wildcards.sample_name_ge,"QC_reads/fastqscreen", name) + "_screen.png")
        files.append(os.path.join(ALIGN_OUTPUT_DIR_GE,wildcards.sample_name_ge,"QC_reads/fastqscreen", name) + "_screen.txt")
    return files

rule multiqc_ge:
    input:
        qc_files2 = multiqc_inputs_ge
    output:
        html_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/QC_reads/{sample_name_ge}_RAW.html")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "multiqc -n {wildcards.sample_name_ge}'_RAW' -i {wildcards.sample_name_ge}' RAW FASTQ' -p -z -f -o {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads {input} && rm -r {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads/fastqc {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads/fastqscreen {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads/{sample_name_ge}_RAW_data.zip {ALIGN_OUTPUT_DIR_GE}/{wildcards.sample_name_ge}/QC_reads/{sample_name_ge}_RAW_plots"


"""
This rule makes the alignment by kallisto.
The function alignment_inputs_ge allows to get all fastq input files for one specific sample (wildcards).
"""
def alignment_inputs_ge(wildcards):
    files=[]
    files=[elem for elem in PATH_ALL_FILES_GE_FQ_GZ if re.search(wildcards.sample_name_ge, elem)]
    return sorted(files)

rule alignment_ge:
    input:
        fq_link = alignment_inputs_ge,
        html_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/QC_reads/{sample_name_ge}_RAW.html")
    output:
        output_bus_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/output.bus"),
        transcripts_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/transcripts.txt"),
        matrix_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/matrix.ec"),
        run_info_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/run_info.json")
    params:
        kbusdir = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS")
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: min(6144 + attempt * 2048, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "mkdir -p {params.kbusdir} && kallisto bus -i {KINDEX_GE} -o {params.kbusdir} -x {SCTECH} -t {threads} {input.fq_link}"

"""
This rule correct UMI from the sorted results of alignment, by bustools.
"""
rule correct_UMIs_ge:
    input:
        output_bus_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/output.bus")
    output:
        corrected_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}_corrected.bus")
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    shell:
        "bustools correct -w {WHITELISTNAME} -o {output} {input} && rm {input}"

"""
This rule sort the results of alignment, by bustools.
"""
rule sort_file_ge:
    input:
        corrected_bus_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}_corrected.bus")
    output:
        sorted_bus_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}_sorted.bus")
    params:
        tmp_dir=os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/tmp")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(12288 + attempt * 2048, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "mkdir -p {params.tmp_dir} && bustools sort -T {params.tmp_dir}/tmp -t {threads} -m 12G -o {output} {input} && rm -r {input} {params.tmp_dir}"

"""
This rule count UMI from the corrected sorted results of alignment, by bustools.
"""
rule build_count_matrix_ge:
    input:
        sorted_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}_sorted.bus"),
        transcripts_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/transcripts.txt"),
        matrix_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/matrix.ec")
    output:
        mtx_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.mtx"),
        barcodes_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.barcodes.txt"),
        genes_file = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.genes.txt"),
	    MandM = os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/Materials_and_Methods.txt")
    params:
        os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        """
        bustools count --genecounts -o {params}/{wildcards.sample_name_ge} -g {TR2GFILE_GE} -e {input.matrix_file} -t {input.transcripts_file} {input.sorted_file} && rm {input}
        FASTQC_V=$(conda list "fastqc" | grep "^fastqc " | sed -e "s/fastqc *//g" | sed -e "s/ .*//g")
        FASTQSCREEN_V=$(conda list "fastq-screen" | grep "^fastq-screen " | sed -e "s/fastq-screen *//g" | sed -e "s/ .*//g")
        KALLISTO_V=$(conda list "kallisto" | grep "^kallisto " | sed -e "s/kallisto *//g" | sed -e "s/ .*//g")
        KBPYTHON_V=$(conda list "kb-python" | grep "^kb-python " | sed -e "s/kb-python *//g" | sed -e "s/ .*//g")
        BUSTOOLS_V=$(conda list "bustools" | grep "^bustools " | sed -e "s/bustools *//g" | sed -e "s/ .*//g")
        if [[ {SCTECH} = '10xv3' ]];then
            CR="10X Chromium 3' scRNA-Seq v3 chemistry"
        elif [[ {SCTECH} = '10xv2' ]];then
            CR="10X Chromium 5' scRNA-Seq v2 chemistry"
        fi
        echo "Raw BCL-files were demultiplexed and converted to Fastq format using bcl2fastq (version 2.20.0.422 from Illumina).
Reads quality control was performed using fastqc (version $FASTQC_V) and assignment to the expected genome species evaluated with fastq-screen (version $FASTQSCREEN_V).
Reads were pseudo-mapped to the {REF_TXT_GE} with kallisto (version $KALLISTO_V) using its 'bus' subcommand and parameters corresponding to the $CR. The index was made with the kb-python (version $KBPYTHON_V) wrapper of kallisto. Barcode correction using whitelist provided by the manufacturer (10X Genomics) and gene-based reads quantification was performed with BUStools (version $BUSTOOLS_V)." > {output.MandM}
        """
