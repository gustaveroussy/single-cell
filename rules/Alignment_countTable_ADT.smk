"""
##########################################################################
These rules make the alignment of cell surface protein in single-cell RNA-seq.
##########################################################################
"""
wildcard_constraints:
    sample_name_adt = ".+_ADT"

"""
This rule makes the symbolic links of fastq files with the good sample name.
"""
def symlink_rename_inputs_adt(wildcards):
    for i in range(0,len(ALIGN_SAMPLE_NAME_ADT ),1):
        if ALIGN_SAMPLE_NAME_ADT[i] == wildcards.sample_name_adt :
            return os.path.normpath(ALIGN_INPUT_DIR_ADT_RAW + "/" + ALIGN_SAMPLE_NAME_ADT_RAW[i] + str("{lane_R_complement}.fastq.gz"))

rule symlink_rename_fq_adt:
    input:
        fq = symlink_rename_inputs_adt
    output:
        fq_link = temp(os.path.normpath(ALIGN_INPUT_DIR_ADT + "/{sample_name_adt}{lane_R_complement}.fastq.gz"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 256, 2048)),
        time_min = (lambda wildcards, attempt: min(attempt * 5, 50))
    run:
        os.environ["TMPDIR"] = GLOBAL_TMP
        sys.stderr.write("\t Create symbolic link: \n")
        sys.stderr.write("\t From :" + "\t" + str(input.fq) + "\n")
        sys.stderr.write("\t To :" + "\t" + str(output.fq_link) + "\n")
        os.symlink(str(input.fq), str(output.fq_link))

"""
This rule makes the fastqc control-quality.
"""
rule fastqc_adt:
    input:
        fq = os.path.normpath(ALIGN_INPUT_DIR_ADT + "/{sample_name_adt}{lane_R_complement}.fastq.gz")
    output:
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/QC_reads/fastqc/{sample_name_adt}{lane_R_complement}_fastqc.html"),
        zip_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/QC_reads/fastqc/{sample_name_adt}{lane_R_complement}_fastqc.zip")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 512, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        mkdir -p {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads/fastqc && \
        fastqc --quiet -o {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads/fastqc -d $TMP_DIR -t {threads} {input} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
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
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/" + wildcards.sample_name_adt + "/QC_reads/fastqc/" + name + "_fastqc.html"))
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/" + wildcards.sample_name_adt + "/QC_reads/fastqc/" + name + "_fastqc.zip"))
    return files

rule multiqc_adt:
    input:
        qc_files2 = multiqc_inputs_adt
    output:
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/QC_reads/{sample_name_adt}_RAW.html")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 256, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        export TMPDIR=$TMP_DIR && \
        multiqc -n {wildcards.sample_name_adt}'_RAW' -i {wildcards.sample_name_adt}' RAW FASTQ' -p -z -f -o {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads {input} && \
        rm -r {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads/fastqc {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads/{wildcards.sample_name_adt}_RAW_data.zip {ALIGN_OUTPUT_DIR_ADT}/{wildcards.sample_name_adt}/QC_reads/{wildcards.sample_name_adt}_RAW_plots && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """


"""
This rule makes the alignment by kallisto.
The function alignment_inputs_adt allows to get all fastq input files for one specific sample (wildcards).
"""
def alignment_inputs_adt(wildcards):
    return sorted([elem for elem in ALIGN_SYMLINK_FILES_ADT if re.search(wildcards.sample_name_adt, elem)])

rule alignment_adt:
    input:
        fq_link = alignment_inputs_adt,
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/QC_reads/{sample_name_adt}_RAW.html")
    output:
        output_bus_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/output.bus"),
        transcripts_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/transcripts.txt"),
        matrix_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/matrix.ec"),
        run_info_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/run_info.json")
    params:
        kbusdir = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: ALIGN_MEM_ADT if (ALIGN_MEM_ADT is not None) else min(attempt * 256, 10240)),
        time_min = (lambda wildcards, attempt: ALIGN_TIME_ADT if (ALIGN_TIME_ADT is not None) else min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "mkdir -p {params.kbusdir} && kallisto bus -i {KINDEX_ADT} -o {params.kbusdir} -x {SCTECH} -t {threads} {input.fq_link}"

"""
This rule correct UMI from the sorted results of alignment, by bustools.
"""
rule correct_UMIs_adt:
    input:
        output_bus_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/output.bus")
    output:
        corrected_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}_corrected.bus")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: ALIGN_MEM_ADT if (ALIGN_MEM_ADT is not None) else min(attempt * 256, 10240)),
        time_min = (lambda wildcards, attempt: ALIGN_TIME_ADT if (ALIGN_TIME_ADT is not None) else min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "bustools correct -w {WHITELISTNAME} -o {output} {input} && rm {input}"

"""
This rule sort the results of alignment, by bustools.
"""
rule sort_file_adt:
    input:
        corrected_bus_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}_corrected.bus")
    output:
        sorted_bus_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}_sorted.bus")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: ALIGN_MEM_ADT if (ALIGN_MEM_ADT is not None) else min(12288 + attempt * 3072, 30720)),
        time_min = (lambda wildcards, attempt: ALIGN_TIME_ADT if (ALIGN_TIME_ADT is not None) else min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        res=$(({resources.mem_mb}-512)) && \
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        bustools sort -T $TMP_DIR -t {threads} -m $res"M" -o {output} {input} && \
        rm -r {input} $TMP_DIR || rm -r $TMP_DIR
        """

"""
This rule count UMI from the corrected sorted results of alignment, by bustools.
"""
rule build_count_matrix_adt:
    input:
        sorted_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}_sorted.bus"),
        transcripts_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/transcripts.txt"),
        matrix_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/matrix.ec")
    output:
        mtx_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.mtx"),
        barcodes_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.barcodes.txt"),
        genes_file = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.genes.txt"),
	    MandM = os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS/Materials_and_Methods.txt")
    params:
        os.path.normpath(ALIGN_OUTPUT_DIR_ADT + "/{sample_name_adt}/KALLISTOBUS")
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: ALIGN_MEM_ADT if (ALIGN_MEM_ADT is not None) else min(attempt * 256, 10240)),
        time_min = (lambda wildcards, attempt: ALIGN_TIME_ADT if (ALIGN_TIME_ADT is not None) else min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        """
        bustools count --genecounts -o {params}/{wildcards.sample_name_adt} -g {TR2GFILE_ADT} -e {input.matrix_file} -t {input.transcripts_file} {input.sorted_file} && rm {input}
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
A customized index, with the correspondance between synthetic DNA-transcripts and tagged protein names, was made with the kb-python (version $KBPYTHON_V) wrapper of kallisto. Reads were pseudo-mapped to the customized index with kallisto (version $KALLISTO_V) using its 'bus' subcommand and parameters corresponding to the $CR. Barcode correction using whitelist provided by the manufacturer (10X Genomics) and gene-based reads quantification was performed with BUStools (version $BUSTOOLS_V)." > {output.MandM}
        """
