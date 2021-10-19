"""
##########################################################################
These rules make the alignment of TCR and BCR in single-cell RNA-seq.
##########################################################################
"""
wildcard_constraints:
    sample_name_tcr_bcr=".+(_TCR|_BCR)"

"""
This rule makes a copy of fastq files with the good sample name.
"""
def rsync_rename_inputs_tcr_bcr(wildcards):
    res = "no result find!"
    to_test = os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + wildcards.lane_R_complement + ".fastq.gz")
    for i in range(0,len(ALIGN_SYMLINK_FILES_TCR_BCR),1):
        if os.path.normpath(ALIGN_SYMLINK_FILES_TCR_BCR[i]) == to_test :
            res = str(ALIGN_ORIG_FILES_TCR_BCR[i])
    return res

rule rsync_rename_fq_tcr_bcr:
    input:
        fq = rsync_rename_inputs_tcr_bcr
    output:
        fq_link = temp(os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}{lane_R_complement}.fastq.gz")),
        md5_validation = temp(os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/md5_validation_{sample_name_tcr_bcr}{lane_R_complement}.txt"))
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(attempt * 256, 2048)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    # run:
    #     sys.stderr.write("\t Create symbolic link: \n")
    #     sys.stderr.write("\t From :" + "\t" + str(input.fq) + "\n")
    #     sys.stderr.write("\t To :" + "\t" + str(output.fq_link) + "\n")
    #     os.symlink(str(input.fq), str(output.fq_link))
    shell:
        """
        echo "From: {input.fq}"
        echo "To: {output.fq_link}"
        rsync -avz -c --quiet {input.fq} {output.fq_link}
        md5_fq=$(md5sum {input.fq} | cut -d' ' -f 1)
        md5_symlink=$(md5sum {output.fq_link} | cut -d' ' -f 1)
        [[ $md5_fq == $md5_symlink ]] && echo "md5sum OK" && touch {output.md5_validation}
        """

"""
This rule makes the fastqc control-quality.
"""
rule fastqc_tcr_bcr:
    input:
        fq = os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}{lane_R_complement}.fastq.gz")
    output:
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/QC_reads/fastqc/{sample_name_tcr_bcr}{lane_R_complement}_fastqc.html"),
        zip_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/QC_reads/fastqc/{sample_name_tcr_bcr}{lane_R_complement}_fastqc.zip")
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
        mkdir -p {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads/fastqc && \
        fastqc --quiet -o {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads/fastqc -d $TMP_DIR -t {threads} {input} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """


"""
This rule makes the fastq-screen control-quality on R2 files.
"""
rule fastqscreen_tcr_bcr:
    input:
        R2_fq = os.path.normpath(ALIGN_INPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}{lane_R_complement}.fastq.gz")
    output:
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR +"/{sample_name_tcr_bcr}/QC_reads/fastqscreen/{sample_name_tcr_bcr}{lane_R_complement}_screen.html"),
        #png_file = os.path.join(OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/QC_reads/fastqscreen/{sample_name_tcr_bcr}{lane_R_complement}_screen.png"),
        txt_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR +"/{sample_name_tcr_bcr}/QC_reads/fastqscreen/{sample_name_tcr_bcr}{lane_R_complement}_screen.txt")
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: min(2048 + attempt * 2048, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 30, 200))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        "mkdir -p {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads/fastqscreen && fastq_screen --quiet --threads {threads} --force --outdir {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads/fastqscreen --subset 100000 --conf {FASTQSCREEN_INDEX} {input.R2_fq}"


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
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + "/QC_reads/fastqc/" + name) + "_fastqc.html")
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + "/QC_reads/fastqc/" + name) + "_fastqc.zip")
    for name in name_R2:
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + "/QC_reads/fastqscreen/" + name) + "_screen.html")
        files.append(os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/" + wildcards.sample_name_tcr_bcr + "/QC_reads/fastqscreen/" + name) + "_screen.txt")
    return files

rule multiqc_tcr_bcr:
    input:
        qc_files = multiqc_inputs_tcr_bcr
    output:
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/QC_reads/{sample_name_tcr_bcr}_RAW.html")
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
        multiqc -n {wildcards.sample_name_tcr_bcr}'_RAW' -i {wildcards.sample_name_tcr_bcr}' RAW FASTQ' -p -z -f -o {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads {input} && \
        rm -r {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads/fastqc {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads/fastqscreen {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads/{wildcards.sample_name_tcr_bcr}_RAW_data.zip {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/QC_reads/{wildcards.sample_name_tcr_bcr}_RAW_plots && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """

"""
This function allows to determine the singularity binding parameters.
"""
def alignment_annotations_tcr_bcr_params_sing(wildcards):
    input_folder = os.path.dirname(ALIGN_SYMLINK_FILES_TCR_BCR[0])
    output_folder = os.path.dirname(ALIGN_OUTPUT_DIR_TCR_BCR + "/wildcards.sample_name_tcr_bcr")
    ref_folder = CRINDEX_TCR_BCR
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + input_folder + ":" + os.path.normpath("/WORKDIR/" + input_folder) + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder) + " -B " + ref_folder + ":" + os.path.normpath("/WORKDIR/" + ref_folder)
    return concat

"""
This rule makes the alignment and annotation by cellranger.
The function alignment_annotations_inputs_tcr_bcr allows to get all fastq input files for one specific sample (wildcards).
"""
def alignment_annotations_inputs_tcr_bcr(wildcards):
    files=[]
    files=[elem for elem in ALIGN_SYMLINK_FILES_TCR_BCR if re.search(wildcards.sample_name_tcr_bcr, elem)]
    return sorted(files)

rule alignment_annotations_tcr_bcr:
    input:
        fq = alignment_annotations_inputs_tcr_bcr,
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/QC_reads/{sample_name_tcr_bcr}_RAW.html")
    output:
        csv_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/filtered_contig_annotations.csv"),
        html_file = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/web_summary.html"),
	    MandM = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/Materials_and_Methods.txt")
    params:
        sing_bind = alignment_annotations_tcr_bcr_params_sing,
        sample_folder = os.path.normpath(ALIGN_OUTPUT_DIR_TCR_BCR + "/{sample_name_tcr_bcr}")
    threads:
        3
    resources:
        mem_mb = (lambda wildcards, attempt: min(1024 + attempt * 4096, 30720)),
        time_min = (lambda wildcards, attempt: min(1*24*60 + attempt * 24*60, 7*24*60))
    conda:
        CONDA_ENV_QC_ALIGN_GE_ADT
    shell:
        """
        #source /mnt/beegfs/software/cellranger/3.1.0/cellranger-3.1.0/sourceme.bash
        rm -r {params.sample_folder}/{wildcards.sample_name_tcr_bcr}_CellRanger
        res=$(({resources.mem_mb}/1000))

        echo 'cd /WORKDIR/{params.sample_folder} && \
        /Softwares/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger vdj \
                 --id={wildcards.sample_name_tcr_bcr}_CellRanger \
                 --reference=/WORKDIR/{CRINDEX_TCR_BCR} \
                 --fastqs=/WORKDIR/{ALIGN_INPUT_DIR_TCR_BCR} \
                 --sample={wildcards.sample_name_tcr_bcr} \
                 --localmem=$res \
                 --localcores={threads}' | singularity exec --contain {params.sing_bind} {SINGULARITY_ENV_TCR_BCR} bash
        rm -r {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/{wildcards.sample_name_tcr_bcr}_CellRanger/SC_VDJ_ASSEMBLER_CS
        rm {ALIGN_OUTPUT_DIR_TCR_BCR}/{wildcards.sample_name_tcr_bcr}/{wildcards.sample_name_tcr_bcr}_CellRanger/_*
        FASTQC_V=$(conda list "fastqc" | grep "^fastqc " | sed -e "s/fastqc *//g" | sed -e "s/ .*//g")
        FASTQSCREEN_V=$(conda list "fastq-screen" | grep "^fastq-screen " | sed -e "s/fastq-screen *//g" | sed -e "s/ .*//g")
        #CELLRANGER_V=`/Softwares/cellranger-3.1.0/cellranger-cs/3.1.0/bin/cellranger vdj --version | grep "cellranger vdj (" | sed -e "s/cellranger vdj (//g" | sed -e "s/)//g"`
        CELLRANGER_V="3.1.0"
        echo "Raw BCL-files were demultiplexed and converted to Fastq format using bcl2fastq (version 2.20.0.422 from Illumina).
Reads quality control was performed using fastqc (version $FASTQC_V) and assignment to the expected genome species evaluated with fastq-screen (version $FASTQSCREEN_V).
CellRanger (version $CELLRANGER_V from 10X Genomics) was used to generate single-cell V(D)J sequences and annotations." > {output.MandM}
        """
