"""
##########################################################################
This rule make the droplets control-quality of genes expression in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input alignment folder.
"""
def QC_input_folder(w_outputqc_droplets_dir_ge):
    for key, value in dic_SAMPLE_NAME_GE_INFO.items():
        if value["QC_OUTPUT_DIR"] == w_outputqc_droplets_dir_ge :
            return key
            #key is the input value in the yaml parameter file.

"""
This function allows to determine the input alignment folder or files.
"""
def QC_droplets_input_ge(wildcards):
    kallisto_folder = QC_input_folder(wildcards.outputqc_droplets_dir_ge)
    if "Alignment_countTable_GE" in STEPS:
        mtx_file = os.path.normpath(kallisto_folder + "/" + wildcards.sample_name_ge + ".mtx")
        barcodes_file = os.path.normpath(kallisto_folder + "/" + wildcards.sample_name_ge + ".barcodes.txt")
        genes_file = os.path.normpath(kallisto_folder + "/" + wildcards.sample_name_ge + ".genes.txt")
        files=[mtx_file, barcodes_file, genes_file]
    else:
        files=[kallisto_folder]
    return files

"""
This function allows to determine the input alignment folder for params section.
"""
def QC_params_input_folder(wildcards):

    input_folder = os.path.normpath(QC_input_folder(wildcards.outputqc_droplets_dir_ge)) + "/"
    return input_folder

"""
This function allows to determine the singularity binding parameters.
"""
def QC_params_sing(wildcards):
    kallisto_folder = QC_input_folder(wildcards.outputqc_droplets_dir_ge)
    output_folder = wildcards.outputqc_droplets_dir_ge + "/"
    concat = " -B " + PIPELINE_FOLDER + "," + kallisto_folder + "," + output_folder
    if QC_MT_FILE != "NULL": concat = concat + "," + QC_MT_FILE
    if QC_RB_FILE != "NULL": concat = concat + "," + QC_RB_FILE
    if QC_ST_FILE != "NULL": concat = concat + "," + QC_ST_FILE
    if QC_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(QC_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + "," + metadatafile
    return concat

"""
This rule launches R scipt to read count matrix and perform droplets control-quality.
"""
rule QC_droplets_ge:
    input:
        QC_droplets_input_ge
    output:
        kneeplot_saturation_file = os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_kneeplot_saturation.png") if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_kneeplot_saturation.png"),
        QC_hist_unfiltred_file =  os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_QChist.png") if str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" +  "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_QChist.png"),
        unfiltred_non_norm_rda = os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_QC_NON-NORMALIZED.rda") if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_QC_NON-NORMALIZED.rda")
    params:
        sing_bind = QC_params_sing,
        input_folder = QC_params_input_folder,
    log:
        "logs/QC_droplets_ge{outputqc_droplets_dir_ge}{sample_name_ge}.log"
    benchmark:
        "benchmark/QC_droplets_ge{outputqc_droplets_dir_ge}{sample_name_ge}.tsv"
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: min(3072 + attempt * 3072, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 90, 200))
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        singularity exec --no-home -B $TMPDIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV} \
        Rscript {PIPELINE_FOLDER}/scripts/pipeline_part1.R \
        --input.dir.ge {params.input_folder} \
        --output.dir.ge {wildcards.outputqc_droplets_dir_ge}/ \
        --sample.name.ge {wildcards.sample_name_ge} \
        --species {SPECIES} \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --nthreads {threads} \
        --pipeline.path {PIPELINE_FOLDER} \
        --emptydrops.fdr {QC_EMPTYDROPS_FDR} \
        --droplets.limit {QC_DROPLETS_LIMIT} \
        --emptydrops.retain {QC_EMPTYDROPS_RETAIN} \
        --pcmito.min {QC_PCMITO_MIN} \
        --pcmito.max {QC_PCMITO_MAX} \
        --pcribo.min {QC_PCRIBO_MIN} \
        --pcribo.max {QC_PC_RIBO_MAX} \
        --min.features {QC_MIN_FEATURES} \
        --min.counts {QC_MIN_COUNTS} \
        --min.cells {QC_MIN_CELLS} \
        --mt.genes.file {QC_MT_FILE} \
        --crb.genes.file {QC_RB_FILE} \
        --str.genes.file {QC_ST_FILE} \
        --metadata.file {QC_METADATA_FILE} &> {log}
        """
