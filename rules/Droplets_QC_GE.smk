"""
##########################################################################
This rule make the droplets control-quality of genes expression in single-cell RNA-seq.
##########################################################################
"""

wildcard_constraints:
    sample_name_ge=".+_GE"

"""
This function allows to determine the input alignment folder/files.
"""
def QC_droplets_input_ge(wildcards):
    kallisto_folder = dic_SAMPLE_NAME_GE_INFO[wildcards.sample_name_ge]['QC_INPUT_DIR']
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
    input_folder = os.path.normpath("/WORKDIR/" + dic_SAMPLE_NAME_GE_INFO[wildcards.sample_name_ge]['QC_INPUT_DIR']) + "/"
    return input_folder

"""
This function allows to determine the singularity binding parameters.
"""
def QC_params_sing(wildcards):
    kallisto_folder = dic_SAMPLE_NAME_GE_INFO[wildcards.sample_name_ge]['QC_INPUT_DIR']
    output_folder = wildcards.outputqc_droplets_dir_ge + "/"
    concat = " -B " + PIPELINE_FOLDER + ":/WORKDIR/" + PIPELINE_FOLDER + " -B " + kallisto_folder + ":" + os.path.normpath("/WORKDIR/" + kallisto_folder) + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if QC_MT_FILE != "NULL": concat = concat + " -B " + QC_MT_FILE + ":" + os.path.normpath("/WORKDIR/" + QC_MT_FILE)
    if QC_RB_FILE != "NULL": concat = concat + " -B " + QC_RB_FILE + ":" + os.path.normpath("/WORKDIR/" + QC_RB_FILE)
    if QC_ST_FILE != "NULL": concat = concat + " -B " + QC_ST_FILE + ":" + os.path.normpath("/WORKDIR/" + QC_ST_FILE)
    if QC_TRANSLATION_FILE != "NULL": concat = concat + " -B " + QC_TRANSLATION_FILE + ":" + os.path.normpath("/WORKDIR/" + QC_TRANSLATION_FILE)
    if QC_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(QC_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches R scipt to read count matrix and perform droplets control-quality.
"""
rule QC_droplets_ge:
    input:
        QC_droplets_input_ge
    output:
        kneeplot_file = os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_kneeplot.png") if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_kneeplot.png"),
        saturation_file = os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_saturation_plot.png") if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_saturation_plot.png"),
        QC_hist_unfiltred_file =  os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_QChist.png") if str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" +  "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_QChist.png"),
        unfiltred_non_norm_rda = os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets/" + "{sample_name_ge}_QC_NON-NORMALIZED.rda") if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else os.path.normpath("{outputqc_droplets_dir_ge}" + "/QC_droplets_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_QC_NON-NORMALIZED.rda")
    params:
        sing_bind = QC_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        # input_folder = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]) + "/",
        input_folder = QC_params_input_folder,
        output_folder = os.path.normpath("/WORKDIR/" + "{outputqc_droplets_dir_ge}") + "/",
        SING_QC_MT_FILE = os.path.normpath("/WORKDIR/" + QC_MT_FILE) if QC_MT_FILE != "NULL" else "NULL",
        SING_QC_RB_FILE = os.path.normpath("/WORKDIR/" + QC_RB_FILE) if QC_RB_FILE != "NULL" else "NULL",
        SING_QC_ST_FILE = os.path.normpath("/WORKDIR/" + QC_ST_FILE) if QC_ST_FILE != "NULL" else "NULL",
        SING_QC_TRANSLATION_FILE = os.path.normpath("/WORKDIR", QC_TRANSLATION_FILE) if QC_TRANSLATION_FILE != "NULL" else "NULL",
        SING_QC_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in QC_METADATA_FILE.split(',')]) if QC_METADATA_FILE != "NULL" else "NULL"
    threads:
        2
    resources:
        mem_mb = (lambda wildcards, attempt: min(3072 + attempt * 3072, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 90, 200))
    shell:
        """
        singularity exec --contain {params.sing_bind} \
        {SINGULARITY_ENV} \
        Rscript {params.pipeline_folder}/scripts/pipeline_part1.R \
        --input.dir.ge {params.input_folder} \
        --output.dir.ge {params.output_folder} \
        --sample.name.ge {wildcards.sample_name_ge} \
        --species {QC_SPECIES} \
        --author.name {QC_AUTHOR_NAME} \
        --author.mail {QC_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --emptydrops.fdr {QC_EMPTYDROPS_FDR} \
        --droplets.limit {QC_DROPLETS_LIMIT} \
        --emptydrops.retain {QC_EMPTYDROPS_RETAIN} \
        --translation {QC_TRANSLATION_BOOL} \
        --pcmito.min {QC_PCMITO_MIN} \
        --pcmito.max {QC_PCMITO_MAX} \
        --pcribo.min {QC_PCRIBO_MIN} \
        --pcribo.max {QC_PC_RIBO_MAX} \
        --min.features {QC_MIN_FEATURES} \
        --min.counts {QC_MIN_COUNTS} \
        --min.cells {QC_MIN_CELLS} \
        --mt.genes.file {params.SING_QC_MT_FILE} \
        --crb.genes.file {params.SING_QC_RB_FILE} \
        --str.genes.file {params.SING_QC_ST_FILE} \
        --translation.file {params.SING_QC_TRANSLATION_FILE} \
        --metadata.file {QC_METADATA_FILE} \
        --metadata.file {params.SING_QC_METADATA_FILE}
        """
