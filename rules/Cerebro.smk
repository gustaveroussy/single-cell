"""
##########################################################################
This rule make the translation of seurat file into cerebro file in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file.
"""
def cerebro_input(wildcards):
    return wildcards.cerebro_input_rda_no_extention + ".rda"

"""
This function allows to determine the singularity binding parameters.
"""
def cerebro_params_sing(wildcards):
    rda_crb_folder = os.path.dirname(wildcards.cerebro_input_rda_no_extention)
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_crb_folder + ":" + os.path.normpath("/WORKDIR/" + rda_crb_folder)
    if CEREBRO_GMT_FILE != "NULL":
        gmt_folder = os.path.dirname(CEREBRO_GMT_FILE)
        concat = concat + " -B " + gmt_folder + ":" + os.path.normpath("/WORKDIR/" + gmt_folder)
    if CEREBRO_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(CEREBRO_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches the R script to translate seurat file into cerebro file.
"""
rule cerebro:
    input:
        cerebro_rda_file = cerebro_input
    output:
        cerebro_crb_file = expand("{{cerebro_input_rda_no_extention}}{cerebro_complement}", cerebro_complement = CEREBRO_COMPLEMENT_CRB)
    params:
        sing_bind = cerebro_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        SING_CEREBRO_GMT_FILE = os.path.normpath("/WORKDIR/" + CEREBRO_GMT_FILE) if CEREBRO_GMT_FILE != "NULL" else "NULL",
        SING_CEREBRO_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in CEREBRO_METADATA_FILE.split(',')]) if CEREBRO_METADATA_FILE != "NULL" else "NULL"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(5120 * attempt , 102400)),
        time_min = (lambda wildcards, attempt: min(attempt * 60, 200))
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --contain --home $TMP_DIR:$HOME -B $TMP_DIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV_CEREBRO} \
        Rscript {params.pipeline_folder}/scripts/pipeline_CEREBRO.R \
        --input.rda.ge {params.input_rda} \
        --author.name {CEREBRO_AUTHOR_NAME} \
        --author.mail {CEREBRO_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --version {CEREBRO_VERSION} \
        --groups {CEREBRO_GROUPS} \
        --remove.other.reductions {CEREBRO_REMOVE_OTHER_RED} \
        --remove.other.idents {CEREBRO_REMOVE_OTHER_IDENT} \
        --remove.mt.genes {CEREBRO_REMOVE_MT} \
        --remove.crb.genes {CEREBRO_REMOVE_CRB} \
        --remove.str.genes {CEREBRO_REMOVE_STR} \
        --only.pos.DE {CEREBRO_ONLY_POS_DE} \
        --remove.custom.DE {CEREBRO_REMOVE_CUSTOM_DE} \
        --gmt.file {params.SING_CEREBRO_GMT_FILE} \
        --metadata.file {params.SING_CEREBRO_METADATA_FILE}
        """
