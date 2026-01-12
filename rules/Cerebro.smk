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
    concat = " -B " + PIPELINE_FOLDER + "," + rda_crb_folder
    if CEREBRO_GMT_FILE != "NULL":
        gmt_folder = os.path.dirname(CEREBRO_GMT_FILE)
        concat = concat + "," + gmt_folder
    if CEREBRO_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(CEREBRO_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + "," + metadatafile
    return concat




"""
This rule launches the R script to translate seurat V5 file into seurat V3 file.
"""
rule cerebro_seurat_conversion:
    input:
        rda_file = cerebro_input
    output:
        rda_file = temp("{cerebro_input_rda_no_extention}_SeuratV3.rda")
    params:
        sing_bind = cerebro_params_sing
    log:
        "logs/cerebro_seurat_conversion{cerebro_input_rda_no_extention}.log"
    benchmark:
        "benchmark/cerebro_seurat_conversion{cerebro_input_rda_no_extention}.tsv"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: 5120 * attempt,
        time_min = lambda wildcards, attempt: min(attempt * 60, 200)
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        singularity exec --contain -B $TMPDIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV} \
        Rscript {PIPELINE_FOLDER}/scripts/pipeline_Seurat_conversion_v5_to_v3.R \
        --input.rda.ge {input[0]} \
        --pipeline.path {PIPELINE_FOLDER} &> {log}
        """

"""
This rule launches the R script to translate seurat file into cerebro file.
"""
rule cerebro:
    input:
        cerebro_rda_file = "{cerebro_input_rda_no_extention}_SeuratV3.rda"
    output:
        cerebro_crb_file = expand("{{cerebro_input_rda_no_extention}}{cerebro_complement}", cerebro_complement = CEREBRO_COMPLEMENT_CRB)
    params:
        sing_bind = cerebro_params_sing
    log:
        "logs/cerebro{cerebro_input_rda_no_extention}.log"
    benchmark:
        "benchmark/cerebro{cerebro_input_rda_no_extention}.tsv"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: 5120 * attempt,
        time_min = lambda wildcards, attempt: min(attempt * 60, 360)
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        singularity exec --contain --home $TMPDIR:$HOME -B $TMPDIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV_CEREBRO} \
        Rscript {PIPELINE_FOLDER}/scripts/pipeline_CEREBRO.R \
        --input.rda.ge {input[0]} \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --nthreads {threads} \
        --pipeline.path {PIPELINE_FOLDER} \
        --version {CEREBRO_VERSION} \
        --groups {CEREBRO_GROUPS} \
        --remove.other.reductions {CEREBRO_REMOVE_OTHER_RED} \
        --remove.other.idents {CEREBRO_REMOVE_OTHER_IDENT} \
        --remove.mt.genes {CEREBRO_REMOVE_MT} \
        --remove.crb.genes {CEREBRO_REMOVE_CRB} \
        --remove.str.genes {CEREBRO_REMOVE_STR} \
        --only.pos.DE {CEREBRO_ONLY_POS_DE} \
        --remove.custom.DE {CEREBRO_REMOVE_CUSTOM_DE} \
        --add.surface.prot.info {CEREBRO_ADD_SURFACE_PROT_INFO} \
        --gmt.file {CEREBRO_GMT_FILE} \
        --metadata.file {CEREBRO_METADATA_FILE} &> {log}
        """
