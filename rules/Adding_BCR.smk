"""
##########################################################################
This rule add tcr information to expression gene analysis in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file and csv file from cellranger vdj.
"""
def add_bcr_input(wildcards):
    rda_file = dic_ADD_BCR_INFO[wildcards.add_bcr_output]['ADD_BCR_INPUT_RDA_GE']
    csv_file = dic_ADD_BCR_INFO[wildcards.add_bcr_output]['ADD_BCR_INPUT_CSV_BCR']
    return [rda_file,csv_file]

"""
This function allows to determine the singularity binding parameters.
"""
def add_bcr_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_ADD_BCR_INFO[wildcards.add_bcr_output]['ADD_BCR_INPUT_RDA_GE']) # output_folder too
    csv_folder = os.path.dirname(dic_ADD_BCR_INFO[wildcards.add_bcr_output]['ADD_BCR_INPUT_CSV_BCR'])
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder) + " -B " + csv_folder + ":" + os.path.normpath("/WORKDIR/" + csv_folder)
    return concat

"""
This function allows to determine the output folder for params (os.path.dirname() not allowed in params slot).
"""
def add_bcr_params_output_folder(wildcards):
    return os.path.normpath("/WORKDIR/" + os.path.dirname(wildcards.add_bcr_output)) + "/"

"""
This rule launches the R script to add adt information to expression gene analysis.
"""
rule add_bcr_ge:
    input:
        add_bcr_file = add_bcr_input
    output:
        add_bcr_rda_file = "{add_bcr_output}" + "_BCR.rda"
    params:
        sing_bind = add_bcr_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        input_csv = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[1]),
        output_folder = add_bcr_params_output_folder
    #conda:
    #    CONDA_ENV_SING
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(3072 + attempt * 1024, 10240)),
        time_min = (lambda wildcards, attempt: min(attempt * 60, 200))
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home -B $TMP_DIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV_TCR_BCR} \
        Rscript {params.pipeline_folder}/scripts/pipeline_BCR.R \
        --input.rda {params.input_rda} \
        --output.dir {params.output_folder} \
        --vdj.input.file.bcr {params.input_csv} \
        --author.name {ADD_BCR_AUTHOR_NAME} \
        --author.mail {ADD_BCR_AUTHOR_MAIL} \
        --pipeline.path {params.pipeline_folder} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
