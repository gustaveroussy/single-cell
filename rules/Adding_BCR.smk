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
    concat = " -B " + PIPELINE_FOLDER + "," + rda_folder + "," + csv_folder
    return concat

"""
This function allows to determine the output folder for params (os.path.dirname() not allowed in params slot).
"""
def add_bcr_params_output_folder(wildcards):
    return os.path.dirname(wildcards.add_bcr_output)

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
        output_folder = add_bcr_params_output_folder
    log:
        "logs/add_bcr_ge{add_bcr_output}.log"
    benchmark:
        "benchmark/add_bcr_ge{add_bcr_output}.tsv"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5120,
        time_min = lambda wildcards, attempt: min(attempt * 60, 200)
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        singularity exec --no-home -B  $TMPDIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV_TCR_BCR} \
        Rscript {PIPELINE_FOLDER}/scripts/pipeline_BCR.R \
        --input.rda {input[0]} \
        --output.dir {params.output_folder}/ \
        --vdj.input.file.bcr {input[1]} \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --pipeline.path {PIPELINE_FOLDER} &> {log}
        """
