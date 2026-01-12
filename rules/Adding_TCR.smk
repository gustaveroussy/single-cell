"""
##########################################################################
This rule add tcr information to expression gene analysis in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file and csv file from cellranger vdj.
"""
def add_tcr_input(wildcards):
    rda_file = dic_ADD_TCR_INFO[wildcards.add_tcr_output]['ADD_TCR_INPUT_RDA_GE']
    csv_file = dic_ADD_TCR_INFO[wildcards.add_tcr_output]['ADD_TCR_INPUT_CSV_TCR']
    return [rda_file,csv_file]

"""
This function allows to determine the singularity binding parameters.
"""
def add_tcr_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_ADD_TCR_INFO[wildcards.add_tcr_output]['ADD_TCR_INPUT_RDA_GE']) # output_folder too
    csv_folder = os.path.dirname(dic_ADD_TCR_INFO[wildcards.add_tcr_output]['ADD_TCR_INPUT_CSV_TCR'])
    concat = " -B " + PIPELINE_FOLDER + "," + rda_folder + "," + csv_folder
    return concat

"""
This function allows to determine the output folder for params (os.path.dirname() not allowed in params slot).
"""
def add_tcr_params_output_folder(wildcards):
    return os.path.dirname(wildcards.add_tcr_output)

"""
This rule launches the R script to add adt information to expression gene analysis.
"""
rule add_tcr_ge:
    input:
        add_tcr_file = add_tcr_input
    output:
        add_tcr_rda_file = "{add_tcr_output}" + "_TCR.rda"
    log:
        "logs/add_tcr_ge{add_tcr_output}.log"
    benchmark:
        "benchmark/add_tcr_ge{add_tcr_output}.tsv"
    params:
        sing_bind = add_tcr_params_sing,
        output_folder = add_tcr_params_output_folder
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5120,
        time_min = lambda wildcards, attempt: min(attempt * 60, 200)
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        singularity exec --no-home -B $TMPDIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV_TCR_BCR} \
        Rscript {PIPELINE_FOLDER}/scripts/pipeline_TCR.R \
        --input.rda {input[0]} \
        --output.dir {params.output_folder}/ \
        --vdj.input.file.tcr {input[1]} \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --pipeline.path {PIPELINE_FOLDER} &> {log}
        """
