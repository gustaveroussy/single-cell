"""
##########################################################################
This rule add tcr information to expression gene analysis in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file and csv file from cellranger vdj.
"""
def grp_add_tcr_input(wildcards):
    rda_file = dic_GRP_ADD_TCR_INFO[wildcards.grp_add_tcr_output]['GRP_ADD_TCR_INPUT_RDA']
    csv_file = list(dict.fromkeys(dic_GRP_ADD_TCR_INFO[wildcards.grp_add_tcr_output]['GRP_ADD_TCR_INPUT_CSV_TCR'].split(",")))
    csv_file.insert(0, rda_file)
    return csv_file

"""
This function allows to determine the singularity binding parameters.
"""
def grp_add_tcr_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_GRP_ADD_TCR_INFO[wildcards.grp_add_tcr_output]['GRP_ADD_TCR_INPUT_RDA']) # output_folder too
    concat = " -B " + PIPELINE_FOLDER + "," + rda_folder
    for tcrfile in list(dict.fromkeys(dic_GRP_ADD_TCR_INFO[wildcards.grp_add_tcr_output]['GRP_ADD_TCR_INPUT_CSV_TCR'].split(","))):
        tcrfile = os.path.dirname(tcrfile)
        concat = concat + "," + tcrfile
    return concat

"""
This function allows to determine the tcr files folders for params.
"""
def grp_add_tcr_params_tcr_files(wildcards):
    return ",".join([ os.path.normpath(tcrfile) for tcrfile in list(dict.fromkeys(dic_GRP_ADD_TCR_INFO[wildcards.grp_add_tcr_output]['GRP_ADD_TCR_INPUT_CSV_TCR'].split(","))) ])

"""
This function allows to determine the output folder for params (os.path.dirname() not allowed in params slot).
"""
def grp_add_tcr_params_output_folder(wildcards):
    return os.path.dirname(wildcards.grp_add_tcr_output) + "/"

"""
This rule launches the R script to add adt information to expression gene analysis.
"""
rule grp_add_tcr_ge:
    input:
        grp_add_tcr_file = grp_add_tcr_input
    output:
        grp_add_tcr_rda_file = "{grp_add_tcr_output}" + "_TCR.rda"
    params:
        sing_bind = grp_add_tcr_params_sing,
        input_csv = grp_add_tcr_params_tcr_files,
        output_folder = grp_add_tcr_params_output_folder
    log:
        "logs/grp_add_tcr_ge{grp_add_tcr_output}.log"
    benchmark:
        "benchmark/grp_add_tcr_ge{grp_add_tcr_output}.tsv"
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: 5120 + attempt * 5120,
        time_min = lambda wildcards, attempt: min(attempt * 120, 200)
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        singularity exec --no-home -B $TMPDIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV_TCR_BCR} \
        Rscript {PIPELINE_FOLDER}/scripts/Int_Grp_pipeline_TCR.R \
        --input.rda {input[0]} \
        --output.dir {params.output_folder}/ \
        --vdj.input.files.tcr {params.input_csv} \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --pipeline.path {PIPELINE_FOLDER} &> {log}
        """
