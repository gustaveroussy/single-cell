"""
##########################################################################
This rule add tcr information to expression gene analysis in single-cell RNA-seq.
##########################################################################
"""
wildcard_constraints:
    grp_add_tcr_output = "|".join(GRP_ADD_TCR_OUTPUT)

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
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder)
    for tcrfile in list(dict.fromkeys(dic_GRP_ADD_TCR_INFO[wildcards.grp_add_tcr_output]['GRP_ADD_TCR_INPUT_CSV_TCR'].split(","))):
        tcrfile = os.path.dirname(tcrfile)
        concat = concat + " -B " + tcrfile + ":" + os.path.normpath("/WORKDIR/" + tcrfile)
    return concat

"""
This function allows to determine the tcr files folders for params.
"""
def grp_add_tcr_params_tcr_files(wildcards):
    return ",".join([ os.path.normpath("/WORKDIR/" + tcrfile) for tcrfile in list(dict.fromkeys(dic_GRP_ADD_TCR_INFO[wildcards.grp_add_tcr_output]['GRP_ADD_TCR_INPUT_CSV_TCR'].split(","))) ])

"""
This function allows to determine the output folder for params (os.path.dirname() not allowed in params slot).
"""
def grp_add_tcr_params_output_folder(wildcards):
    return os.path.normpath("/WORKDIR/" + os.path.dirname(wildcards.grp_add_tcr_output)) + "/"

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
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        input_csv = grp_add_tcr_params_tcr_files,
        output_folder = grp_add_tcr_params_output_folder
    threads:
        1
    resources:
        mem_mb = lambda wildcards, attempt: min(5120 + attempt * 3072, 20480),
        time_min = lambda wildcards, attempt: min(attempt * 120, 200)
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home -B $TMP_DIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV_TCR_BCR} \
        Rscript {params.pipeline_folder}/scripts/Int_Grp_pipeline_TCR.R \
        --input.rda {params.input_rda} \
        --output.dir {params.output_folder} \
        --vdj.input.files.tcr {params.input_csv} \
        --author.name {GRP_ADD_TCR_AUTHOR_NAME} \
        --author.mail {GRP_ADD_TCR_AUTHOR_MAIL} \
        --pipeline.path {params.pipeline_folder} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
