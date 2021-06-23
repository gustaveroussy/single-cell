"""
##########################################################################
This rule add adt information to expression gene analysis in grouped single-cell RNA-seq.
##########################################################################
"""

wildcard_constraints:
    grp_add_adt_output = "|".join(GRP_ADD_ADT_OUTPUT)

"""
This function allows to determine the input .rda ge file and kallisto adt folder.
"""
def grp_add_adt_input(wildcards):
    ge_rda_file = dic_GRP_ADD_ADT_INFO[wildcards.grp_add_adt_output]['GRP_ADD_ADT_INPUT_RDA']
    kallisto_folder = list(dict.fromkeys(dic_GRP_ADD_ADT_INFO[wildcards.grp_add_adt_output]['GRP_ADD_ADT_INPUT_DIR_ADT'].split(",")))
    kallisto_folder.insert(0,ge_rda_file)
    return kallisto_folder

"""
This function allows to determine the singularity binding parameters.
"""
def grp_add_adt_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_GRP_ADD_ADT_INFO[wildcards.grp_add_adt_output]['GRP_ADD_ADT_INPUT_RDA']) # output_folder too
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder)
    for kallisto_folder in list(dict.fromkeys(dic_GRP_ADD_ADT_INFO[wildcards.grp_add_adt_output]['GRP_ADD_ADT_INPUT_DIR_ADT'].split(","))):
        kallisto_folder = os.path.dirname(kallisto_folder)
        concat = concat + " -B " + kallisto_folder + ":" + os.path.normpath("/WORKDIR/" + kallisto_folder)
    return concat

"""
This function allows to determine the input alignment folder for params section.
"""
def grp_add_adt_params_input_folder(wildcards):
    return ",".join([ os.path.normpath("/WORKDIR/" + kallisto_folder + "/") for kallisto_folder in list(dict.fromkeys(dic_GRP_ADD_ADT_INFO[wildcards.grp_add_adt_output]['GRP_ADD_ADT_INPUT_DIR_ADT'].split(","))) ])

"""
This function allows to determine the output folder for params (os.path.dirname() not allowed in params slot).
"""
def grp_add_adt_params_output_folder(wildcards):
    return os.path.normpath("/WORKDIR/" + os.path.dirname(wildcards.grp_add_adt_output)) + "/"

"""
This function allows to determine the sample.name.adt for params.
"""
def grp_add_adt_params_sample_name_adt(wildcards):
    return dic_GRP_ADD_ADT_INFO[wildcards.grp_add_adt_output]['GRP_ADD_ADT_SAMPLE_NAME_ADT']


"""
This rule launches the R script to add adt information to expression gene analysis.
"""
rule grp_add_adt_ge:
    input:
        grp_add_adt_file = grp_add_adt_input
    output:
        grp_add_adt_rda_file = "{grp_add_adt_output}" + "_ADT.rda"
    params:
        sing_bind = grp_add_adt_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        kallisto_folder = grp_add_adt_params_input_folder,
        output_folder = grp_add_adt_params_output_folder,
        sample_name_adt = grp_add_adt_params_sample_name_adt
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
        {SINGULARITY_ENV} \
        Rscript {params.pipeline_folder}/scripts/Int_Grp_pipeline_ADT.R \
        --samples.name.adt {params.sample_name_adt} \
        --input.rda {params.input_rda} \
        --output.dir {params.output_folder} \
        --input.dirs.adt {params.kallisto_folder} \
        --author.name {GRP_ADD_ADT_AUTHOR_NAME} \
        --author.mail {GRP_ADD_ADT_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --gene.names  {GRP_ADD_ADT_GENE_NAMES} \
        --ADT.min.cutoff {GRP_ADD_ADT_MIN_CUTOFF} \
        --ADT.max.cutoff {GRP_ADD_ADT_MAX_CUTOFF} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
