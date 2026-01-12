"""
##########################################################################
This rule add adt information to expression gene analysis in grouped single-cell RNA-seq.
##########################################################################
"""

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
    concat = " -B " + PIPELINE_FOLDER + "," + rda_folder
    for kallisto_folder in list(dict.fromkeys(dic_GRP_ADD_ADT_INFO[wildcards.grp_add_adt_output]['GRP_ADD_ADT_INPUT_DIR_ADT'].split(","))):
        kallisto_folder = os.path.dirname(kallisto_folder)
        concat = concat + "," + kallisto_folder
    return concat

"""
This function allows to determine the input alignment folder for params section.
"""
def grp_add_adt_params_input_folder(wildcards):
    return ",".join([ os.path.normpath(kallisto_folder + "/") for kallisto_folder in list(dict.fromkeys(dic_GRP_ADD_ADT_INFO[wildcards.grp_add_adt_output]['GRP_ADD_ADT_INPUT_DIR_ADT'].split(","))) ])

"""
This function allows to determine the output folder for params (os.path.dirname() not allowed in params slot).
"""
def grp_add_adt_params_output_folder(wildcards):
    return os.path.dirname(wildcards.grp_add_adt_output)

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
        kallisto_folder = grp_add_adt_params_input_folder,
        output_folder = grp_add_adt_params_output_folder,
        sample_name_adt = grp_add_adt_params_sample_name_adt
    log:
        "logs/grp_add_adt_ge{grp_add_adt_output}.log"
    benchmark:
        "benchmark/grp_add_adt_ge{grp_add_adt_output}.tsv"
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
        {SINGULARITY_ENV} \
        Rscript {PIPELINE_FOLDER}/scripts/Int_Grp_pipeline_ADT.R \
        --samples.name.adt {params.sample_name_adt} \
        --input.rda {input[0]} \
        --output.dir {params.output_folder}/ \
        --input.dirs.adt {params.kallisto_folder} \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --nthreads {threads} \
        --pipeline.path {PIPELINE_FOLDER} \
        --gene.names {GRP_ADD_ADT_GENE_NAMES} &> {log}
        """
