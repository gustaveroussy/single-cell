"""
##########################################################################
This rule add adt information to expression gene analysis in single-cell RNA-seq.
##########################################################################
"""

wildcard_constraints:
    sample_name_ge=".+_GE"

"""
This function allows to determine the input .rda ge file and kallisto adt folder.
"""
def add_adt_input(wildcards):
    ge_rda_file = dic_ADD_ADT_INFO[wildcards.add_adt_output]['ADD_ADT_INPUT_RDA_GE']
    kallisto_folder = dic_ADD_ADT_INFO[wildcards.add_adt_output]['ADD_ADT_INPUT_DIR_ADT']
    if "Alignment_countTable_ADT" in STEPS:
        MandM = os.path.normpath(kallisto_folder + "/Materials_and_Methods.txt")
        files=[MandM]
    else:
        files=[kallisto_folder]
    files.insert(0,ge_rda_file)
    return files

"""
This function allows to determine the singularity binding parameters.
"""
def add_adt_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_ADD_ADT_INFO[wildcards.add_adt_output]['ADD_ADT_INPUT_RDA_GE']) # output_folder too
    kallisto_folder = os.path.dirname(dic_ADD_ADT_INFO[wildcards.add_adt_output]['ADD_ADT_INPUT_DIR_ADT'])
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder) + " -B " + kallisto_folder + ":" + os.path.normpath("/WORKDIR/" + kallisto_folder)
    return concat

"""
This function allows to determine the input alignment folder for params section.
"""
def add_adt_params_input_folder(wildcards):
    return os.path.normpath("/WORKDIR/" + dic_ADD_ADT_INFO[wildcards.add_adt_output]['ADD_ADT_INPUT_DIR_ADT']) + "/"

"""
This function allows to determine the output folder for params (os.path.dirname() not allowed in params slot).
"""
def add_adt_params_output_folder(wildcards):
    return os.path.normpath("/WORKDIR/" + os.path.dirname(wildcards.add_adt_output)) + "/"

"""
This function allows to determine the sample.name.adt for params.
"""
def add_adt_params_sample_name_adt(wildcards):
    return dic_ADD_ADT_INFO[wildcards.add_adt_output]['ADD_ADT_SAMPLE_NAME_ADT']


"""
This rule launches the R script to add adt information to expression gene analysis.
"""
rule add_adt_ge:
    input:
        add_adt_file = add_adt_input
    output:
        add_adt_rda_file = "{add_adt_output}" + "_ADT.rda"
    params:
        sing_bind = add_adt_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        kallisto_folder = add_adt_params_input_folder,
        output_folder = add_adt_params_output_folder,
        sample_name_adt = add_adt_params_sample_name_adt
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
        {SINGULARITY_ENV} \
        Rscript {params.pipeline_folder}/scripts/pipeline_ADT.R \
        --sample.name.adt {params.sample_name_adt} \
        --input.rda.ge {params.input_rda} \
        --output.dir {params.output_folder} \
        --input.dir.adt {params.kallisto_folder} \
        --author.name {ADD_ADT_AUTHOR_NAME} \
        --author.mail {ADD_ADT_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --gene.names  {ADD_ADT_GENE_NAMES} \
        --ADT.min.cutoff  {ADD_ADT_MIN_CUTOFF} \
        --ADT.max.cutoff  {ADD_ADT_MAX_CUTOFF} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
