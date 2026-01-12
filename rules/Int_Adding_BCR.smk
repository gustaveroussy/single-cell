"""
##########################################################################
This rule add bcr information to expression gene analysis in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file and csv file from cellranger vdj.
"""
def int_add_bcr_input(wildcards):
    rda_file = dic_INT_ADD_BCR_INFO[wildcards.int_add_bcr_output]['INT_ADD_BCR_INPUT_RDA']
    csv_file = list(dict.fromkeys(dic_INT_ADD_BCR_INFO[wildcards.int_add_bcr_output]['INT_ADD_BCR_INPUT_CSV_BCR'].split(",")))
    csv_file.insert(0, rda_file)
    return csv_file

"""
This function allows to determine the singularity binding parameters.
"""
def int_add_bcr_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_INT_ADD_BCR_INFO[wildcards.int_add_bcr_output]['INT_ADD_BCR_INPUT_RDA']) # output_folder too
    concat = " -B " + PIPELINE_FOLDER + "," + rda_folder
    for bcrfile in list(dict.fromkeys(dic_INT_ADD_BCR_INFO[wildcards.int_add_bcr_output]['INT_ADD_BCR_INPUT_CSV_BCR'].split(","))):
        bcrfile = os.path.dirname(bcrfile)
        concat = concat + "," + bcrfile
    return concat

"""
This function allows to determine the bcr files folders for params.
"""
def int_add_bcr_params_bcr_files(wildcards):
    return ",".join([ os.path.normpath(bcrfile) for bcrfile in list(dict.fromkeys(dic_INT_ADD_BCR_INFO[wildcards.int_add_bcr_output]['INT_ADD_BCR_INPUT_CSV_BCR'].split(","))) ])

"""
This function allows to determine the output folder for params (os.path.dirname() not allowed in params slot).
"""
def int_add_bcr_params_output_folder(wildcards):
    return os.path.dirname(wildcards.int_add_bcr_output) + "/"

"""
This rule launches the R script to add adt information to expression gene analysis.
"""
rule int_add_bcr_ge:
    input:
        int_add_bcr_file = int_add_bcr_input
    output:
        int_add_bcr_rda_file = "{int_add_bcr_output}" + "_BCR.rda"
    params:
        sing_bind = int_add_bcr_params_sing,
        input_csv = int_add_bcr_params_bcr_files,
        output_folder = int_add_bcr_params_output_folder
    log:
        "logs/int_add_bcr_ge{int_add_bcr_output}.log"
    benchmark:
        "benchmark/int_add_bcr_ge{int_add_bcr_output}.tsv"
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
        Rscript {PIPELINE_FOLDER}/scripts/Int_Grp_pipeline_BCR.R \
        --input.rda {input[0]} \
        --output.dir {params.output_folder}/ \
        --vdj.input.files.bcr {params.input_csv} \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --pipeline.path {PIPELINE_FOLDER} &> {log}
        """
