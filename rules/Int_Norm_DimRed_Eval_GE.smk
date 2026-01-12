"""
##########################################################################
This rule make the integration, normalization, dimensions reduction and evaluation in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file and format it for R script.
"""
def int_norm_dimred_input_ge_R(wildcards):
    for key, value in dic_INT_NDRE_INFO.items():
        if value["INT_NDRE_OUTPUT_DIR"] == (wildcards.output_int_norm_dimred_dir_ge + "/GROUPED_ANALYSIS/INTEGRATED/" + wildcards.name_int) :
            return ','.join([os.path.normpath(x) for x in key.split(',')])
            #key is the input value in the yaml parameter file.

"""
This function allows to determine the input .rda file.
"""
def int_norm_dimred_input_ge_snakemake(wildcards):
    return (int_norm_dimred_input_ge_R(wildcards)).split(',')

"""
This function allows to determine the singularity binding parameters.
"""
def int_norm_dimred_params_sing(wildcards):
    for key, value in dic_INT_NDRE_INFO.items():
        if value["INT_NDRE_OUTPUT_DIR"] == (wildcards.output_int_norm_dimred_dir_ge + "/GROUPED_ANALYSIS/INTEGRATED/" + wildcards.name_int) :
            sing_rda_folder = ','.join([os.path.dirname(x) for x in key.split(',')])
            #key is the input value in the yaml parameter file.
    output_folder = wildcards.output_int_norm_dimred_dir_ge
    concat = " -B " + PIPELINE_FOLDER + "," + sing_rda_folder + "," + output_folder
    if INT_NDRE_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(INT_NDRE_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + "," + metadatafile
    return concat

"""
These functions allows to find and format the biases for the R script from wilcards.
"""
def int_norm_dimred_biases_norm(wildcards):
    #recherche si c'est dans int_ndre_norm_vtr ou s'il n'y en a pas.
    for prefix in ("SCTransform_", "LogNormalize_"):
        if wildcards.int_ndre_norm_vtr.startswith(prefix):
            biases = wildcards.int_ndre_norm_vtr.replace(prefix, "")
            break
    else:
        return "NULL"
    #rechercher parmi tous les ensembles de correcction de biais
    for x in INT_NDRE_VTR_BIASES_NORM:
        #le group de biais dont tous les élements sont dans biases
        to_search = "_".join(sorted(list(dict.fromkeys(x.split(",")))))
        if to_search == biases:
            return x

def int_norm_dimred_biases_dimred(wildcards):
    #recherche si c'est dans int_ndre_dimred_vtr ou s'il n'y en a pas.
    for prefix in ("scbfa_", "bpca_", "pca_", "ica_", "mds_"):
        if wildcards.int_ndre_dimred_vtr.startswith(prefix):
            biases = wildcards.int_ndre_dimred_vtr.replace(prefix, "")
            for i in INT_NDRE_VTR_BATCH.split(","):
                biases = biases.replace(i, "").replace("__","_")
                if biases == "":
                    return "NULL"
            break
    else:
        return "NULL"
    #rechercher parmi tous les ensembles de correcction de biais
    for x in INT_NDRE_VTR_BIASES_DIMRED:
        #le group de biais dont tous les élements sont dans biases
        to_search = "_".join(sorted(list(dict.fromkeys(x.split(",")))))
        if to_search == biases:
            return x

"""
This rule launches the R script to apply normalization, dimensions reduction and evaluation of parameters.
"""
rule int_norm_dimred_ge:
    input:
        int_rda_file = int_norm_dimred_input_ge_snakemake
    output:
        int_ndre_Eval_rda_file = os.path.normpath("{output_int_norm_dimred_dir_ge}/GROUPED_ANALYSIS/INTEGRATED/{name_int}/{int_ndre_norm_vtr}/{int_ndre_dimred_vtr}/{name_int}_{int_ndre_norm_vtr}_{int_ndre_dimred_vtr}.rda")
    params:
        sing_int_bind = int_norm_dimred_params_sing,
        int_input_rda = int_norm_dimred_input_ge_R,
        int_ndre_vtr_biases_norm = int_norm_dimred_biases_norm,
        int_ndre_vtr_biases_dimred = int_norm_dimred_biases_dimred
    log:
        "logs/int_norm_dimred_ge{output_int_norm_dimred_dir_ge}/GROUPED_ANALYSIS/INTEGRATED/{name_int}/{int_ndre_norm_vtr}/{int_ndre_dimred_vtr}/{name_int}_{int_ndre_norm_vtr}_{int_ndre_dimred_vtr}.log"
    benchmark:
        "benchmark/int_norm_dimred_ge{output_int_norm_dimred_dir_ge}/GROUPED_ANALYSIS/INTEGRATED/{name_int}/{int_ndre_norm_vtr}/{int_ndre_dimred_vtr}/{name_int}_{int_ndre_norm_vtr}_{int_ndre_dimred_vtr}.tsv"
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, input, attempt: 30720 * len(input) + (attempt-1) * 40960),
        time_min = (lambda wildcards, input, attempt: 180 * len(input) + (attempt-1) * 360)
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        singularity exec --contain -B $TMPDIR:/$HOME -B $TMPDIR:/tmp {params.sing_int_bind} \
        {SINGULARITY_ENV_INT} \
        Rscript {PIPELINE_FOLDER}/scripts/Integration_part1.R \
        --input.list.rda {params.int_input_rda} \
        --output.dir.int {wildcards.output_int_norm_dimred_dir_ge}/ \
        --name.int {wildcards.name_int} \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --nthreads {threads} \
        --pipeline.path {PIPELINE_FOLDER} \
        --eval.markers {INT_NDRE_EVAL_MARKERS} \
        --min.cells {INT_NDRE_MIN_CELLS} \
        --integration.method {INT_NDRE_INT_METHOD}  \
        --vtr.batch {INT_NDRE_VTR_BATCH} \
        --features.n {INT_NDRE_FEATURES_N} \
        --norm.method {INT_NDRE_NORM_METHOD} \
        --HVG.FindVariableFeaturesMix {INT_NDRE_HVG_FINDVARIABLEFEATURESMIX} \
        --regex.genes.to.remove.from.HVG {INT_NDRE_REGEX_REMOVE_HVG} \
        --dimred.method {INT_NDRE_DIMRED_METHOD} \
        --vtr.biases.norm {params.int_ndre_vtr_biases_norm} \
        --vtr.biases.dimred {params.int_ndre_vtr_biases_dimred} \
        --vtr.scale {INT_NDRE_VTR_SCALE} \
        --dims.max {INT_NDRE_DIM_MAX} \
        --skip.eval_dims_res {INT_NDRE_SKIP_EVALDIMRES} \
        --eval.dims.max {INT_NDRE_EVAL_DIM_MAX} \
        --eval.dims.min {INT_NDRE_EVAL_DIM_MIN} \
        --eval.dims.steps {INT_NDRE_EVAL_DIM_STEPS} \
        --eval.res.max {INT_NDRE_EVAL_RES_MAX} \
        --eval.res.min {INT_NDRE_EVAL_RES_MIN} \
        --eval.res.steps {INT_NDRE_EVAL_RES_STEPS} \
        --eval.pt.size {INT_NDRE_EVAL_PTSIZE} \
        --metadata.file {INT_NDRE_METADATA_FILE} &> {log}
        """
