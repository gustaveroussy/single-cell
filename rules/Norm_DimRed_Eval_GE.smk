"""
##########################################################################
This rule make the normalization, dimensions reduction and evaluation in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file.
"""
def norm_dimred_input_ge(wildcards):
    for key, value in dic_NDRE_INFO.items():
        if value["NDRE_OUTPUT_DIR"] == wildcards.output_norm_dimred_dir_ge :
            return key
            #key is the input value in the yaml parameter file.

"""
This function allows to determine the singularity binding parameters.
"""
def norm_dimred_params_sing(wildcards):
    for key, value in dic_NDRE_INFO.items():
        if value["NDRE_OUTPUT_DIR"] == wildcards.output_norm_dimred_dir_ge :
            rda_folder = os.path.dirname(key)
    output_folder = wildcards.output_norm_dimred_dir_ge
    concat = " -B " + PIPELINE_FOLDER + ","+ rda_folder + "," + output_folder
    if NDRE_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(NDRE_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + "," + metadatafile
    return concat

"""
These functions allows to find and format the biases for the R script from wilcards.
"""
def norm_dimred_biases_norm(wildcards):
    #recherche si c'est dans ndre_norm_vtr ou s'il n'y en a pas.
    for prefix in ("SCTransform_", "LogNormalize_"):
        if wildcards.ndre_norm_vtr.startswith(prefix):
            biases = wildcards.ndre_norm_vtr.replace(prefix, "")
            break
    else:
        return "NULL"
    #rechercher parmi tous les ensembles de correcction de biais
    for x in NDRE_VTR_BIASES_NORM:
        #le group de biais dont tous les élements sont dans biases
        to_search = "_".join(sorted(list(dict.fromkeys(x.split(",")))))
        if to_search == biases:
            return x

def norm_dimred_biases_dimred(wildcards):
    #recherche si c'est dans ndre_dimred_vtr ou s'il n'y en a pas.
    for prefix in ("scbfa_", "bpca_"):
        if wildcards.ndre_dimred_vtr.startswith(prefix):
            biases = wildcards.ndre_dimred_vtr.replace(prefix, "")
            break
    else:
        return "NULL"
    #rechercher parmi tous les ensembles de correcction de biais
    for x in NDRE_VTR_BIASES_DIMRED:
        #le group de biais dont tous les élements sont dans biases
        to_search = "_".join(sorted(list(dict.fromkeys(x.split(",")))))
        if to_search == biases:
            return x
"""
This rule launches the R script to apply normalization, dimensions reduction and evaluation of parameters.
"""
rule norm_dimred_ge:
    input:
        rda_file = norm_dimred_input_ge
    output:
        ndre_Eval_rda_file = os.path.normpath("{output_norm_dimred_dir_ge}/{ndre_norm_vtr}/{ndre_dimred_vtr}/{sample_name_ge}_{ndre_norm_vtr}_{ndre_dimred_vtr}.rda")
    params:
        sing_bind = norm_dimred_params_sing,
        ndre_vtr_biases_norm = norm_dimred_biases_norm,
        ndre_vtr_biases_dimred = norm_dimred_biases_dimred
    log:
        "logs/norm_dimred_ge{output_norm_dimred_dir_ge}/{ndre_norm_vtr}/{ndre_dimred_vtr}/{sample_name_ge}_{ndre_norm_vtr}_{ndre_dimred_vtr}.log"
    benchmark:
        "benchmark/norm_dimred_ge{output_norm_dimred_dir_ge}/{ndre_norm_vtr}/{ndre_dimred_vtr}/{sample_name_ge}_{ndre_norm_vtr}_{ndre_dimred_vtr}.tsv"
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: 5210 + attempt * 10240),
        time_min = (lambda wildcards, attempt: attempt * 360)
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        singularity exec --contain -B $TMPDIR:/$HOME -B $TMPDIR:/tmp {params.sing_bind} \
        {SINGULARITY_ENV} \
        Rscript {PIPELINE_FOLDER}/scripts/pipeline_part3.R \
        --input.rda.ge {input[0]} \
        --output.dir.ge {wildcards.output_norm_dimred_dir_ge}/ \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --nthreads {threads} \
        --pipeline.path {PIPELINE_FOLDER} \
        --eval.markers {NDRE_EVAL_MARKERS} \
        --features.n {NDRE_FEATURES_N} \
        --norm.method {NDRE_NORM_METHOD} \
        --HVG.FindVariableFeaturesMix {NDRE_HVG_FINDVARIABLEFEATURESMIX} \
        --regex.genes.to.remove.from.HVG {NDRE_REGEX_REMOVE_HVG} \
        --dimred.method {NDRE_DIMRED_METHOD} \
        --vtr.biases.norm {params.ndre_vtr_biases_norm} \
        --vtr.biases.dimred {params.ndre_vtr_biases_dimred} \
        --vtr.scale {NDRE_VTR_SCALE} \
        --dims.max {NDRE_DIM_MAX} \
        --skip.eval_dims_res {NDRE_SKIP_EVALDIMRES} \
        --eval.dims.max {NDRE_EVAL_DIM_MAX} \
        --eval.dims.min {NDRE_EVAL_DIM_MIN} \
        --eval.dims.steps {NDRE_EVAL_DIM_STEPS} \
        --eval.res.max {NDRE_EVAL_RES_MAX} \
        --eval.res.min {NDRE_EVAL_RES_MIN} \
        --eval.res.steps {NDRE_EVAL_RES_STEPS} \
        --eval.pt.size {NDRE_EVAL_PTSIZE} \
        --metadata.file {NDRE_METADATA_FILE} &> {log}
        """
