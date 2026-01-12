"""
##########################################################################
This rule make the grpegration, normalization, dimensions reduction and evaluation in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file and format it for R script.
"""
def grp_norm_dimred_input_ge_R(wildcards):
    for key, value in dic_GRP_NDRE_INFO.items():
        if value["GRP_NDRE_OUTPUT_DIR"] == (wildcards.output_grp_norm_dimred_dir_ge + "/GROUPED_ANALYSIS/NO_INTEGRATED/" + wildcards.name_grp) :
            return key
            #key is the input value in the yaml parameter file.

"""
This function allows to determine the input .rda file.
"""
def grp_norm_dimred_input_ge_snakemake(wildcards):
    return (grp_norm_dimred_input_ge_R(wildcards)).split(',')

"""
This function allows to determine the singularity binding parameters.
"""
def grp_norm_dimred_params_sing(wildcards):
    for key, value in dic_GRP_NDRE_INFO.items():
        if value["GRP_NDRE_OUTPUT_DIR"] == (wildcards.output_grp_norm_dimred_dir_ge + "/GROUPED_ANALYSIS/NO_INTEGRATED/" + wildcards.name_grp) :
            sing_rda_folder = ','.join([os.path.dirname(x) for x in key.split(',')])
            #key is the input value in the yaml parameter file.
    output_folder = wildcards.output_grp_norm_dimred_dir_ge
    concat = " -B " + PIPELINE_FOLDER + "," + sing_rda_folder + "," + output_folder
    if GRP_NDRE_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(GRP_NDRE_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + "," + metadatafile
    return concat

"""
These functions allows to find and format the biases for the R script from wilcards.
"""
def grp_norm_dimred_biases_norm(wildcards):
    #recherche si c'est dans grp_ndre_norm_vtr ou s'il n'y en a pas.
    for prefix in ("SCTransform_", "LogNormalize_"):
        if wildcards.grp_ndre_norm_vtr.startswith(prefix):
            biases = wildcards.grp_ndre_norm_vtr.replace(prefix, "")
            break
    else:
        return "NULL"
    #rechercher parmi tous les ensembles de correcction de biais
    for x in GRP_NDRE_VTR_BIASES_NORM:
        #le group de biais dont tous les élements sont dans biases
        to_search = "_".join(sorted(list(dict.fromkeys(x.split(",")))))
        if to_search == biases:
            return x

def grp_norm_dimred_biases_dimred(wildcards):
    #recherche si c'est dans grp_ndre_dimred_vtr ou s'il n'y en a pas.
    for prefix in ("scbfa_", "bpca_"):
        if wildcards.grp_ndre_dimred_vtr.startswith(prefix):
            biases = wildcards.grp_ndre_dimred_vtr.replace(prefix, "")
            break
    else:
        return "NULL"
    #rechercher parmi tous les ensembles de correcction de biais
    for x in GRP_NDRE_VTR_BIASES_DIMRED:
        #le group de biais dont tous les élements sont dans biases
        to_search = "_".join(sorted(list(dict.fromkeys(x.split(",")))))
        if to_search == biases:
            return x

"""
This rule launches the R script to apply normalization, dimensions reduction and evaluation of parameters.
"""
rule grp_norm_dimred_ge:
    input:
        grp_rda_file = grp_norm_dimred_input_ge_snakemake
    output:
        grp_ndre_Eval_rda_file = os.path.normpath("{output_grp_norm_dimred_dir_ge}/GROUPED_ANALYSIS/NO_INTEGRATED/{name_grp}/{grp_ndre_norm_vtr}/{grp_ndre_dimred_vtr}/{name_grp}_{grp_ndre_norm_vtr}_{grp_ndre_dimred_vtr}.rda")
    params:
        sing_grp_bind = grp_norm_dimred_params_sing,
        grp_input_rda = grp_norm_dimred_input_ge_R,
        grp_ndre_vtr_biases_norm = grp_norm_dimred_biases_norm,
        grp_ndre_vtr_biases_dimred = grp_norm_dimred_biases_dimred
    log:
        "logs/grp_norm_dimred_ge{output_grp_norm_dimred_dir_ge}/GROUPED_ANALYSIS/NO_INTEGRATED/{name_grp}/{grp_ndre_norm_vtr}/{grp_ndre_dimred_vtr}/{name_grp}_{grp_ndre_norm_vtr}_{grp_ndre_dimred_vtr}.log"
    benchmark:
        "benchmark/grp_norm_dimred_ge{output_grp_norm_dimred_dir_ge}/GROUPED_ANALYSIS/NO_INTEGRATED/{name_grp}/{grp_ndre_norm_vtr}/{grp_ndre_dimred_vtr}/{name_grp}_{grp_ndre_norm_vtr}_{grp_ndre_dimred_vtr}.tsv"
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, input, attempt: 30720 * len(input) + (attempt-1) * 40960),
        time_min = (lambda wildcards, input, attempt: 180 * len(input) + (attempt-1) * 360)
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        singularity exec --no-home -B $TMPDIR:/tmp {params.sing_grp_bind} \
        {SINGULARITY_ENV} \
        Rscript {PIPELINE_FOLDER}/scripts/Grouped_analysis_part1.R \
        --input.list.rda {params.grp_input_rda} \
        --output.dir.grp {wildcards.output_grp_norm_dimred_dir_ge}/ \
        --name.grp {wildcards.name_grp} \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --nthreads {threads} \
        --pipeline.path {PIPELINE_FOLDER} \
        --eval.markers {GRP_NDRE_EVAL_MARKERS} \
        --min.cells {GRP_NDRE_MIN_CELLS} \
        --individual.norm {GRP_NDRE_INDIV_NORM}  \
        --features.n {GRP_NDRE_FEATURES_N} \
        --norm.method {GRP_NDRE_NORM_METHOD} \
        --HVG.FindVariableFeaturesMix {GRP_NDRE_HVG_FINDVARIABLEFEATURESMIX} \
        --regex.genes.to.remove.from.HVG {GRP_NDRE_REGEX_REMOVE_HVG} \
        --dimred.method {GRP_NDRE_DIMRED_METHOD} \
        --vtr.biases.norm {params.grp_ndre_vtr_biases_norm} \
        --vtr.biases.dimred {params.grp_ndre_vtr_biases_dimred} \
        --vtr.scale {GRP_NDRE_VTR_SCALE} \
        --dims.max {GRP_NDRE_DIM_MAX} \
        --skip.eval_dims_res {GRP_NDRE_SKIP_EVALDIMRES} \
        --eval.dims.max {GRP_NDRE_EVAL_DIM_MAX} \
        --eval.dims.min {GRP_NDRE_EVAL_DIM_MIN} \
        --eval.dims.steps {GRP_NDRE_EVAL_DIM_STEPS} \
        --eval.res.max {GRP_NDRE_EVAL_RES_MAX} \
        --eval.res.min {GRP_NDRE_EVAL_RES_MIN} \
        --eval.res.steps {GRP_NDRE_EVAL_RES_STEPS} \
        --eval.pt.size {GRP_NDRE_EVAL_PTSIZE} \
        --metadata.file {GRP_NDRE_METADATA_FILE} &> {log}
        """
