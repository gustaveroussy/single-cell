"""
##########################################################################
This rule make the grpegration, normalization, dimensions reduction and evaluation in single-cell RNA-seq.
##########################################################################
"""

# wildcard_constraints:
#    name_grp = "|".join(GRP_NDRE_NAME_GRP)

"""
This function allows to determine the input .rda file.
"""
def grp_norm_dimred_input_ge_snakemake(wildcards):
    return dic_GRP_NDRE_INFO[wildcards.name_grp]['GRP_NDRE_INPUT_LIST_RDA'].split(',')

def grp_norm_dimred_input_ge_R(wildcards):
    return ','.join([os.path.normpath("/WORKDIR/" + x) for x in dic_GRP_NDRE_INFO[wildcards.name_grp]['GRP_NDRE_INPUT_LIST_RDA'].split(',')])

"""
This function allows to determine the singularity binding parameters.
"""
def grp_norm_dimred_params_sing(wildcards):
    sing_rda_folder = ''.join([" -B " + os.path.dirname(x) + ":" + os.path.normpath("/WORKDIR/" + os.path.dirname(x)) for x in dic_GRP_NDRE_INFO[wildcards.name_grp]['GRP_NDRE_INPUT_LIST_RDA'].split(',')])
    output_folder = wildcards.output_grp_norm_dimred_dir_ge
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + sing_rda_folder + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if GRP_NDRE_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(GRP_NDRE_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches the R script to apply normalization, dimensions reduction and evaluation of parameters.
"""
rule grp_norm_dimred_ge:
    input:
        grp_rda_file = grp_norm_dimred_input_ge_snakemake
    output:
        grp_ndre_Eval_rda_file = os.path.normpath("{output_grp_norm_dimred_dir_ge}" + "/GROUPED_ANALYSIS/NO_INTEGRATED/{name_grp}/" + GRP_NDRE_NORM_VTR + "/" + GRP_NDRE_DIMRED_VTR + "/" + "{name_grp}_" + str(GRP_NDRE_NORM_VTR) + "_" + str(GRP_NDRE_DIMRED_VTR) + ".rda")
    params:
        sing_grp_bind = grp_norm_dimred_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        grp_input_rda = grp_norm_dimred_input_ge_R,
        grp_output_folder = os.path.normpath("/WORKDIR/" + "{output_grp_norm_dimred_dir_ge}") + "/",
        SING_GRP_NDRE_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in GRP_NDRE_METADATA_FILE.split(',')]) if GRP_NDRE_METADATA_FILE != "NULL" else "NULL"
    threads:
        4
    resources:
        mem_mb = (lambda wildcards, attempt: min(5120 * len(dic_GRP_NDRE_INFO[wildcards.name_grp]['GRP_NDRE_INPUT_LIST_RDA'].split(',')) + (attempt-1) * 10240, 184320)),
        time_min = (lambda wildcards, attempt: min(60 * len(dic_GRP_NDRE_INFO[wildcards.name_grp]['GRP_NDRE_INPUT_LIST_RDA'].split(',')) + attempt * 120, 4320))
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home -B $TMP_DIR:/tmp {params.sing_grp_bind} \
        {SINGULARITY_ENV} \
        Rscript {params.pipeline_folder}/scripts/Grouped_analysis_part1.R \
        --input.list.rda {params.grp_input_rda} \
        --output.dir.grp {params.grp_output_folder} \
        --name.grp {wildcards.name_grp} \
        --author.name {GRP_NDRE_AUTHOR_NAME} \
        --author.mail {GRP_NDRE_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --eval.markers {GRP_NDRE_EVAL_MARKERS} \
        --min.cells {GRP_NDRE_MIN_CELLS} \
        --keep.norm {GRP_NDRE_KEEP_NORM}  \
        --features.n {GRP_NDRE_FEATURES_N} \
        --norm.method {GRP_NDRE_NORM_METHOD} \
        --dimred.method {GRP_NDRE_DIMRED_METHOD} \
        --vtr.biases {GRP_NDRE_VTR_BIASES} \
        --vtr.scale {GRP_NDRE_VTR_SCALE} \
        --dims.max {GRP_NDRE_DIM_MAX} \
        --dims.min {GRP_NDRE_DIM_MIN} \
        --dims.steps {GRP_NDRE_DIM_STEPS} \
        --res.max {GRP_NDRE_RES_MAX} \
        --res.min {GRP_NDRE_RES_MIN} \
        --res.steps {GRP_NDRE_RES_STEPS} \
        --metadata.file {params.SING_GRP_NDRE_METADATA_FILE} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
