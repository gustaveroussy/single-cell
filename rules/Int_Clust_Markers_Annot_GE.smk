"""
##########################################################################
This rule make the clustering, to find markers genes and to apply annotations of cell types in single-cell integrated RNA-seq.
##########################################################################
"""

wildcard_constraints:
    name_int = "|".join(INT_NDRE_NAME_INT)

"""
This function allows to determine the input .rda file.
"""
def int_clust_markers_annot_input_ge(wildcards):
    return dic_INT_CMA_INFO[wildcards.name_int]['INT_CMA_INPUT_RDA']

"""
This function allows to determine the singularity binding parameters.
"""
def int_clust_markers_annot_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_INT_CMA_INFO[wildcards.name_int]['INT_CMA_INPUT_RDA'])
    output_folder = wildcards.output_int_clust_markers_annot_dir_ge
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder) + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if INT_CMA_MARKFILE != "NULL":
        for markfile in list(dict.fromkeys(INT_CMA_MARKFILE.split(","))):
            markfile = os.path.dirname(markfile)
            concat = concat + " -B " + markfile + ":" + os.path.normpath("/WORKDIR/" + markfile)
    if INT_CMA_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(INT_CMA_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    if INT_CMA_CUSTOM_SCE_REF != "NULL":
        for custom_cse_ref in list(dict.fromkeys(INT_CMA_CUSTOM_SCE_REF.split(","))):
            custom_cse_ref = os.path.dirname(custom_cse_ref)
            concat = concat + " -B " + custom_cse_ref + ":" + os.path.normpath("/WORKDIR/" + custom_cse_ref)
    if INT_CMA_CUSTOM_MARKERS_REF != "NULL":
        for custom_marker_ref in list(dict.fromkeys(INT_CMA_CUSTOM_MARKERS_REF.split(","))):
            custom_marker_ref = os.path.dirname(custom_marker_ref)
            concat = concat + " -B " + custom_marker_ref + ":" + os.path.normpath("/WORKDIR/" + custom_marker_ref)
    return concat

"""
This rule launches the R script to make the clustering, to find markers genes and to apply annotations of cell types.
"""
rule int_clust_markers_annot_ge:
    input:
        int_cma_file = int_clust_markers_annot_input_ge
    output:
        int_cma_rda_file = os.path.normpath("{output_int_clust_markers_annot_dir_ge}" + "/" + INT_CMA_CLUST_FOLDER + "/" + "{name_int}" + "{int_complement}" + str(INT_CMA_KEEP_DIM) + "_" + str(INT_CMA_KEEP_RES) + ".rda")
    params:
        sing_int_bind = int_clust_markers_annot_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        int_input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        int_output_folder = os.path.normpath("/WORKDIR/" + "{output_int_clust_markers_annot_dir_ge}") + "/",
        SING_INT_CMA_MARKFILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in INT_CMA_MARKFILE.split(',')]) if INT_CMA_MARKFILE != "NULL" else "NULL",
        SING_INT_CMA_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in INT_CMA_METADATA_FILE.split(',')]) if INT_CMA_METADATA_FILE != "NULL" else "NULL",
        SING_INT_CMA_CUSTOM_SCE_REF = ','.join([os.path.normpath("/WORKDIR/" + x) for x in INT_CMA_CUSTOM_SCE_REF.split(',')]) if INT_CMA_CUSTOM_SCE_REF != "NULL" else "NULL",
        SING_INT_CMA_CUSTOM_MARKERS_REF = ','.join([os.path.normpath("/WORKDIR/" + x) for x in INT_CMA_CUSTOM_MARKERS_REF.split(',')]) if INT_CMA_CUSTOM_MARKERS_REF != "NULL" else "NULL"

    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(10240 + attempt * 5120, 102400)),
        time_min = (lambda wildcards, attempt: min(attempt * 120, 200))
    shell:
        """
        export TMPDIR={GLOBAL_TMP}
        TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
        singularity exec --no-home -B $TMP_DIR:/tmp -B $TMP_DIR:$HOME {params.sing_int_bind} \
        {SINGULARITY_ENV} \
        Rscript {params.pipeline_folder}/scripts/Integration_part2.R \
        --input.rda.int {params.int_input_rda} \
        --output.dir.int {params.int_output_folder} \
        --author.name {INT_CMA_AUTHOR_NAME} \
        --author.mail {INT_CMA_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --markfile {params.SING_INT_CMA_MARKFILE} \
        --custom.sce.ref {params.SING_INT_CMA_CUSTOM_SCE_REF} \
        --custom.markers.ref {params.SING_INT_CMA_CUSTOM_MARKERS_REF} \
        --keep.dims {INT_CMA_KEEP_DIM} \
        --keep.res {INT_CMA_KEEP_RES} \
        --cfr.minscore {INT_CMA_CFR_MINSCORE} \
        --sr.minscore {INT_CMA_SR_MINSCORE} \
        --metadata.file {params.SING_INT_CMA_METADATA_FILE} && \
        rm -r $TMP_DIR || rm -r $TMP_DIR
        """
