"""
##########################################################################
This rule make the clustering, to find markers genes and to apply annotations of cell types in single-cell integrated RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file.
"""
def int_clust_markers_annot_input_ge(wildcards):
    #return dic_INT_CMA_INFO[wildcards.name_int]['INT_CMA_INPUT_RDA']
    for key, value in dic_INT_CMA_INFO.items():
        if value["INT_CMA_OUTPUT_DIR"] == wildcards.output_int_clust_markers_annot_dir_ge :
            return key
            #key is the input value in the yaml parameter file.

"""
This function allows to determine the singularity binding parameters.
"""
def int_clust_markers_annot_params_sing(wildcards):
    #rda_folder = os.path.dirname(dic_INT_CMA_INFO[wildcards.name_int]['INT_CMA_INPUT_RDA'])
    for key, value in dic_INT_CMA_INFO.items():
        if value["INT_CMA_OUTPUT_DIR"] == wildcards.output_int_clust_markers_annot_dir_ge :
            rda_folder = os.path.dirname(key)
    output_folder = wildcards.output_int_clust_markers_annot_dir_ge
    concat = " -B " + PIPELINE_FOLDER + "," + rda_folder + "," + output_folder
    if INT_CMA_MARKFILE != "NULL":
        for markfile in list(dict.fromkeys(INT_CMA_MARKFILE.split(","))):
            markfile = os.path.dirname(markfile)
            concat = concat + "," + markfile
    if INT_CMA_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(INT_CMA_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + "," + metadatafile
    if INT_CMA_CUSTOM_SCE_REF != "NULL":
        for custom_cse_ref in list(dict.fromkeys(INT_CMA_CUSTOM_SCE_REF.split(","))):
            custom_cse_ref = os.path.dirname(custom_cse_ref)
            concat = concat + "," + custom_cse_ref
    if INT_CMA_CUSTOM_MARKERS_REF != "NULL":
        for custom_marker_ref in list(dict.fromkeys(INT_CMA_CUSTOM_MARKERS_REF.split(","))):
            custom_marker_ref = os.path.dirname(custom_marker_ref)
            concat = concat + "," + custom_marker_ref
    return concat

"""
This rule launches the R script to make the clustering, to find markers genes and to apply annotations of cell types.
"""
rule int_clust_markers_annot_ge:
    input:
        int_cma_file = int_clust_markers_annot_input_ge
    output:
        int_cma_rda_file = os.path.normpath("{output_int_clust_markers_annot_dir_ge}/{int_clust_folders}/{name_int}{int_complement}{int_keep_dim_res}.rda")
    params:
        sing_int_bind = int_clust_markers_annot_params_sing,
        keep_dim = lambda wildcards: wildcards.int_keep_dim_res.split('_')[0],
        keep_res = lambda wildcards: wildcards.int_keep_dim_res.split('_')[1]
    log:
        "logs/int_clust_markers_annot_ge{output_int_clust_markers_annot_dir_ge}/{int_clust_folders}/{name_int}{int_complement}{int_keep_dim_res}.log"
    benchmark:
        "benchmark/int_clust_markers_annot_ge{output_int_clust_markers_annot_dir_ge}/{int_clust_folders}/{name_int}{int_complement}{int_keep_dim_res}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(5120 + attempt * 51200, 716800)),
        time_min = (lambda wildcards, attempt: min(attempt * 1440, 10080))
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        rsync -az -c {PIPELINE_FOLDER}/resources/DATABASE/ref_annot_cache $TMPDIR/ref_annot_cache  && \
        singularity exec --no-home -B $TMPDIR:/tmp -B $TMPDIR/ref_annot_cache:$HOME {params.sing_int_bind} \
        {SINGULARITY_ENV} \
        Rscript {PIPELINE_FOLDER}/scripts/Integration_part2.R \
        --input.rda.int {input[0]} \
        --output.dir.int {wildcards.output_int_clust_markers_annot_dir_ge}/ \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --nthreads {threads} \
        --pipeline.path {PIPELINE_FOLDER} \
        --markfile {INT_CMA_MARKFILE} \
        --markers.pt.size {INT_CMA_MARKERS_PTSIZE} \
        --markers.order {INT_CMA_MARKERS_ORDER} \
        --custom.sce.ref {INT_CMA_CUSTOM_SCE_REF} \
        --custom.markers.ref {INT_CMA_CUSTOM_MARKERS_REF} \
        --skip.technical_plots {INT_CMA_SKIP_TECHNICAL} \
        --skip.annotation {INT_CMA_SKIP_ANNOT} \
        --skip.markers_identification {INT_CMA_SKIP_MARKERS_IDENT} \
        --keep.dims {params.keep_dim} \
        --keep.res {params.keep_res} \
        --cfr.minscore {INT_CMA_CFR_MINSCORE} \
        --sr.minscore {INT_CMA_SR_MINSCORE} \
        --metadata.file {INT_CMA_METADATA_FILE} &> {log}
        """
