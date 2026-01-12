"""
##########################################################################
This rule make the clustering, to find markers genes and to apply annotations of cell types in single-cell RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file.
"""
def clust_markers_annot_input_ge(wildcards):
    #return dic_CMA_INFO[wildcards.sample_name_ge]['CMA_INPUT_RDA']
    for key, value in dic_CMA_INFO.items():
        if value["CMA_OUTPUT_DIR"] == wildcards.output_clust_markers_annot_dir_ge :
            return key
            #key is the input value in the yaml parameter file.

"""
This function allows to determine the singularity binding parameters.
"""
def clust_markers_annot_params_sing(wildcards):
    #rda_folder = os.path.dirname(dic_CMA_INFO[wildcards.sample_name_ge]['CMA_INPUT_RDA'])
    for key, value in dic_CMA_INFO.items():
        if value["CMA_OUTPUT_DIR"] == wildcards.output_clust_markers_annot_dir_ge :
            rda_folder = os.path.dirname(key)
    output_folder = wildcards.output_clust_markers_annot_dir_ge
    concat = " -B " + PIPELINE_FOLDER + ","+ rda_folder + "," + output_folder
    if CMA_MARKFILE != "NULL":
        for markfile in list(dict.fromkeys(CMA_MARKFILE.split(","))):
            markfile = os.path.dirname(markfile)
            concat = concat + "," + markfile
    if CMA_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(CMA_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + "," + metadatafile
    if CMA_CUSTOM_SCE_REF != "NULL":
        for custom_cse_ref in list(dict.fromkeys(CMA_CUSTOM_SCE_REF.split(","))):
            custom_cse_ref = os.path.dirname(custom_cse_ref)
            concat = concat + "," + custom_cse_ref
    if CMA_CUSTOM_MARKERS_REF != "NULL":
        for custom_marker_ref in list(dict.fromkeys(CMA_CUSTOM_MARKERS_REF.split(","))):
            custom_marker_ref = os.path.dirname(custom_marker_ref)
            concat = concat + "," + custom_marker_ref
    return concat

"""
This rule launches the R script to make the clustering, to find markers genes and to apply annotations of cell types.
"""
rule clust_markers_annot_ge:
    input:
        cma_file = clust_markers_annot_input_ge
    output:
        cma_rda_file = os.path.normpath("{output_clust_markers_annot_dir_ge}/{clust_folders}/{sample_name_ge}{complement}{keep_dim_res}.rda")
    params:
        sing_bind = clust_markers_annot_params_sing,
        keep_dim = lambda wildcards: wildcards.keep_dim_res.split('_')[0],
        keep_res = lambda wildcards: wildcards.keep_dim_res.split('_')[1]
    log:
        "logs/clust_markers_annot_ge{output_clust_markers_annot_dir_ge}/{clust_folders}/{sample_name_ge}{complement}{keep_dim_res}.log"
    benchmark:
        "benchmark/clust_markers_annot_ge{output_clust_markers_annot_dir_ge}/{clust_folders}/{sample_name_ge}{complement}{keep_dim_res}.tsv"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(10240 + attempt * 5120, 81920)),
        time_min = (lambda wildcards, attempt: min(attempt * 120, 420))
    shell:
        """
        TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
        trap "rm -r $TMPDIR" EXIT
        rsync -az -c {PIPELINE_FOLDER}/resources/DATABASE/ref_annot_cache $TMPDIR/ref_annot_cache  && \
        singularity exec --no-home -B $TMPDIR:/tmp -B $TMPDIR/ref_annot_cache:$HOME {params.sing_bind} \
        {SINGULARITY_ENV} \
        Rscript {PIPELINE_FOLDER}/scripts/pipeline_part4.R \
        --input.rda.ge {input[0]} \
        --output.dir.ge {wildcards.output_clust_markers_annot_dir_ge}/ \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --nthreads {threads} \
        --pipeline.path {PIPELINE_FOLDER} \
        --markfile  {CMA_MARKFILE} \
        --markers.pt.size  {CMA_MARKERS_PTSIZE} \
        --markers.order  {CMA_MARKERS_ORDER} \
        --custom.sce.ref {CMA_CUSTOM_SCE_REF} \
        --custom.markers.ref {CMA_CUSTOM_MARKERS_REF} \
        --skip.technical_plots {CMA_SKIP_TECHNICAL} \
        --skip.annotation {CMA_SKIP_ANNOT} \
        --skip.markers_identification {CMA_SKIP_MARKERS_IDENT} \
        --keep.dims {params.keep_dim} \
        --keep.res {params.keep_res} \
        --cfr.minscore {CMA_CFR_MINSCORE} \
        --sr.minscore {CMA_SR_MINSCORE} \
        --metadata.file {CMA_METADATA_FILE} &> {log}
        """
