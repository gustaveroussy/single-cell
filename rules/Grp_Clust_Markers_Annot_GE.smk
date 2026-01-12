"""
##########################################################################
This rule make the clustering, to find markers genes and to apply annotations of cell types in single-cell grpegrated RNA-seq.
##########################################################################
"""

"""
This function allows to determine the input .rda file.
"""
def grp_clust_markers_annot_input_ge(wildcards):
    #return dic_GRP_CMA_INFO[wildcards.name_grp]['GRP_CMA_INPUT_RDA']
    for key, value in dic_GRP_CMA_INFO.items():
        if value["GRP_CMA_OUTPUT_DIR"] == wildcards.output_grp_clust_markers_annot_dir_ge :
            return key
            #key is the input value in the yaml parameter file.

"""
This function allows to determine the singularity binding parameters.
"""
def grp_clust_markers_annot_params_sing(wildcards):
    #rda_folder = os.path.dirname(dic_GRP_CMA_INFO[wildcards.name_grp]['GRP_CMA_INPUT_RDA'])
    for key, value in dic_GRP_CMA_INFO.items():
        if value["GRP_CMA_OUTPUT_DIR"] == wildcards.output_grp_clust_markers_annot_dir_ge :
            rda_folder = os.path.dirname(key)
            #key is the input value in the yaml parameter file.
    output_folder = wildcards.output_grp_clust_markers_annot_dir_ge
    concat = " -B " + PIPELINE_FOLDER + "," + rda_folder + "," + output_folder
    if GRP_CMA_MARKFILE != "NULL":
        for markfile in list(dict.fromkeys(GRP_CMA_MARKFILE.split(","))):
            markfile = os.path.dirname(markfile)
            concat = concat + "," + markfile
    if GRP_CMA_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(GRP_CMA_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + "," + metadatafile
    if GRP_CMA_CUSTOM_SCE_REF != "NULL":
        for custom_cse_ref in list(dict.fromkeys(GRP_CMA_CUSTOM_SCE_REF.split(","))):
            custom_cse_ref = os.path.dirname(custom_cse_ref)
            concat = concat + "," + custom_cse_ref
    if GRP_CMA_CUSTOM_MARKERS_REF != "NULL":
        for custom_marker_ref in list(dict.fromkeys(GRP_CMA_CUSTOM_MARKERS_REF.split(","))):
            custom_marker_ref = os.path.dirname(custom_marker_ref)
            concat = concat + "," + custom_marker_ref
    return concat

"""
This rule launches the R script to make the clustering, to find markers genes and to apply annotations of cell types.
"""
rule grp_clust_markers_annot_ge:
    input:
        grp_cma_file = grp_clust_markers_annot_input_ge
    output:
        grp_cma_rda_file = os.path.normpath("{output_grp_clust_markers_annot_dir_ge}/{grp_clust_folders}/{name_grp}{grp_complement}{grp_keep_dim_res}.rda")
    params:
        sing_grp_bind = grp_clust_markers_annot_params_sing,
        keep_dim = lambda wildcards: wildcards.grp_keep_dim_res.split('_')[0],
        keep_res = lambda wildcards: wildcards.grp_keep_dim_res.split('_')[1]
    log:
        "logs/grp_clust_markers_annot_ge{output_grp_clust_markers_annot_dir_ge}/{grp_clust_folders}/{name_grp}{grp_complement}{grp_keep_dim_res}.log"
    benchmark:
        "benchmark/grp_clust_markers_annot_ge{output_grp_clust_markers_annot_dir_ge}/{grp_clust_folders}/{name_grp}{grp_complement}{grp_keep_dim_res}.tsv"
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
        singularity exec --no-home -B $TMPDIR:/tmp -B $TMPDIR/ref_annot_cache:$HOME {params.sing_grp_bind} \
        {SINGULARITY_ENV} \
        Rscript {PIPELINE_FOLDER}/scripts/Grouped_analysis_part2.R \
        --input.rda.grp {input[0]} \
        --output.dir.grp {wildcards.output_grp_clust_markers_annot_dir_ge}/ \
        --author.name "{AUTHOR_NAME}" \
        --author.mail "{AUTHOR_MAIL}" \
        --nthreads {threads} \
        --pipeline.path {PIPELINE_FOLDER} \
        --markfile  {GRP_CMA_MARKFILE} \
        --markers.pt.size  {GRP_CMA_MARKERS_PTSIZE} \
        --markers.order  {GRP_CMA_MARKERS_ORDER} \
        --custom.sce.ref {GRP_CMA_CUSTOM_SCE_REF} \
        --custom.markers.ref {GRP_CMA_CUSTOM_MARKERS_REF} \
        --skip.technical_plots {GRP_CMA_SKIP_TECHNICAL} \
        --skip.annotation {GRP_CMA_SKIP_ANNOT} \
        --skip.markers_identification {GRP_CMA_SKIP_MARKERS_IDENT} \
        --keep.dims {params.keep_dim} \
        --keep.res {params.keep_res} \
        --cfr.minscore {GRP_CMA_CFR_MINSCORE} \
        --sr.minscore {GRP_CMA_SR_MINSCORE} \
        --metadata.file {GRP_CMA_METADATA_FILE} &> {log}
        """
