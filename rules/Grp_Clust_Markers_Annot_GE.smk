"""
##########################################################################
This rule make the clustering, to find markers genes and to apply annotations of cell types in single-cell grpegrated RNA-seq.
##########################################################################
"""

wildcard_constraints:
    name_grp = "|".join(GRP_NDRE_NAME_GRP)

"""
This function allows to determine the input .rda file.
"""
def grp_clust_markers_annot_input_ge(wildcards):
    return dic_GRP_CMA_INFO[wildcards.name_grp]['GRP_CMA_INPUT_RDA']

"""
This function allows to determine the singularity binding parameters.
"""
def grp_clust_markers_annot_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_GRP_CMA_INFO[wildcards.name_grp]['GRP_CMA_INPUT_RDA'])
    output_folder = wildcards.output_grp_clust_markers_annot_dir_ge
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder) + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if GRP_CMA_MARKFILE != "NULL":
        for markfile in list(dict.fromkeys(GRP_CMA_MARKFILE.split(","))):
            markfile = os.path.dirname(markfile)
            concat = concat + " -B " + markfile + ":" + os.path.normpath("/WORKDIR/" + markfile)
    if GRP_CMA_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(GRP_CMA_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches the R script to make the clustering, to find markers genes and to apply annotations of cell types.
"""
rule grp_clust_markers_annot_ge:
    input:
        grp_cma_file = grp_clust_markers_annot_input_ge
    output:
        grp_cma_rda_file = os.path.normpath("{output_grp_clust_markers_annot_dir_ge}" + "/" + GRP_CMA_CLUST_FOLDER + "/" + "{name_grp}" + "{grp_complement}" + str(GRP_CMA_KEEP_DIM) + "_" + str(GRP_CMA_KEEP_RES) + ".rda")
    params:
        sing_grp_bind = grp_clust_markers_annot_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        grp_input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        grp_output_folder = os.path.normpath("/WORKDIR/" + "{output_grp_clust_markers_annot_dir_ge}") + "/",
        SING_GRP_CMA_MARKFILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in GRP_CMA_MARKFILE.split(',')]) if GRP_CMA_MARKFILE != "NULL" else "NULL",
        SING_GRP_CMA_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in GRP_CMA_METADATA_FILE.split(',')]) if GRP_CMA_METADATA_FILE != "NULL" else "NULL"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(5120 + attempt * 3072, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 120, 200))
    shell:
        """
        singularity exec --contain {params.sing_grp_bind} \
        {SINGULARITY_ENV} \
        Rscript {params.pipeline_folder}/scripts/Grouped_analysis_part2.R \
        --input.rda.grp {params.grp_input_rda} \
        --output.dir.grp {params.grp_output_folder} \
        --author.name {GRP_CMA_AUTHOR_NAME} \
        --author.mail {GRP_CMA_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --markfile  {params.SING_GRP_CMA_MARKFILE} \
        --keep.dims {GRP_CMA_KEEP_DIM} \
        --keep.res {GRP_CMA_KEEP_RES} \
        --cfr.minscore {GRP_CMA_CFR_MINSCORE} \
        --sr.minscore {GRP_CMA_SR_MINSCORE} \
        --metadata.file {params.SING_GRP_CMA_METADATA_FILE}
        """