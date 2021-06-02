"""
##########################################################################
This rule make the clustering, to find markers genes and to apply annotations of cell types in single-cell RNA-seq.
##########################################################################
"""

wildcard_constraints:
    sample_name_ge=".+_GE"

"""
This function allows to determine the input .rda file.
"""
def clust_markers_annot_input_ge(wildcards):
    return dic_CMA_INFO[wildcards.sample_name_ge]['CMA_INPUT_RDA']

"""
This function allows to determine the singularity binding parameters.
"""
def clust_markers_annot_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_CMA_INFO[wildcards.sample_name_ge]['CMA_INPUT_RDA'])
    output_folder = wildcards.output_clust_markers_annot_dir_ge
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder) + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if CMA_MARKFILE != "NULL":
        CMA_MARKFILE_LIST = list(dict.fromkeys(CMA_MARKFILE.split(",")))
        for markfile in CMA_MARKFILE_LIST:
            markfile = os.path.dirname(markfile)
            concat = concat + " -B " + markfile + ":" + os.path.normpath("/WORKDIR/" + markfile)
    if CMA_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(CMA_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches the R script to make the clustering, to find markers genes and to apply annotations of cell types.
"""
rule clust_markers_annot_ge:
    input:
        cma_file = clust_markers_annot_input_ge
    output:
        cma_rda_file = os.path.normpath("{output_clust_markers_annot_dir_ge}" + "/" + CMA_CLUST_FOLDER + "/" + "{sample_name_ge}" + "{complement}" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda")
    params:
        sing_bind = clust_markers_annot_params_sing,
        pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
        input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
        output_folder = os.path.normpath("/WORKDIR/" + "{output_clust_markers_annot_dir_ge}") + "/",
        SING_CMA_MARKFILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in CMA_MARKFILE.split(',')]) if CMA_MARKFILE != "NULL" else "NULL",
        SING_CMA_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in CMA_METADATA_FILE.split(',')]) if CMA_METADATA_FILE != "NULL" else "NULL"
    threads:
        1
    resources:
        mem_mb = (lambda wildcards, attempt: min(3072 + attempt * 1024, 20480)),
        time_min = (lambda wildcards, attempt: min(attempt * 60, 200))
    shell:
        """
        singularity exec --contain {params.sing_bind} \
        {SINGULARITY_ENV} \
        Rscript {params.pipeline_folder}/scripts/pipeline_part4.R \
        --input.rda.ge {params.input_rda} \
        --output.dir.ge {params.output_folder} \
        --author.name {CMA_AUTHOR_NAME} \
        --author.mail {CMA_AUTHOR_MAIL} \
        --nthreads {threads} \
        --pipeline.path {params.pipeline_folder} \
        --markfile  {params.SING_CMA_MARKFILE} \
        --keep.dims {CMA_KEEP_DIM} \
        --keep.res {CMA_KEEP_RES} \
        --cfr.minscore {CMA_CFR_MINSCORE} \
        --sr.minscore {CMA_SR_MINSCORE} \
        --metadata.file {params.SING_CMA_METADATA_FILE}
        """
