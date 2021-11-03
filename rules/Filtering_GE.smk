"""
##########################################################################
These rules filter droplets according to quality-control and whether to filter doublets or not, in single-cell RNA-seq.
##########################################################################
"""

wildcard_constraints:
    sample_name_ge=".+_GE"

"""
This function allows to determine the input .rda file.
"""
def filtering_input_ge(wildcards):
    return dic_FILTER_INFO[wildcards.sample_name_ge]['FILTER_INPUT_RDA']

"""
This function allows to determine the singularity binding parameters.
"""
def filtering_params_sing(wildcards):
    rda_folder = os.path.dirname(dic_FILTER_INFO[wildcards.sample_name_ge]['FILTER_INPUT_RDA'])
    output_folder = wildcards.output_filtering_dir_ge + "/"
    concat = " -B " + PIPELINE_FOLDER + ":" + os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER) + " -B " + rda_folder + ":" + os.path.normpath("/WORKDIR/" + rda_folder) + " -B " + output_folder + ":" + os.path.normpath("/WORKDIR/" + output_folder)
    if FILERING_CC_SEURAT_FILE != "NULL": concat = concat + " -B " + os.path.dirname(FILERING_CC_SEURAT_FILE) + ":" + os.path.normpath("/WORKDIR/" + os.path.dirname(FILERING_CC_SEURAT_FILE))
    if FILERING_CC_CYCLONE_FILE != "NULL": concat = concat + " -B " + os.path.dirname(FILERING_CC_CYCLONE_FILE) + ":" + os.path.normpath("/WORKDIR/" + os.path.dirname(FILERING_CC_CYCLONE_FILE))
    if FILERING_METADATA_FILE != "NULL":
        for metadatafile in list(dict.fromkeys(FILERING_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + " -B " + metadatafile + ":" + os.path.normpath("/WORKDIR/" + metadatafile)
    return concat

"""
This rule launches the R script to filter droplets according to quality control and to filter the doublets.
"""
if FILERING_DOUBLET_FILTER_METHOD_NAME != "none":
    rule filtering_ge:
        input:
            rda_file = filtering_input_ge
        output:
            dbk_QC_hist_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_QChist.png"),
            dbk_rda_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_DOUBLETSKEPT_NON-NORMALIZED.rda"),
            dbk_technical_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/LogNormalize/pca/dims15_res0.8/technical/" + "{sample_name_ge}_technical_MULTI_ALL_uMAPs.png"),
            dbf_QC_hist_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_QChist.png"),
            dbf_rda_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_FILTERED_NON-NORMALIZED.rda"),
            dbf_stat_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_stat.txt")
        params:
            sing_bind = filtering_params_sing,
            pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
            input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
            output_folder = os.path.normpath("/WORKDIR/" + "{output_filtering_dir_ge}") + "/",
            SING_FILERING_CC_SEURAT_FILE = os.path.normpath("/WORKDIR/" + FILERING_CC_SEURAT_FILE) if FILERING_CC_SEURAT_FILE != "NULL" else "NULL",
            SING_FILERING_CC_CYCLONE_FILE = os.path.normpath("/WORKDIR/" + FILERING_CC_CYCLONE_FILE) if FILERING_CC_CYCLONE_FILE != "NULL" else "NULL",
            SING_FILERING_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in FILERING_METADATA_FILE.split(',')]) if FILERING_METADATA_FILE != "NULL" else "NULL"
        threads:
            4
        resources:
            mem_mb = (lambda wildcards, attempt: min(5121 + attempt * 5121, 51200)),
            time_min = (lambda wildcards, attempt: min(attempt * 180, 1000))
        shell:
            """
            export TMPDIR={GLOBAL_TMP}
            TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
            singularity exec --no-home -B $TMP_DIR:/tmp {params.sing_bind} \
            {SINGULARITY_ENV} \
            Rscript {params.pipeline_folder}/scripts/pipeline_part2.R \
            --input.rda.ge {params.input_rda} \
            --output.dir.ge {params.output_folder} \
            --author.name {FILERING_AUTHOR_NAME} \
            --author.mail {FILERING_AUTHOR_MAIL} \
            --nthreads {threads} \
            --pipeline.path {params.pipeline_folder} \
            --pcmito.min {FILERING_PCMITO_MIN} \
            --pcmito.max {FILERING_PCMITO_MAX} \
            --pcribo.min {FILERING_PCRIBO_MIN} \
            --pcribo.max {FILERING_PC_RIBO_MAX} \
            --min.features {FILERING_MIN_FEATURES} \
            --min.counts {FILERING_MIN_COUNTS} \
            --min.cells {FILERING_MIN_CELLS} \
            --doublets.filter.method {FILERING_DOUBLET_FILTER_METHOD} \
            --cc.seurat.file {params.SING_FILERING_CC_SEURAT_FILE} \
            --cc.cyclone.file {params.SING_FILERING_CC_CYCLONE_FILE} \
            --metadata.file {params.SING_FILERING_METADATA_FILE} && \
            rm -r $TMP_DIR || rm -r $TMP_DIR
            """

"""
This rule launches the R script to filter droplets according to quality control but NOT to filter the doublets.
"""
if FILERING_DOUBLET_FILTER_METHOD_NAME == "none":
    rule filtering_ge_none:
        input:
            rda_file = filtering_input_ge
        output:
            dbk_QC_hist_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_QChist.png"),
            dbk_rda_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_DOUBLETSKEPT_NON-NORMALIZED.rda"),
            dbk_technical_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/LogNormalize/pca/pca_15_0.8/technical/" + "{sample_name_ge}_technical_MULTI_ALL_uMAPs.png"),
            dbk_stat_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "/{sample_name_ge}_stat.txt")
        params:
            sing_bind = filtering_params_sing,
            pipeline_folder = os.path.normpath("/WORKDIR/" + PIPELINE_FOLDER),
            input_rda = lambda wildcards, input: os.path.normpath("/WORKDIR/" + input[0]),
            output_folder = os.path.normpath("/WORKDIR/" + "{output_filtering_dir_ge}") + "/",
            SING_FILERING_CC_SEURAT_FILE = os.path.normpath("/WORKDIR/" + FILERING_CC_SEURAT_FILE) if FILERING_CC_SEURAT_FILE != "NULL" else "NULL",
            SING_FILERING_CC_CYCLONE_FILE = os.path.normpath("/WORKDIR/" + FILERING_CC_CYCLONE_FILE) if FILERING_CC_CYCLONE_FILE != "NULL" else "NULL",
            SING_FILERING_METADATA_FILE = ','.join([os.path.normpath("/WORKDIR/" + x) for x in FILERING_METADATA_FILE.split(',')]) if FILERING_METADATA_FILE != "NULL" else "NULL"
        threads:
            4
        resources:
            mem_mb = (lambda wildcards, attempt: min(5121 + attempt * 5121, 51200)),
            time_min = (lambda wildcards, attempt: min(attempt * 180, 1000))
        shell:
            """
            export TMPDIR={GLOBAL_TMP}
            TMP_DIR=$(mktemp -d -t sc_pipeline-XXXXXXXXXX) && \
            singularity exec --no-home -B $TMP_DIR:/tmp {params.sing_bind} \
            {SINGULARITY_ENV} \
            Rscript {params.pipeline_folder}/scripts/pipeline_part2.R \
            --input.rda.ge {params.input_rda} \
            --output.dir.ge {params.output_folder} \
            --author.name {FILERING_AUTHOR_NAME} \
            --author.mail {FILERING_AUTHOR_MAIL} \
            --nthreads {threads} \
            --pipeline.path {params.pipeline_folder} \
            --pcmito.min {FILERING_PCMITO_MIN} \
            --pcmito.max {FILERING_PCMITO_MAX} \
            --pcribo.min {FILERING_PCRIBO_MIN} \
            --pcribo.max {FILERING_PC_RIBO_MAX} \
            --min.features {FILERING_MIN_FEATURES} \
            --min.counts {FILERING_MIN_COUNTS} \
            --min.cells {FILERING_MIN_CELLS} \
            --doublets.filter.method {FILERING_DOUBLET_FILTER_METHOD} \
            --cc.seurat.file {params.SING_FILERING_CC_SEURAT_FILE} \
            --cc.cyclone.file {params.SING_FILERING_CC_CYCLONE_FILE} \
            --metadata.file {params.SING_FILERING_METADATA_FILE} && \
            rm -r $TMP_DIR || rm -r $TMP_DIR
            """
