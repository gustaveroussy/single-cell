"""
##########################################################################
These rules filter droplets according to quality-control and whether to filter doublets or not, in single-cell RNA-seq.
##########################################################################
"""

"""
This functions allows to determine the input .rda file.
"""
def filtering_input_ge(wildcards):
    for key, value in dic_FILTER_INFO.items():
        if value["FILTER_OUTPUT_DIR"] == wildcards.output_filtering_dir_ge :
            return key
            #key is the input value in the yaml parameter file.

"""
This function allows to determine the singularity binding parameters.
"""
def filtering_params_sing(wildcards):
    for key, value in dic_FILTER_INFO.items():
        if value["FILTER_OUTPUT_DIR"] == wildcards.output_filtering_dir_ge :
            rda_folder = os.path.dirname(key)
    output_folder = wildcards.output_filtering_dir_ge + "/"
    concat = " -B " + PIPELINE_FOLDER + "," + rda_folder + "," + output_folder 
    if os.path.exists(FILTERING_CC_SEURAT_FILE): concat = concat + "," + os.path.dirname(FILTERING_CC_SEURAT_FILE)
    if os.path.exists(FILTERING_CC_CYCLONE_FILE): concat = concat + "," + os.path.dirname(FILTERING_CC_CYCLONE_FILE)
    if os.path.exists(FILTERING_METADATA_FILE):
        for metadatafile in list(dict.fromkeys(FILTERING_METADATA_FILE.split(","))):
            metadatafile = os.path.dirname(metadatafile)
            concat = concat + "," + metadatafile
    return concat

"""
This rule launches the R script to filter droplets according to quality control and to filter the doublets.
"""
if FILTERING_DOUBLET_FILTER_METHOD_NAME != "none":
    rule filtering_ge:
        input:
            rda_file = filtering_input_ge
        output:
            dbk_QC_hist_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_QChist.png"),
            dbk_rda_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_DOUBLETSKEPT_NON-NORMALIZED.rda"),
            dbk_technical_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/LogNormalize/pca/dims15_res0.8/technical/" + "{sample_name_ge}_technical_MULTI_ALL_uMAPs.png"),
            dbf_QC_hist_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILTERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_QChist.png"),
            dbf_rda_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILTERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_FILTERED_NON-NORMALIZED.rda"),
            dbf_stat_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILTERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_stat.csv")
        params:
            sing_bind = filtering_params_sing
        log:
            "logs/filtering_ge{output_filtering_dir_ge}{sample_name_ge}.log"
        benchmark:
            "benchmark/filtering_ge{output_filtering_dir_ge}{sample_name_ge}.tsv"
        threads:
            4
        resources:
            mem_mb = (lambda wildcards, attempt: min(5121 + attempt * 5121, 51200)),
            time_min = (lambda wildcards, attempt: min(attempt * 180, 1000))
        shell:
            """
            TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
            trap "rm -r $TMPDIR" EXIT
            singularity exec --no-home -B $TMPDIR:/tmp {params.sing_bind} \
            {SINGULARITY_ENV} \
            Rscript {PIPELINE_FOLDER}/scripts/pipeline_part2.R \
            --input.rda.ge {input[0]} \
            --output.dir.ge {wildcards.output_filtering_dir_ge}/ \
            --author.name "{AUTHOR_NAME}" \
            --author.mail "{AUTHOR_MAIL}" \
            --nthreads {threads} \
            --pipeline.path {PIPELINE_FOLDER} \
            --pcmito.min {FILTERING_PCMITO_MIN} \
            --pcmito.max {FILTERING_PCMITO_MAX} \
            --pcribo.min {FILTERING_PCRIBO_MIN} \
            --pcribo.max {FILTERING_PC_RIBO_MAX} \
            --min.features {FILTERING_MIN_FEATURES} \
            --min.counts {FILTERING_MIN_COUNTS} \
            --min.cells {FILTERING_MIN_CELLS} \
            --doublets.filter.method {FILTERING_DOUBLET_FILTER_METHOD} \
            --cc.seurat.file {FILTERING_CC_SEURAT_FILE} \
            --cc.cyclone.file {FILTERING_CC_CYCLONE_FILE} \
            --metadata.file {FILTERING_METADATA_FILE} &> {log}
            """

"""
This rule launches the R script to filter droplets according to quality control but NOT to filter the doublets.
"""
if FILTERING_DOUBLET_FILTER_METHOD_NAME == "none":
    rule filtering_ge_none:
        input:
            rda_file = filtering_input_ge
        output:
            dbk_QC_hist_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_QChist.png"),
            dbk_rda_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_DOUBLETSKEPT_NON-NORMALIZED.rda"),
            dbk_technical_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/LogNormalize/pca/dims15_res0.8/technical/" + "{sample_name_ge}_technical_MULTI_ALL_uMAPs.png"),
            dbk_stat_file = os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "/{sample_name_ge}_stat.csv")
        params:
            sing_bind = filtering_params_sing,
        log:
            "logs/filtering_ge/{output_filtering_dir_ge}{sample_name_ge}.log"
        benchmark:
            "benchmark/filtering_ge/{output_filtering_dir_ge}{sample_name_ge}.tsv"
        threads:
            4
        resources:
            mem_mb = (lambda wildcards, attempt: min(5121 + attempt * 5121, 51200)),
            time_min = (lambda wildcards, attempt: min(attempt * 180, 1000))
        shell:
            """
            TMPDIR=$(mktemp -d {resources.tmpdir}/XXXXXX)
            trap "rm -r $TMPDIR" EXIT
            singularity exec --no-home -B $TMPDIR:/tmp {params.sing_bind} \
            {SINGULARITY_ENV} \
            Rscript {PIPELINE_FOLDER}/scripts/pipeline_part2.R \
            --input.rda.ge {input[0]} \
            --output.dir.ge {wildcards.output_filtering_dir_ge}/ \
            --author.name "{AUTHOR_NAME}" \
            --author.mail "{AUTHOR_MAIL}" \
            --nthreads {threads} \
            --pipeline.path {PIPELINE_FOLDER} \
            --pcmito.min {FILTERING_PCMITO_MIN} \
            --pcmito.max {FILTERING_PCMITO_MAX} \
            --pcribo.min {FILTERING_PCRIBO_MIN} \
            --pcribo.max {FILTERING_PC_RIBO_MAX} \
            --min.features {FILTERING_MIN_FEATURES} \
            --min.counts {FILTERING_MIN_COUNTS} \
            --min.cells {FILTERING_MIN_CELLS} \
            --doublets.filter.method {FILTERING_DOUBLET_FILTER_METHOD} \
            --cc.seurat.file {FILTERING_CC_SEURAT_FILE} \
            --cc.cyclone.file {FILTERING_CC_CYCLONE_FILE} \
            --metadata.file {FILTERING_METADATA_FILE} &> {log}
            """
