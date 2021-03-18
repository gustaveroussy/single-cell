"""
##########################################################################
This function make the input of rule all
##########################################################################
"""
def get_targets(STEPS):
    targets = {}
    if "Alignment_countTable_GE" in STEPS:
        targets["Alignment_countTable_GE"]=[
        #multiqc
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/QC/multiqc/{sample_name_ge}_RAW.html"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/QC/multiqc/{sample_name_ge}_RAW_data.zip"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        #alignment
        # expand(os.path.join(OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/output.bus"), zip, sample_name_ge=SAMPLE_NAME_GE),
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/transcripts.txt"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/matrix.ec"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/run_info.json"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        #sort
        # expand(os.path.join(OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}_sorted.bus"), zip, sample_name_ge=SAMPLE_NAME_GE),
        #correct_UMIs
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}_corrected.bus"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        #build_count_matrix
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.mtx"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.barcodes.txt"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/{sample_name_ge}.genes.txt"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE),
        expand(os.path.join(ALIGN_OUTPUT_DIR_GE,"{sample_name_ge}/KALLISTOBUS/Materials_and_Methods.txt"), zip, sample_name_ge=ALIGN_SAMPLE_NAME_GE)
        ]
    if "Alignment_countTable_ADT" in STEPS:
        targets["Alignment_countTable_ADT"]=[
        #multiqc
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/QC/multiqc/{sample_name_adt}_RAW.html"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/QC/multiqc/{sample_name_adt}_RAW_data.zip"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        #alignment
        # expand(os.path.join(OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/output.bus"), zip, sample_name_adt=SAMPLE_NAME_ADT),
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/transcripts.txt"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/matrix.ec"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/run_info.json"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        #sort
        # expand(os.path.join(OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/{sample_name_adt}_sorted.bus"), zip, sample_name_adt=SAMPLE_NAME_ADT),
        #correct_UMIs
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/{sample_name_adt}_corrected.bus"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        #build_count_matrix
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.mtx"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.barcodes.txt"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/{sample_name_adt}.genes.txt"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT),
        expand(os.path.join(ALIGN_OUTPUT_DIR_ADT,"{sample_name_adt}/KALLISTOBUS/Materials_and_Methods.txt"), zip, sample_name_adt=ALIGN_SAMPLE_NAME_ADT)
        ]
    if "Alignment_annotations_TCR_BCR" in STEPS:
        targets["Alignment_annotations_TCR_BCR"]=[
        #multiqc
        expand(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/QC/multiqc/{sample_name_tcr_bcr}_RAW.html"), zip, sample_name_tcr_bcr=ALIGN_SAMPLE_NAME_TCR_BCR),
        expand(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/QC/multiqc/{sample_name_tcr_bcr}_RAW_data.zip"), zip, sample_name_tcr_bcr=ALIGN_SAMPLE_NAME_TCR_BCR),
        #alignment_annotations_
        expand(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/filtered_contig_annotations.csv"), zip, sample_name_tcr_bcr=ALIGN_SAMPLE_NAME_TCR_BCR),
        expand(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/web_summary.html"), zip, sample_name_tcr_bcr=ALIGN_SAMPLE_NAME_TCR_BCR),
        expand(os.path.join(ALIGN_OUTPUT_DIR_TCR_BCR,"{sample_name_tcr_bcr}/{sample_name_tcr_bcr}_CellRanger/outs/Materials_and_Methods.txt"), zip, sample_name_tcr_bcr=ALIGN_SAMPLE_NAME_TCR_BCR)
        ]
    if "Droplets_QC_GE" in STEPS:
        targets["Droplets_QC_GE"]=[
        expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/UNFILTERED/" + "{sample_name_ge}_kneeplot.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE) if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/UNFILTERED_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_kneeplot.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE),
        expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/UNFILTERED/" + "{sample_name_ge}_saturation_plot.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE) if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/UNFILTERED_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_saturation_plot.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE),
        expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/UNFILTERED/" + "{sample_name_ge}_QChist.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE) if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/UNFILTERED_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_QChist.png"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE),
        expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/UNFILTERED/" + "{sample_name_ge}_UNFILTERED_NON-NORMALIZED.rda"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE) if  str(QC_EMPTYDROPS_RETAIN) == "NULL" else expand(os.path.normpath("{outputqc_droplets_dir_ge}" + "/UNFILTERED_retain" + str(QC_EMPTYDROPS_RETAIN) + "/{sample_name_ge}_UNFILTERED_NON-NORMALIZED.rda"), zip, outputqc_droplets_dir_ge=QC_OUTPUT_DIR_GE, sample_name_ge=QC_SAMPLE_NAME_GE)
        ]
    if "Filtering_GE" in STEPS:
        targets["Filtering_GE"]=[
            expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_QChist.png"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE),
            expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "{sample_name_ge}_DOUBLETSKEPT_NON-NORMALIZED.rda"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE),
            expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/LogNormalize/pca/pca_15_0.8/technical/" + "{sample_name_ge}_technical_MULTI_ALL_uMAPs.png"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE)
        ]
        if FILERING_DOUBLET_FILTER_METHOD_NAME == "none":
            targets["Filtering_GE"].append(expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSKEPT/" + "/{sample_name_ge}_stat.txt"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE))
        else:
            targets["Filtering_GE"].append(expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_QChist.png"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE))
            targets["Filtering_GE"].append(expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_FILTERED_NON-NORMALIZED.rda"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE))
            targets["Filtering_GE"].append(expand(os.path.normpath("{output_filtering_dir_ge}/" + FILTERS_FOLDER + "/DOUBLETSFILTER_" + FILERING_DOUBLET_FILTER_METHOD_NAME + "/{sample_name_ge}_stat.txt"), zip, output_filtering_dir_ge=FILERING_OUTPUT_DIR_GE, sample_name_ge=FILERING_SAMPLE_NAME_GE))
    if "Norm_DimRed_Eval_GE" in STEPS:
        targets["Norm_DimRed_Eval_GE"]=[
        expand(os.path.normpath("{output_norm_dimred_dir_ge}" + "/" + NDRE_NORM_VTR + "/" + NDRE_DIMRED_VTR + "/" + "{sample_name_ge}_" + NDRE_NORM_VTR + "_" + NDRE_DIMRED_VTR + ".rda"), zip, output_norm_dimred_dir_ge=NDRE_OUTPUT_DIR_GE, sample_name_ge=NDRE_SAMPLE_NAME_GE),
        expand(os.path.normpath("{output_norm_dimred_dir_ge}" + "/" + NDRE_NORM_VTR + "/" + NDRE_DIMRED_VTR + "/" + "{sample_name_ge}_" + ASSAY + "_" + NDRE_DIMRED_METHOD + "_dims.bias.cor.png"), zip, output_norm_dimred_dir_ge=NDRE_OUTPUT_DIR_GE, sample_name_ge=NDRE_SAMPLE_NAME_GE),
        expand(os.path.normpath("{output_norm_dimred_dir_ge}" + "/" + NDRE_NORM_VTR + "/" + NDRE_DIMRED_VTR + "/clustree_" + ASSAY + "_" + NDRE_DIMRED_METHOD + "/dimensions/{sample_name_ge}_" + ASSAY + "_" + NDRE_DIMRED_METHOD + "{dims}.png"), zip, output_norm_dimred_dir_ge=NDRE_OUTPUT_DIR_GE*len(POSSIBLE_DIM), sample_name_ge=NDRE_SAMPLE_NAME_GE*len(POSSIBLE_DIM), dims=POSSIBLE_DIM*len(NDRE_SAMPLE_NAME_GE)),
        expand(os.path.normpath("{output_norm_dimred_dir_ge}" + "/" + NDRE_NORM_VTR + "/" + NDRE_DIMRED_VTR + "/clustree_" + ASSAY + "_" + NDRE_DIMRED_METHOD + "/louvain_resolution/{sample_name_ge}_" + ASSAY + "_res{res}.png"), zip, output_norm_dimred_dir_ge=NDRE_OUTPUT_DIR_GE*len(POSSIBLE_DIM), sample_name_ge=NDRE_SAMPLE_NAME_GE*len(POSSIBLE_DIM), res=POSSIBLE_RES*len(NDRE_SAMPLE_NAME_GE)),
        expand(os.path.normpath("{output_norm_dimred_dir_ge}" + "/" + NDRE_NORM_VTR + "/" + NDRE_DIMRED_VTR + "/clustree_" + ASSAY + "_" + NDRE_DIMRED_METHOD + "/uMAPs/{sample_name_ge}_uMAPs_" + ASSAY + "_" + NDRE_DIMRED_METHOD + "{dims}_ALLres.png"), zip, output_norm_dimred_dir_ge=NDRE_OUTPUT_DIR_GE*len(POSSIBLE_DIM), sample_name_ge=NDRE_SAMPLE_NAME_GE*len(POSSIBLE_DIM), dims=POSSIBLE_DIM*len(NDRE_SAMPLE_NAME_GE))
        ]
    if "Clust_Markers_Annot_GE" in STEPS:
        targets["Clust_Markers_Annot_GE"]=[
        expand(os.path.normpath("{output_clust_markers_annot_dir_ge}" + "/" + CMA_CLUST_FOLDER + "/" + "{sample_name_ge}" + "{complement}" + "_" + str(CMA_KEEP_DIM) + "_" + str(CMA_KEEP_RES) + ".rda"), zip, output_clust_markers_annot_dir_ge = CMA_OUTPUT_DIR_GE, sample_name_ge = CMA_SAMPLE_NAME_GE, complement = CMA_COMPLEMENT)
        ]
    if "Adding_ADT" in STEPS:
        targets["Adding_ADT"]=[
        expand("{add_adt_output}_ADT.rda", add_adt_output = ADD_ADT_OUTPUT)
        ]
    if "Adding_TCR" in STEPS:
        targets["Adding_TCR"]=[
        expand("{add_tcr_output}_TCR.rda", add_tcr_output = ADD_TCR_OUTPUT)
        ]
    if "Adding_BCR" in STEPS:
        targets["Adding_BCT"]=[
        expand("{add_bcr_output}_BCR.rda", add_bcr_output = ADD_BCR_OUTPUT)
        ]
    if "Cerebro" in STEPS:
        targets["Cerebro"]=[
        expand("{cerebro_input_rda_no_extention}{cerebro_complement}", cerebro_input_rda_no_extention = CEREBRO_INPUT_RDA_NO_EXTENTION, cerebro_complement = CEREBRO_COMPLEMENT_CRB)
        ]
    return targets
