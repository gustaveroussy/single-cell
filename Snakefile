from pathlib import Path
import snakemake.utils
import os
import datetime
import re
import glob
import zipfile
import json
import numpy

__author__ = "Marine AGLAVE & Bastien JOB"

#using: snakemake --profile /mnt/beegfs02/pipelines/bigr_single-cell/<version>/profiles/slurm -s /mnt/beegfs02/pipelines/bigr_single-cell/<version>/Snakefile --configfile /path/to/Params.yaml

pipeline_version = "2.0.0"
gitub_pipeline = "https://github.com/gustaveroussy/single-cell"

sys.stderr.write("\n############################################################# \n")
sys.stderr.write("\n\n\t Single-cell RNA-seq pipeline (version " + pipeline_version + ")\n\n")
sys.stderr.write("\n############################################################# \n\n")
sys.stderr.write("\nSee " + gitub_pipeline + " for more information about the pipeline. \n\n")

PIPELINE_FOLDER = workflow.snakefile
PIPELINE_FOLDER = PIPELINE_FOLDER.replace("/Snakefile", "")

include: "rules/Check_config_and_Set_variables.smk"

### rule all ###################################################################################################################################
sys.stderr.write("\n########################### Run ############################\n\n")

include: "rules/Rule_all.smk"
rule all:
    input:
        **get_targets(STEPS)
    message:
        "Single-cell RNA-seq pipeline done!"


### other rules ###################################################################################################################################
if "Alignment_countTable_GE" in STEPS:
    include: "rules/Alignment_countTable_GE.smk"
    wildcard_constraints:
        sample_name_ge = ".+_GE"

if "Alignment_countTable_ADT" in STEPS:
    include: "rules/Alignment_countTable_ADT.smk"
    wildcard_constraints:
        sample_name_adt = ".+_ADT"

if "Alignment_annotations_TCR_BCR" in STEPS:
    include: "rules/Alignment_annotations_TCR_BCR.smk"
    wildcard_constraints:
        sample_name_tcr_bcr=".+(_TCR|_BCR)"

if "Droplets_QC_GE" in STEPS:
    include: "rules/Droplets_QC_GE.smk"
    wildcard_constraints:
        sample_name_ge = ".+_GE"

if "Filtering_GE" in STEPS:
    include: "rules/Filtering_GE.smk"
    wildcard_constraints:
        sample_name_ge = ".+_GE"

if "Norm_DimRed_Eval_GE" in STEPS:
    include: "rules/Norm_DimRed_Eval_GE.smk"
    wildcard_constraints:
        sample_name_ge = ".+_GE",
        ndre_norm_vtr='|'.join([x for x in NDRE_NORM_VTR]),
        ndre_dimred_vtr='|'.join([x for x in NDRE_DIMRED_VTR])

if "Clust_Markers_Annot_GE" in STEPS:
    include: "rules/Clust_Markers_Annot_GE.smk"
    wildcard_constraints:
        sample_name_ge = ".+_GE",
        keep_dim_res='|'.join([x for x in CMA_KEEP_DIM_RES]),
        clust_folders='|'.join([x for x in CMA_CLUST_FOLDERS])

if "Adding_ADT" in STEPS:
    include: "rules/Adding_ADT.smk"

if "Adding_TCR" in STEPS:
    include: "rules/Adding_TCR.smk"

if "Adding_BCR" in STEPS:
    include: "rules/Adding_BCR.smk"

if "Cerebro" in STEPS:
    include: "rules/Cerebro.smk"

if "Int_Norm_DimRed_Eval_GE" in STEPS:
    include: "rules/Int_Norm_DimRed_Eval_GE.smk"
    wildcard_constraints:
        name_int = "|".join(INT_NDRE_NAME_INT),
        int_ndre_norm_vtr='|'.join([x for x in INT_NDRE_NORM_VTR]),
        int_ndre_dimred_vtr='|'.join([x for x in INT_NDRE_DIMRED_VTR])

if "Int_Clust_Markers_Annot_GE" in STEPS:
    include: "rules/Int_Clust_Markers_Annot_GE.smk"
    wildcard_constraints:
        name_int = "|".join(INT_CMA_NAME_INT),
        int_keep_dim_res='|'.join([x for x in INT_CMA_KEEP_DIM_RES]),
        int_clust_folders='|'.join([x for x in INT_CMA_CLUST_FOLDERS])

if "Int_Adding_ADT" in STEPS:
    include: "rules/Int_Adding_ADT.smk"
    wildcard_constraints:
        int_add_adt_output = "|".join(INT_ADD_ADT_OUTPUT)

if "Int_Adding_TCR" in STEPS:
    include: "rules/Int_Adding_TCR.smk"
    wildcard_constraints:
        int_add_tcr_output = "|".join(INT_ADD_TCR_OUTPUT)

if "Int_Adding_BCR" in STEPS:
    include: "rules/Int_Adding_BCR.smk"
    wildcard_constraints:
        int_add_bcr_output = "|".join(INT_ADD_BCR_OUTPUT)

if "Grp_Norm_DimRed_Eval_GE" in STEPS:
    include: "rules/Grp_Norm_DimRed_Eval_GE.smk"
    wildcard_constraints:
        name_grp = "|".join(GRP_NDRE_NAME_GRP),
        grp_ndre_norm_vtr='|'.join([x for x in GRP_NDRE_NORM_VTR]),
        grp_ndre_dimred_vtr='|'.join([x for x in GRP_NDRE_DIMRED_VTR])

if "Grp_Clust_Markers_Annot_GE" in STEPS:
    include: "rules/Grp_Clust_Markers_Annot_GE.smk"
    wildcard_constraints:
        name_grp = "|".join(GRP_CMA_NAME_GRP),
        grp_keep_dim_res='|'.join([x for x in GRP_CMA_KEEP_DIM_RES]),
        grp_clust_folders='|'.join([x for x in GRP_CMA_CLUST_FOLDERS])

if "Grp_Adding_ADT" in STEPS:
    include: "rules/Grp_Adding_ADT.smk"
    wildcard_constraints:
        grp_add_adt_output = "|".join(GRP_ADD_ADT_OUTPUT)

if "Grp_Adding_TCR" in STEPS:
    include: "rules/Grp_Adding_TCR.smk"
    wildcard_constraints:
        grp_add_tcr_output = "|".join(GRP_ADD_TCR_OUTPUT)

if "Grp_Adding_BCR" in STEPS:
    include: "rules/Grp_Adding_BCR.smk"
    wildcard_constraints:
        grp_add_bcr_output = "|".join(GRP_ADD_BCR_OUTPUT)

