#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda", help="Input seurat object (in .rda format)."),
  make_option("--output.dir", help="Output path"),
  make_option("--vdj.input.file.tcr", help="File filtered_contig_annotations.csv from CellRanger aligment pipeline."),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Yaml parameters file to remplace all parameters before (usefull to use R script without snakemake)
  make_option("--yaml", help="Patho to yaml file with all parameters")
)
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)

#### Formatting Parameters ####
#convert "NULL"/"FALSE"/"TRUE" (in character) into NULL/FALSE/TRUE
for (i in names(args$options)){
  if (toupper(args$options[i]) == "NULL") { args$options[i] <- NULL
  } else if (toupper(args$options[i]) == "FALSE") { args$options[i] <- FALSE
  } else if (toupper(args$options[i]) == "TRUE") { args$options[i] <- TRUE
  }
}

#### Get Paramaters ####
### Project
input.rda <- args$options$input.rda
output.dir <- args$options$output.dir
vdj.input.file.tcr <- args$options$vdj.input.file.tcr
author.name <- args$options$author.name
author.mail <- args$options$author.mail
### Computational Parameters
pipeline.path <- args$options$pipeline.path
### Yaml parameters file to remplace all parameters before (usefull to use R script without snakemake)
if (!is.null(args$options$yaml)){
  yaml_options <- yaml::yaml.load_file(args$options$yaml)
  for(i in names(yaml_options)) assign(i, yaml_options[[i]])
  rm(yaml_options)
}
### Clean
rm(args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(input.rda)) stop("input.rda parameter can't be empty!")
if (is.null(output.dir)) stop("output.dir parameter can't be empty!")
if (is.null(vdj.input.file.tcr)) stop("vdj.input.file.tcr parameter can't be empty!")

### Load data
load(input.rda)

### Save project parameters
if (!is.null(author.name) && !tolower(author.name) %in% tolower(sobj@misc$params$author.name)) sobj@misc$params$author.name <- c(sobj@misc$params$author.name, author.name)
if (!is.null(author.mail) && !tolower(author.mail) %in% tolower(sobj@misc$params$author.mail)) sobj@misc$params$author.mail <- c(sobj@misc$params$author.mail, author.mail)

#### Get Missing Paramaters ####
### Analysis Parameters
## Normalization and dimension reduction
dimred.method <- sobj@misc$params$reductions$method
## Clustering
GE_file <- sub('.rda', '', input.rda)
dimred.method <- sobj@misc$params$reductions$method
ident.name <- sobj@misc$params$clustering$ident
RNA.reduction <- sobj@misc$params$clustering$umap
sample.name <- sub("_GE", "", sobj@misc$params$sample.name.GE)

#### Fixed parameters ####
output_path_TCR <- paste0(output.dir, "/TCR_results/")
list_type_clT <- c("gene+nt", "gene", "nt", "aa")
list_type_contig <- c("nt","aa")
caption <- '“gene” - use the genes comprising the TCR
“nt” - use the nucleotide sequence of the CDR3 region
“aa” - use the amino acid sequence of the CDR3 region
“gene+nt” - use the genes comprising the TCR + the nucleotide sequence of the CDR3 region for T cells. This is the proper definition of clonotype.'

#### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

#########
## MAIN
#########

### printing parameters:
print("###########################################")
print(paste0("input.rda : ",input.rda))
print(paste0("ident.name : ",ident.name))
print(paste0("vdj.input.file.tcr : ",vdj.input.file.tcr))
print("###########################################")

## Load libraries
require(patchwork)
suppressMessages(require(Seurat))
library(dplyr)

## Set the seed
set.seed(sobj@misc$params$seed)

## GLOBAL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
global_output <- paste0(output_path_TCR, "Global_analysis")
dir.create(path = global_output, recursive = TRUE, showWarnings = TRUE)

## Loading input data and Combining contigs
cr_res <- load.sc.tcr.bcr(sobj=sobj, vdj.input.file = vdj.input.file.tcr)
tcr.combined <- scRepertoire::combineTCR(df = list(cr_res), samples = sample.name, ID = "TCR", cells = "T-AB")

## Quantification analysis
QC.tcr.bcr(cr_res=cr_res, out.dir=global_output, type="TCR")

## Quantification of unique contig analysis
Quantif.unique.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name)

## Abundance analysis
### Plots
for(x in list_type_clT) assign(paste0("plot_abundanceContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::abundanceContig(tcr.combined, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
### Save
png(paste0(global_output,'/abundanceContig.png'), width = 800, height = 300)
(plot_abundanceContig_gene_nt | plot_abundanceContig_gene | plot_abundanceContig_nt | plot_abundanceContig_aa ) +
  plot_annotation(title = sample.name, subtitle = paste0("(",dim(sobj)[2]," cells)"), caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
dev.off()

## Contigs Length analysis
for(x in list_type_contig){
  ### This should give multimodal plot
  assign(paste0("plot_lengthContig_",x,"_comb"), scRepertoire::lengthContig(tcr.combined, cloneCall=x, chains = "combined") + Seurat::NoLegend())
  ### Plots the A and B chains distribution separately
  assign(paste0("plot_lengthContig_",x,"_sin"), scRepertoire::lengthContig(tcr.combined, cloneCall=x, chains = "single") + Seurat::NoLegend())
}
### Save
png(paste0(global_output,'/lengthContig.png'), width = 800, height = 800)
(plot_lengthContig_nt_comb | plot_lengthContig_nt_sin) / (plot_lengthContig_aa_comb | plot_lengthContig_aa_sin ) +
  plot_annotation(title = sample.name, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
dev.off()

## Clonal Homeostasis analysis
Homeo.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name)

## Clonal Proportions analysis
Prop.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name)

## Diversity analysis
Div.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name)

## Combine TCR data with seurat object
### Corresponding barcode
tcr.combined[[1]]$barcode <- gsub(pattern = paste0(sample.name, "_TCR_"), replacement = '', tcr.combined[[1]]$barcode)
### Combination
sobj <- scRepertoire::combineSeurat(df = tcr.combined, seurat = sobj, cloneCall="aa")

### Spliting CTstrict (into separate columns for TRA-V/J/C, TRB-V/J/C and corresponding sequences, with 2 possible clonotypes) and save as metadata
### and Adding length of TR sequence to meta.data
sobj <- split.CTstrict.tcr(sobj)

### Plots
for(x in c("TRAV_1","TRAJ_1","TRAC_1","TRBV_1","TRBJ_1","TRBC_1","TRAV_2","TRAJ_2","TRAC_2","TRBV_2","TRBJ_2","TRBC_2")) assign(paste0("dimplot_",x), if (all(is.na(sobj@meta.data[[x]]))) patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x) else Seurat::DimPlot(sobj, group.by = x, reduction = RNA.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x))
for(x in c("TRA_nt_1_len","TRA_nt_2_len","TRB_nt_1_len","TRB_nt_2_len")) {
  assign(paste0("featureplot_",x),
         tryCatch({ Seurat::FeaturePlot(sobj, features = x, reduction = RNA.reduction) },
            error=function(e) { return(patchwork::plot_spacer()) }
         ) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x)))
}
### Save plots
png(paste0(global_output,'/cloneType.png'), width =1400, height = 3000)
( wrap_elements( (dimplot_TRAV_1 / dimplot_TRAJ_1 / dimplot_TRAC_1 / dimplot_TRBV_1 / dimplot_TRBJ_1 / dimplot_TRBC_1 / featureplot_TRA_nt_1_len / featureplot_TRB_nt_1_len) +
     plot_annotation(title = 'TR clonotype 1', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) |
    wrap_elements( (dimplot_TRAV_2 / dimplot_TRAJ_2 /dimplot_TRAC_2 / dimplot_TRBV_2 / dimplot_TRBJ_2 / dimplot_TRBC_2 / featureplot_TRA_nt_2_len / featureplot_TRB_nt_2_len) +
     plot_annotation(title = 'TR clonotype 2', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) )
dev.off()

## Frequency analysis
sobj <- Freq.g(sobj=sobj, out.dir = global_output, sample.name=sample.name, reduction=RNA.reduction, freq_col="Frequency")

## Physicochemical properties of the CDR3
Physicochemical_properties.g(sobj=sobj, list_type_clT = list_type_clT, out.dir = global_output, sample.name=sample.name, type='TCR')




## CLUSTERS LEVEL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
clusters_output <- paste0(output_path_TCR, "Clusters_analysis")
dir.create(path = clusters_output, recursive = TRUE, showWarnings = TRUE)

## Filter cells that have no value for the x concerned + conversion to a list by clusters
## Need to filter, otherwise the functions count the 'NA' as a sequence.
for(x in list_type_clT){
  if(x=="gene+nt") y="strict" else y=x
  filtred_sobj = sobj[,!is.na(sobj@meta.data[paste0("CT", y)])]
  assign(paste0("filtred_metadata_", sub("\\+","_",x)), scRepertoire::seurat2List(filtred_sobj))
}
rm(filtred_sobj)

## Quantification of unique contig analysis
sobj <- Quantif.unique.c(sobj = sobj, ident.name=ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Abundance analysis
### Plots
for(x in list_type_clT) assign(paste0("plot_cluster_abundanceContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::abundanceContig(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
### Save
png(paste0(clusters_output,'/abundanceContig.png'), width = 2000, height = 600)
(plot_cluster_abundanceContig_gene_nt | plot_cluster_abundanceContig_gene | plot_cluster_abundanceContig_nt | plot_cluster_abundanceContig_aa ) +
  plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
dev.off()

## Clonal Homeostasis analysis
sobj <- Homeo.c(sobj = sobj, ident.name=ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Clonal Proportions analysis
sobj <- Prop.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Diversity analysis
sobj <- Div.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Frequency analysis
Freq.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, ident.name=ident.name, reduction=, freq_col="Frequency", filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Clonal Overlap analysis (si plus de 1)
if(length(levels(Seurat::Idents(sobj)))!=1) Overlap.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Physico-chemical properties of the CDR3
Physicochemical_properties.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, ident.name=ident.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt, type='TCR')

#renamme TCR columns with 'TCR_' prefix
toMatch <- c("^CTgene","^CTnt","^CTaa","^CTstrict","^Frequency","^cloneType","^TRAV_1","^TRAJ_1","^TRAC_1","^TRAV_2","^TRAJ_2","^TRAC_2","^TRA_nt_1","^TRA_nt_2","^TRBV_1","^TRBJ_1","^None_1","^TRBC_1","^TRBV_2","^TRBJ_2","^None_2","^TRBC_2","^TRB_nt_1","^TRB_nt_2","^TRA_nt_1_len","^TRA_nt_2_len","^TRB_nt_1_len","^TRB_nt_2_len","^highlight_aa")
matches <- grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data))
colnames(sobj@meta.data)[matches] <- paste0("TCR_", grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data), value=TRUE))

## Save packages versions
sobj@misc$technical_info$scRepertoire <- utils::packageVersion('scRepertoire')
sobj@misc$technical_info$alakazam <- utils::packageVersion('alakazam')

### Materials and Methods
if(file.exists(paste0(dirname(vdj.input.file.tcr), "/../../Materials_and_Methods.txt"))){
  tmp <- readr::read_tsv(paste0(dirname(vdj.input.file.tcr), "/../../Materials_and_Methods.txt"), col_names = FALSE)$X1
  tmp2 <- ""
  for (i in 1:length(tmp)) tmp2=paste(tmp2,tmp[i], sep="")
  sobj@misc$parameters$Materials_and_Methods$TCR <- tmp2
} else sobj@misc$parameters$Materials_and_Methods$TCR <- NULL
sobj@misc$parameters$Materials_and_Methods$TCR <- paste0(sobj@misc$parameters$Materials_and_Methods$TCR, " The annotation was merged with corresponding cell barcode of 5’ gene expression. The scRepertoire package (",sobj@misc$technical_info$scRepertoire,") was used to process annotation to assign clonotype based on TCR chains. scRepertoire allows to study contig quantification, contig abundance, contig length, clonal space homeostasis, clonal proportion, clonal overlap beetween clusters and diversity. Physicochemical properties of the CDR3, based on amino-acid sequences, was determined by the alakazam R package (",sobj@misc$technical_info$alakazam,").")
sobj@misc$parameters$Materials_and_Methods$packages_references <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)

### Saving GE_ADT_TCR object
GE_TCR_file <- paste0(output.dir, basename(GE_file), '_TCR')
save(sobj, file = paste0(GE_TCR_file, '.rda'), compress = "bzip2")
