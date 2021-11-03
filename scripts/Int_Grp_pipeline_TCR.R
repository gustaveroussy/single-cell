#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda", help="Input seurat object (in .rda format)."),
  make_option("--output.dir", help="Output path"),
  make_option("--vdj.input.files.tcr", help="Files filtered_contig_annotations.csv from CellRanger aligment pipeline (in the same order than seurat objects for integration or grouping)."),
  make_option("--author.name", help="Name of author of the analysis"),
  make_option("--author.mail", help="Email of author of the analysis"),
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
  if ((length(args$options[i]) == 0) || (length(args$options[i]) == 1 && toupper(args$options[i]) == "NULL")) { args$options[i] <- NULL
  } else if ((length(args$options[i]) == 1) && (toupper(args$options[i]) == "FALSE")) { args$options[i] <- FALSE
  } else if ((length(args$options[i]) == 1) && (toupper(args$options[i]) == "TRUE")) { args$options[i] <- TRUE
  }
}

#### Get Paramaters ####
### Project
input.rda <- args$options$input.rda
output.dir <- args$options$output.dir
vdj.input.files.tcr <- unlist(stringr::str_split(args$options$vdj.input.files.tcr, ","))
list.author.name <- if (!is.null(args$options$author.name)) unlist(stringr::str_split(args$options$author.name, ","))
list.author.mail <- if (!is.null(args$options$author.mail)) unlist(stringr::str_split(args$options$author.mail, ","))
### Computational Parameters
pipeline.path <- args$options$pipeline.path
### Yaml parameters file to remplace all parameters before (usefull to use R script without snakemake)
if (!is.null(args$options$yaml)){
  yaml_options <- yaml::yaml.load_file(args$options$yaml)
  for(i in names(yaml_options)) {
    #convert "NULL"/"FALSE"/"TRUE" (in character) into NULL/FALSE/TRUE
    if ((length(yaml_options[[i]]) == 0) || (length(yaml_options[[i]]) == 1 && toupper(yaml_options[[i]]) == "NULL")) { yaml_options[[i]] <- NULL
    } else if ((length(yaml_options[[i]]) == 1) && (toupper(yaml_options[[i]]) == "FALSE")) { yaml_options[[i]] <- FALSE
    } else if ((length(yaml_options[[i]]) == 1) && (toupper(yaml_options[[i]]) == "TRUE")) { yaml_options[[i]] <- TRUE
    }
    #assign values
    assign(i, yaml_options[[i]])
  }
  rm(yaml_options, i)
}
### Clean
rm(args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(input.rda)) stop("input.rda parameter can't be empty!")
if (is.null(output.dir)) stop("output.dir parameter can't be empty!")
if (is.null(vdj.input.files.tcr)) stop("vdj.input.files.tcr parameter can't be empty!")

### Load data
load(input.rda)

### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

### Save project parameters
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)

#### Get Missing Paramaters ####
### Analysis Parameters
## Normalization and dimension reduction
dimred.method <- sobj@misc$params$reductions$method
## Clustering
dimred.method <- sobj@misc$params$reductions$method
ident.name <- sobj@misc$params$clustering$ident
INT_GRP.reduction <- sobj@misc$params$clustering$umap
sample.name.INT_GRP <- Seurat::Project(sobj)
samples.name.GE <- sobj@misc$params$names.ge
samples.name <- sub("_GE","",samples.name.GE)

#### Fixed parameters ####
output_path_TCR <- paste0(output.dir, "/TCR_results/")
list_type_clT <- c("gene+nt", "gene", "nt", "aa")
list_type_contig <- c("nt","aa")
caption <- '“gene” - use the genes comprising the TCR
“nt” - use the nucleotide sequence of the CDR3 region
“aa” - use the amino acid sequence of the CDR3 region
“gene+nt” - use the genes comprising the TCR + the nucleotide sequence of the CDR3 region for T cells. This is the proper definition of clonotype.'

#########
## MAIN
#########

### printing parameters:
print("###########################################")
print(paste0("input.rda : ",input.rda))
print(paste0("ident.name : ",ident.name))
print(paste0("vdj.input.files.tcr : ",vdj.input.files.tcr))
print("###########################################")

## Load libraries
require(patchwork)
suppressMessages(require(Seurat))
library(dplyr)

## Set the seed
set.seed(sobj@misc$params$seed)

## GLOBAL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cat("\nGlobal Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
global_output <- paste0(output_path_TCR, "Global_analysis")
dir.create(path = global_output, recursive = TRUE, showWarnings = FALSE)

## Loading input data and Combining contigs
cat("\nLoading input data and Combining contigs...\n")
cr_res <- lapply(seq_along(vdj.input.files.tcr), load.sc.tcr.bcr, sobj=sobj, vdj.input.file=vdj.input.files.tcr, sample.name=samples.name.GE)
tcr.combined <- scRepertoire::combineTCR(df = cr_res, samples = samples.name, ID = rep("GE", length(samples.name)), cells = "T-AB")

## Quantification of unique contig analysis
cat("\nQuantification analysis...\n")
Quantif.unique.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name.INT_GRP)

## Abundance analysis
cat("\nAbundance analysis...\n")
### Plots
for(x in list_type_clT) assign(paste0("plot_abundanceContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::abundanceContig(tcr.combined, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
### Save
png(paste0(global_output,'/abundanceContig.png'), width = 800, height = 300)
(plot_abundanceContig_gene_nt | plot_abundanceContig_gene | plot_abundanceContig_nt | plot_abundanceContig_aa ) +
  plot_annotation(title = sample.name.INT_GRP, subtitle = paste0("(",dim(sobj)[2]," cells)"), caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
dev.off()

## Contigs Length analysis
cat("\nContigs Length analysis...\n")
for(x in list_type_contig){
  ### This should give multimodal plot
  assign(paste0("plot_lengthContig_",x,"_comb"), scRepertoire::lengthContig(tcr.combined, cloneCall=x, chains = "combined") + Seurat::NoLegend())
  ### Plots the A and B chains distribution separately
  assign(paste0("plot_lengthContig_",x,"_sin"), scRepertoire::lengthContig(tcr.combined, cloneCall=x, chains = "single") + Seurat::NoLegend())
}
### Save
png(paste0(global_output,'/lengthContig.png'), width = 800, height = 800)
(plot_lengthContig_nt_comb | plot_lengthContig_nt_sin) / (plot_lengthContig_aa_comb | plot_lengthContig_aa_sin ) +
  plot_annotation(title = sample.name.INT_GRP, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
dev.off()

## Clonal Homeostasis analysis
cat("\nClonal Homeostasis analysis...\n")
Homeo.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name.INT_GRP)

## Clonal Proportions analysis
cat("\nClonal Proportions analysis...\n")
Prop.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name.INT_GRP)

## Diversity analysis
cat("\nDiversity analysis...\n")
Div.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name.INT_GRP)

## Combine TCR data with seurat object
cat("\nCombine TCR data with seurat object...\n")
### Corresponding barcode
for(i in 1:length(tcr.combined)) tcr.combined[[i]]$barcode <- gsub(pattern = "TCR", replacement = "GE", tcr.combined[[i]]$barcode)
### Combination
sobj <- scRepertoire::combineExpression(df = tcr.combined, sc = sobj, cloneCall="aa")
sobj@meta.data$Frequency_indiv=sobj@meta.data$Frequency
### Combination tout echantillon confondu
tcr.combined_unlist <- do.call("rbind", tcr.combined)
sobj <- scRepertoire::combineExpression(df = tcr.combined_unlist, sc = sobj, cloneCall="aa")
sobj@meta.data$Frequency_all=sobj@meta.data$Frequency
sobj@meta.data$Frequency=NULL
rm(tcr.combined,tcr.combined_unlist)
### Spliting CTstrict (into separate columns for TRA-V/J/C, TRB-V/J/C and corresponding sequences, with 2 possible clonotypes) and save as metadata
### and Adding length of TR sequence to meta.data
sobj <- split.CTstrict.tcr(sobj)
### for all samples
### Plots
for(x in c("TRAV_1","TRAJ_1","TRAC_1","TRBV_1","TRBJ_1","TRBC_1","TRAV_2","TRAJ_2","TRAC_2","TRBV_2","TRBJ_2","TRBC_2")) assign(paste0("dimplot_",x), if (all(is.na(sobj@meta.data[[x]]))) patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x) else Seurat::DimPlot(sobj, group.by = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x))
for(x in c("TRA_nt_1_len","TRA_nt_2_len","TRB_nt_1_len","TRB_nt_2_len")) {
  assign(paste0("featureplot_",x),
         tryCatch({ Seurat::FeaturePlot(sobj, features = x, reduction = INT_GRP.reduction) },
            error=function(e) { return(patchwork::plot_spacer()) }
         ) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x)))
}
### Save plots
png(paste0(global_output,'/cloneType_',sample.name.INT_GRP,'.png'), width =1400, height = 3000)
( wrap_elements( (dimplot_TRAV_1 / dimplot_TRAJ_1 / dimplot_TRAC_1 / dimplot_TRBV_1 / dimplot_TRBJ_1 / dimplot_TRBC_1 / featureplot_TRA_nt_1_len / featureplot_TRB_nt_1_len) +
                  plot_annotation(title = 'TR clonotype 1', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) |
  wrap_elements( (dimplot_TRAV_2 / dimplot_TRAJ_2 /dimplot_TRAC_2 / dimplot_TRBV_2 / dimplot_TRBJ_2 / dimplot_TRBC_2 / featureplot_TRA_nt_2_len / featureplot_TRB_nt_2_len) +
                    plot_annotation(title = 'TR clonotype 2', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) +
  plot_annotation(title = sample.name.INT_GRP, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) )
dev.off()
### by samples
for (i in seq(samples.name)){
  #### selection des data du sample
  cells_sample=rownames(sobj@meta.data[sobj@meta.data$orig.ident == samples.name.GE[i],])
  sub_sobj=subset(sobj, cells = cells_sample)
  #### Plots
  for(x in c("TRAV_1","TRAJ_1","TRAC_1","TRBV_1","TRBJ_1","TRBC_1","TRAV_2","TRAJ_2","TRAC_2","TRBV_2","TRBJ_2","TRBC_2")) assign(paste0("dimplot_",x), tryCatch( {  print(Seurat::DimPlot(sub_sobj, group.by = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x)) },  error=function(err) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x)  } ))
  for(x in c("TRA_nt_1_len","TRA_nt_2_len","TRB_nt_1_len","TRB_nt_2_len")) assign(paste0("featureplot_",x), tryCatch( {  print(Seurat::FeaturePlot(sub_sobj, features = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))) },  error=function(err) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))  } ))
  #### Save plots
  png(paste0(global_output,'/cloneType_',samples.name[i],'.png'), width =1400, height = 3000)
  print( wrap_elements( (dimplot_TRAV_1 / dimplot_TRAJ_1 / dimplot_TRAC_1 / dimplot_TRBV_1 / dimplot_TRBJ_1 / dimplot_TRBC_1 / featureplot_TRA_nt_1_len / featureplot_TRB_nt_1_len) +
                          plot_annotation(title = 'TR clonotype 1', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) |
           wrap_elements( (dimplot_TRAV_2 / dimplot_TRAJ_2 /dimplot_TRAC_2 / dimplot_TRBV_2 / dimplot_TRBJ_2 / dimplot_TRBC_2 / featureplot_TRA_nt_2_len / featureplot_TRB_nt_2_len) +
                            plot_annotation(title = 'TR clonotype 2', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) +
           plot_annotation(title = samples.name[i], theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) )
  dev.off()
}

## Frequency analysis
cat("\nFrequency analysis...\n")
sobj <- Freq.g(sobj=sobj, out.dir = global_output, sample.name = sample.name.INT_GRP, reduction = INT_GRP.reduction, freq_col = "Frequency_all")
for (i in seq(samples.name)){
  #### selection des data du sample
  cells_sample=rownames(sobj@meta.data[sobj@meta.data$orig.ident == samples.name.GE[i],])
  sub_sobj=subset(sobj, cells = cells_sample)
  #### Analysis
  #Top 10 frequencies, and top 10 to top 20 frequencies
  top20_freq = sub_sobj@meta.data %>% select(Frequency_indiv,CTaa,highlight_aa_all) %>% distinct() %>% arrange(desc(Frequency_indiv)) %>% na.omit() %>% top_n(n = 20, wt = Frequency_indiv)
  if (dim(top20_freq)[1]>20) top20_freq = top20_freq[1:20,]
  rownames(top20_freq)=top20_freq$highlight_aa_all
  sub_sobj$highlight_aa_top10_freq <- ifelse(sub_sobj$highlight_aa_all %in% top20_freq$highlight[1:10], sub_sobj$highlight_aa_all, NA)
  sub_sobj$highlight_aa_top11to20_freq <- ifelse(sub_sobj$highlight_aa_all %in% top20_freq$highlight[11:length(top20_freq$highlight)], sub_sobj$highlight_aa_all, NA)
  #UMAP of top 10 frequencies
  png(paste0(global_output,'/Frequency_top_10_umap',samples.name[i],'.png'), width = 800, height = (400+350))
  print(patchwork::wrap_elements( (Seurat::DimPlot(sub_sobj, reduction = INT_GRP.reduction, group.by = "highlight_aa_top10_freq")  + Seurat::DarkTheme()) / gridExtra::tableGrob(top20_freq[1:10,c("Frequency_indiv","CTaa")], theme = gridExtra::ttheme_default(base_size = 10)) +
                                    plot_annotation(title = paste0(samples.name[i],": Top 10 Clonotypes (by frequencies)"),  subtitle = paste0("(",dim(sobj)[2]," cells)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0.5, face="bold"))) +
                                    plot_layout(heights = c(2, 1))))
  dev.off()
  #UMAP of top 11 to 20 frequencies
  png(paste0(global_output,'/Frequency_top11to20_umap',samples.name[i],'.png'), width = 800, height = (400+350))
  print(patchwork::wrap_elements( (Seurat::DimPlot(sub_sobj, reduction = INT_GRP.reduction, group.by = "highlight_aa_top11to20_freq")  + Seurat::DarkTheme()) / gridExtra::tableGrob(top20_freq[11:length(top20_freq$CTaa),c("Frequency_indiv","CTaa")], theme = gridExtra::ttheme_default(base_size = 10)) +
                                    plot_annotation(title = paste0(samples.name[i], ": Top 11 to 20 Clonotypes (by frequencies)"),  subtitle = paste0("(",dim(sobj)[2]," cells)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0.5, face="bold"))) +
                                    plot_layout(heights = c(2, 1))))
  dev.off()
}

## Physicochemical properties of the CDR3
cat("\nPhysicochemical properties of the CDR3 analysis...\n")
Physicochemical_properties.g(sobj=sobj, list_type_clT = list_type_clT, out.dir = global_output, sample.name=sample.name.INT_GRP, type='TCR')




## CLUSTERS LEVEL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cat("\nClusters Level Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
clusters_output <- paste0(output_path_TCR, "Clusters_analysis")
dir.create(path = clusters_output, recursive = TRUE, showWarnings = FALSE)

## By sample
for (i in seq(samples.name)){
  #create directory
  sample_output=paste0(clusters_output, "/", samples.name[i])
  dir.create(path = sample_output, recursive = TRUE, showWarnings = TRUE)
  
  #### selection des data du sample
  cells_sample=rownames(sobj@meta.data[sobj@meta.data$orig.ident == samples.name.GE[i],])
  sub_sobj=subset(sobj, cells = cells_sample)
  
  ## Filter cells that have no value for the x concerned + conversion to a list by clusters
  ## Need to filter, otherwise the functions count the 'NA' as a sequence.
  for(x in list_type_clT){
    if(x=="gene+nt") y="strict" else y=x
    filtred_sobj = sub_sobj[,!is.na(sub_sobj@meta.data[paste0("CT", y)])]
    assign(paste0("filtred_metadata_", sub("\\+","_",x)), scRepertoire::expression2List(sc=filtred_sobj, group=ident.name))
  }
  rm(filtred_sobj)
  
  ## Quantification of unique contig analysis
  cat("\nQuantification analysis...\n")
  sub_sobj <- Quantif.unique.c(sobj = sub_sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], filtred_metadata_aa = filtred_metadata_aa, filtred_metadata_nt = filtred_metadata_nt, filtred_metadata_gene = filtred_metadata_gene, filtred_metadata_gene_nt = filtred_metadata_gene_nt)

  ## Abundance analysis
  cat("\nAbundance analysis...\n")
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_cluster_abundanceContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::abundanceContig(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  ### Save
  png(paste0(clusters_output,'/clust_abundanceContig.png'), width = 2000, height = 600)
  (plot_cluster_abundanceContig_gene_nt | plot_cluster_abundanceContig_gene | plot_cluster_abundanceContig_nt | plot_cluster_abundanceContig_aa ) +
    plot_annotation(title = samples.name[i], caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
  dev.off()
  
  ## Clonal Homeostasis analysis
  cat("\nClonal Homeostasis analysis...\n")
  sub_sobj <- Homeo.c(sobj = sub_sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], filtred_metadata_aa = filtred_metadata_aa, filtred_metadata_nt = filtred_metadata_nt, filtred_metadata_gene = filtred_metadata_gene, filtred_metadata_gene_nt = filtred_metadata_gene_nt)
  
  ## Clonal Proportions analysis
  cat("\nClonal Proportions analysis...\n")
  sub_sobj <- Prop.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], filtred_metadata_aa = filtred_metadata_aa, filtred_metadata_nt = filtred_metadata_nt, filtred_metadata_gene = filtred_metadata_gene, filtred_metadata_gene_nt = filtred_metadata_gene_nt)
  
  ## Diversity analysis
  cat("\nDiversity analysis...\n")
  sub_sobj <- Div.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], filtred_metadata_aa = filtred_metadata_aa, filtred_metadata_nt = filtred_metadata_nt, filtred_metadata_gene = filtred_metadata_gene, filtred_metadata_gene_nt = filtred_metadata_gene_nt)
  
  ## Frequency analysis
  cat("\nFrequency analysis...\n")
  Freq.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], ident.name = ident.name, reduction=INT_GRP.reduction, freq_col = "Frequency_indiv", filtred_metadata_aa = filtred_metadata_aa, filtred_metadata_nt = filtred_metadata_nt, filtred_metadata_gene = filtred_metadata_gene, filtred_metadata_gene_nt = filtred_metadata_gene_nt)
  
  ## Clonal Overlap analysis (si plus de 1)
  cat("\nClonal Overlap analysis...\n")
  if(length(levels(Seurat::Idents(sub_sobj)))!=1 && length(unique(sub_sobj@meta.data[!is.na(sub_sobj@meta.data$CTstrict),ident.name]))!=1) Overlap.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], filtred_metadata_aa = filtred_metadata_aa, filtred_metadata_nt = filtred_metadata_nt, filtred_metadata_gene = filtred_metadata_gene, filtred_metadata_gene_nt = filtred_metadata_gene_nt)
  
  ## Physico-chemical properties of the CDR3
  cat("\nPhysico-chemical properties of the CDR3 analysis...\n")
  Physicochemical_properties.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], ident.name = ident.name, filtred_metadata_aa = filtred_metadata_aa, filtred_metadata_nt = filtred_metadata_nt, filtred_metadata_gene = filtred_metadata_gene, filtred_metadata_gene_nt = filtred_metadata_gene_nt, type='TCR')
}

#All samples
## Filter cells that have no value for the x concerned + conversion to a list by clusters
## Need to filter, otherwise the functions count the 'NA' as a sequence.
for(x in list_type_clT){
  if(x=="gene+nt") y="strict" else y=x
  filtred_sobj = sobj[,!is.na(sobj@meta.data[paste0("CT", y)])]
  assign(paste0("filtred_metadata_", sub("\\+","_",x)), scRepertoire::expression2List(sc=filtred_sobj, group=ident.name))
}
rm(filtred_sobj)

## Quantification of unique contig analysis
cat("\nQuantification analysis...\n")
sobj <- Quantif.unique.c(sobj = sobj, ident.name=ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name.INT_GRP, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Abundance analysis
cat("\nAbundance analysis...\n")
### Plots
for(x in list_type_clT) assign(paste0("plot_cluster_abundanceContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::abundanceContig(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
### Save
png(paste0(clusters_output,'/clust_abundanceContig.png'), width = 2000, height = 600)
(plot_cluster_abundanceContig_gene_nt | plot_cluster_abundanceContig_gene | plot_cluster_abundanceContig_nt | plot_cluster_abundanceContig_aa ) +
  plot_annotation(title = sample.name.INT_GRP, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
dev.off()

## Clonal Homeostasis analysis
cat("\nClonal Homeostasis analysis...\n")
sobj <- Homeo.c(sobj = sobj, ident.name=ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name.INT_GRP, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Clonal Proportions analysis
cat("\nClonal Proportions analysis...\n")
sobj <- Prop.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name.INT_GRP, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Diversity analysis
cat("\nDiversity analysis...\n")
sobj <- Div.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name.INT_GRP, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Frequency analysis
cat("\nFrequency analysis...\n")
Freq.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name.INT_GRP, ident.name=ident.name, reduction=, freq_col="Frequency", filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Clonal Overlap analysis (si plus de 1)
cat("\nClonal Overlap analysis...\n")
if(length(levels(Seurat::Idents(sobj)))!=1 && length(unique(sobj@meta.data[!is.na(sobj@meta.data$CTstrict),ident.name]))!=1) Overlap.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name.INT_GRP, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

## Physico-chemical properties of the CDR3
cat("\nPhysico-chemical properties of the CDR3 analysis...\n")
Physicochemical_properties.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name.INT_GRP, ident.name=ident.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt, type='TCR')


#renamme TCR columns with 'TCR_' prefix
toMatch <- c("^CTgene","^CTnt","^CTaa","^CTstrict","^Frequency","^cloneType","^TRAV_1","^TRAJ_1","^TRAC_1","^TRAV_2","^TRAJ_2","^TRAC_2","^TRA_nt_1","^TRA_nt_2","^TRBV_1","^TRBJ_1","^None_1","^TRBC_1","^TRBV_2","^TRBJ_2","^None_2","^TRBC_2","^TRB_nt_1","^TRB_nt_2","^TRA_nt_1_len","^TRA_nt_2_len","^TRB_nt_1_len","^TRB_nt_2_len","^highlight_aa")
matches <- grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data))
colnames(sobj@meta.data)[matches] <- paste0("TCR_", grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data), value=TRUE))

## Save packages versions
sobj@misc$technical_info$scRepertoire <- utils::packageVersion('scRepertoire')
sobj@misc$technical_info$alakazam <- utils::packageVersion('alakazam')

### Materials and Methods
if(all(file.exists(paste0(dirname(vdj.input.files.tcr), "/Materials_and_Methods.txt")))){
  tmp2=c()
  for (nb_file in 1:length(dirname(vdj.input.files.tcr))){
    tmp <- readr::read_tsv(paste0(dirname(vdj.input.files.tcr[nb_file]), "/Materials_and_Methods.txt"), col_names = FALSE)$X1
    tmp2[nb_file]=""
    for (i in 1:length(tmp)) tmp2[nb_file]=paste(tmp2[nb_file],tmp[i], sep="")
  }
  if(length(unique(tmp2)) == 1) sobj@misc$parameters$Materials_and_Methods$TCR <- tmp2[1]
} else sobj@misc$parameters$Materials_and_Methods$TCR <- NULL
sobj@misc$parameters$Materials_and_Methods$TCR <- paste0(sobj@misc$parameters$Materials_and_Methods$TCR, " The annotation was merged with corresponding cell barcode of gene expression. The scRepertoire package (version ",sobj@misc$technical_info$scRepertoire,") was used to process annotation to assign clonotype based on TCR chains. scRepertoire allows to study contig quantification, contig abundance, contig length, clonal space homeostasis, clonal proportion, clonal overlap beetween clusters and diversity. Physicochemical properties of the CDR3, based on amino-acid sequences, was determined by the alakazam R package (version ",sobj@misc$technical_info$alakazam,").")
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj=sobj, output.dir=output.dir)

### Saving GE_ADT_TCR object
cat("\nSaving object...\n")
GE_TCR_file <- paste0(output.dir, sub("\\.rda$|\\.RData$", "", basename(input.rda)), '_TCR')
save(sobj, file = paste0(GE_TCR_file, '.rda'), compress = "bzip2")
