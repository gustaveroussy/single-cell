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
  make_option("--pipeline.path", help="Path to pipeline folder")
)
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)


### Sourcing functions ####
if(is.null(args$options$pipeline.path)) stop("--pipeline.path parameter must be set!")
source(paste0(args$options$pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))


#### Get Paramaters ####
### Formatting and extract from args
args <- convert_NULL_BOOL(args)
### Project
vdj.input.files.tcr <- splitByComma.ifnotNULL(vdj.input.files.tcr)
list.author.name <- splitByComma.ifnotNULL(author.name)
list.author.mail <- splitByComma.ifnotNULL(author.mail)
#### Fixed parameters
output_path_TCR <- paste0(output.dir, "/TCR_results/")
list_type_clT <- c("strict", "gene", "nt", "aa")
list_type_contig <- c("nt","aa")
caption <- '"gene" - use the genes comprising the TCR
"nt" - use the nucleotide sequence of the CDR3 region
"aa" - use the amino acid sequence of the CDR3 region
"strict" - use the genes comprising the TCR + the nucleotide sequence of the CDR3 region for T cells. This is the proper definition of clone.'
### Clean
rm(args)

#### Check non-optional parameters ####
if (is.null(input.rda)) stop("input.rda parameter can't be empty!")
if (is.null(output.dir)) stop("output.dir parameter can't be empty!")
if (is.null(vdj.input.files.tcr)) stop("vdj.input.files.tcr parameter can't be empty!")

### Load data
load(input.rda)

#### Get Missing Paramaters ####
dimred.method <- sobj@misc$params$reductions$method
assay <- sobj@misc$params$reductions$assay
ident.name <- sobj@misc$params$clustering$ident
INT_GRP.reduction <- sobj@misc$params$clustering$umap
sample.name.INT_GRP <- Seurat::Project(sobj)
samples.name.ge <- sobj@misc$params$names.ge
samples.name <- sub("_GE","",samples.name.ge)

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

### Add authors
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)

## Remove scale.data to allow subset of seurat object
if(assay == "RNA") sobj <- reset_data_matrix(sobj, assay = assay, data = "scale.data", to_matrix = FALSE)

## GLOBAL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cat("\nGlobal Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
global_output <- paste0(output_path_TCR, "Global_analysis")
dir.create(path = global_output, recursive = TRUE, showWarnings = FALSE)

## Loading input data and Combining contigs
cat("\nLoading input data and Combining contigs...\n")
cr_res <- lapply(seq_along(vdj.input.files.tcr), load.sc.tcr.bcr, sobj = sobj, vdj.input.file = vdj.input.files.tcr, sample.name = samples.name.ge)
tcr.combined <- scRepertoire::combineTCR(input.data = cr_res, samples = samples.name, ID = rep("GE", length(samples.name)), filterMulti = TRUE, filterNonproductive = TRUE)

## Quantification of unique contig analysis
cat("\nQuantification analysis...\n")
Quantif.unique.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name.INT_GRP)

## Abundance analysis
cat("\nAbundance analysis...\n")
### Plots
for(x in list_type_clT) assign(paste0("plot_clonalAbundance_",x),patchwork::wrap_elements(scRepertoire::clonalAbundance(tcr.combined, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
for(x in list_type_clT) assign(paste0("plot_clonalAbundance_",x, "_T"),patchwork::wrap_elements(scRepertoire::clonalAbundance(tcr.combined, cloneCall = x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
### Save
png(paste0(global_output,'/clonalAbundance.png'), width = 1000, height = 600)
print(((plot_clonalAbundance_strict | plot_clonalAbundance_gene | plot_clonalAbundance_nt | plot_clonalAbundance_aa ) / (plot_clonalAbundance_strict_T | plot_clonalAbundance_gene_T | plot_clonalAbundance_nt_T | plot_clonalAbundance_aa_T ))+
plot_annotation(title = sample.name.INT_GRP, subtitle = paste0("(",dim(sobj)[2]," cells)"), caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
dev.off()

## Contigs Length analysis
cat("\nContigs Length analysis...\n")
for(x in list_type_contig){
    ### This should give multimodal plot
    assign(paste0("plot_clonalLength_",x,"_both"), scRepertoire::clonalLength(tcr.combined, cloneCall=x, chain = "both") + Seurat::NoLegend() + ggplot2::ggtitle('Both'))
    ### Plots the TRA chain distribution
    assign(paste0("plot_clonalLength_",x,"_TRA"), scRepertoire::clonalLength(tcr.combined, cloneCall=x, chain = "TRA") + Seurat::NoLegend() + ggplot2::ggtitle('TRA'))
    ### Plots the TRB chain distribution
    assign(paste0("plot_clonalLength_",x,"_TRB"), scRepertoire::clonalLength(tcr.combined, cloneCall=x, chain = "TRB") + Seurat::NoLegend() + ggplot2::ggtitle('TRB'))
}
### Save
png(paste0(global_output,'/clonalLength.png'), width = 1200, height = 800)
print((plot_clonalLength_nt_both | plot_clonalLength_nt_TRA | plot_clonalLength_nt_TRB) / (plot_clonalLength_aa_both | plot_clonalLength_aa_TRA | plot_clonalLength_aa_TRB ) +
plot_annotation(title = sample.name.INT_GRP, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
dev.off()

## Clonal Homeostasis analysis
cat("\nClonal Homeostasis analysis...\n")
Homeo.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name.INT_GRP)

## Clonal Proportions analysis
cat("\nClonal Proportions analysis...\n")
Prop.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name.INT_GRP)

# ## Diversity analysis
# cat("\nDiversity analysis...\n")
# dir.create(paste0(global_output,"/Diversity_scores/"), recursive = TRUE, showWarnings = FALSE)
# Div.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = paste0(global_output,"/Diversity_scores/"), caption = caption, sample.name = sample.name.INT_GRP)

## Combine TCR data with seurat object
cat("\nCombine TCR data with seurat object...\n")
### Corresponding barcode
for(i in 1:length(tcr.combined)) tcr.combined[[i]]$barcode <- gsub(pattern = "TCR", replacement = "GE", tcr.combined[[i]]$barcode)
### Combination
sobj <- scRepertoire::combineExpression(input.data = tcr.combined, sc = sobj, cloneCall="aa")
sobj@meta.data$Frequency_indiv <- sobj@meta.data$clonalFrequency
### Combination tout echantillon confondu
tcr.combined_unlist <- do.call("rbind", tcr.combined)
sobj <- scRepertoire::combineExpression(input.data = tcr.combined_unlist, sc = sobj, cloneCall="aa")
sobj@meta.data$Frequency_all <- sobj@meta.data$clonalFrequency
sobj@meta.data$clonalFrequency <- NULL
rm(tcr.combined,tcr.combined_unlist)
### Spliting CTstrict (into separate columns for TRA-V/J/C, TRB-V/J/C and corresponding sequences, with 2 possible clones) and save as metadata
### and Adding length of TR sequence to meta.data
sobj <- split.CTstrict.tcr(sobj)
### for all samples
### Plots
for(x in c("TRAV","TRAJ","TRAD","TRAC","TRBV","TRBD","TRBJ","TRBC")) assign(paste0("dimplot_",x), if (all(is.na(sobj@meta.data[[x]]))) patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x) else Seurat::DimPlot(sobj, group.by = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x))
for(x in c("TRA_nt_len","TRB_nt_len")) {
    assign(paste0("featureplot_",x),
          tryCatch({  print(Seurat::FeaturePlot(sobj, features = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))) },
              error = function(e) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))  } ))
}
### Save plots
png(paste0(global_output,'/cloneType_',sample.name.INT_GRP,'_uMAP.png'), width = 1400, height = 2000)
print(( wrap_elements( (dimplot_TRAV / dimplot_TRAD / dimplot_TRAJ /dimplot_TRAC / featureplot_TRA_nt_len ) +
   plot_annotation(title = 'TRA', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) |
wrap_elements( (dimplot_TRBV / dimplot_TRBD / dimplot_TRBJ / dimplot_TRBC / featureplot_TRB_nt_len) +
   plot_annotation(title = 'TRB', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) ))
dev.off()
### by samples
for (i in seq(samples.name)){
  #### selection des data du sample
  cells_sample <- rownames(sobj@meta.data[sobj@meta.data$orig.ident == samples.name.ge[i],])
  sub_sobj <- subset(sobj, cells = cells_sample)
  ### Plots
  for(x in c("TRAV","TRAJ","TRAD","TRAC","TRBV","TRBD","TRBJ","TRBC")) assign(paste0("dimplot_",x), if (all(is.na(sub_sobj@meta.data[[x]]))) patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x) else Seurat::DimPlot(sub_sobj, group.by = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x))
  for(x in c("TRA_nt_len","TRB_nt_len")) {
      assign(paste0("featureplot_",x),
             tryCatch({  print(Seurat::FeaturePlot(sub_sobj, features = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))) },
                error = function(e) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))  } ))
  }
  ### Save plots
  png(paste0(global_output,'/cloneType_',samples.name[i],'_uMAP.png'), width = 1400, height = 3000)
  print(( wrap_elements( (dimplot_TRAV / dimplot_TRAD / dimplot_TRAJ /dimplot_TRAC / featureplot_TRA_nt_len ) +
       plot_annotation(title = 'TRA', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) |
    wrap_elements( (dimplot_TRBV / dimplot_TRBD / dimplot_TRBJ / dimplot_TRBC / featureplot_TRB_nt_len) +
       plot_annotation(title = 'TRB', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) +
    plot_annotation(title = samples.name[i], theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ))
  dev.off()
}

## Frequency analysis
cat("\nFrequency analysis...\n")
sobj <- Freq.g(sobj = sobj, out.dir = global_output, sample.name = sample.name.INT_GRP, reduction = INT_GRP.reduction, freq_col = "Frequency_all")

for (i in seq(samples.name)){
  #### selection des data du sample
  cells_sample <- rownames(sobj@meta.data[sobj@meta.data$orig.ident == samples.name.ge[i],])
  sub_sobj <- subset(sobj, cells = cells_sample)
  #### Analysis
  #Top 10 frequencies, and top 10 to top 20 frequencies
  top20_freq <- sub_sobj@meta.data %>% select(Frequency_indiv,CTaa,highlight_aa_all) %>% distinct() %>% arrange(desc(Frequency_indiv)) %>% na.omit() %>% top_n(n = 20, wt = Frequency_indiv)
  if (dim(top20_freq)[1]>20) top20_freq <- top20_freq[1:20,]
  rownames(top20_freq) <- top20_freq$highlight_aa_all
  sub_sobj$highlight_aa_top10_freq <- ifelse(sub_sobj$highlight_aa_all %in% top20_freq$highlight[1:10], sub_sobj$highlight_aa_all, NA)
  sub_sobj$highlight_aa_top11to20_freq <- ifelse(sub_sobj$highlight_aa_all %in% top20_freq$highlight[11:length(top20_freq$highlight)], sub_sobj$highlight_aa_all, NA)
  #UMAP of top 10 frequencies
  png(paste0(global_output,'/Frequency_top_10_umap',samples.name[i],'.png'), width = 800, height = (400+350))
  print(patchwork::wrap_elements( (Seurat::DimPlot(sub_sobj, reduction = INT_GRP.reduction, group.by = "highlight_aa_top10_freq")  + Seurat::DarkTheme()) / gridExtra::tableGrob(top20_freq[1:10,c("Frequency_indiv","CTaa")], theme = gridExtra::ttheme_default(base_size = 10)) +
                                    plot_annotation(title = paste0(samples.name[i],": Top 10 Clones (by frequencies)"),  subtitle = paste0("(",dim(sobj)[2]," cells)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0.5, face="bold"))) +
                                    plot_layout(heights = c(2, 1))))
  dev.off()
  #UMAP of top 11 to 20 frequencies
  png(paste0(global_output,'/Frequency_top11to20_umap',samples.name[i],'.png'), width = 800, height = (400+350))
  print(patchwork::wrap_elements( (Seurat::DimPlot(sub_sobj, reduction = INT_GRP.reduction, group.by = "highlight_aa_top11to20_freq")  + Seurat::DarkTheme()) / gridExtra::tableGrob(top20_freq[11:length(top20_freq$CTaa),c("Frequency_indiv","CTaa")], theme = gridExtra::ttheme_default(base_size = 10)) +
                                    plot_annotation(title = paste0(samples.name[i], ": Top 11 to 20 Clones (by frequencies)"),  subtitle = paste0("(",dim(sobj)[2]," cells)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0.5, face="bold"))) +
                                    plot_layout(heights = c(2, 1))))
  dev.off()
}

## Physicochemical properties of the CDR3
cat("\nPhysicochemical properties of the CDR3 analysis...\n")
Physicochemical_properties.g(sobj = sobj, list_type_clT = list_type_clT, out.dir = global_output, sample.name = sample.name.INT_GRP, type = 'TCR')




## CLUSTERS LEVEL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cat("\nClusters Level Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
clusters_output <- paste0(output_path_TCR, "Clusters_analysis")
dir.create(path = clusters_output, recursive = TRUE, showWarnings = FALSE)

## By sample
for (i in seq(samples.name)){
  #create directory
  sample_output <- paste0(clusters_output, "/", samples.name[i])
  dir.create(path = sample_output, recursive = TRUE, showWarnings = FALSE)
  
  #### selection des data du sample
  cells_sample <- rownames(sobj@meta.data[sobj@meta.data$orig.ident == samples.name.ge[i],])
  sub_sobj <- subset(sobj, cells = cells_sample)
  
  ## Quantification of unique contig analysis
  cat("\nQuantification analysis...\n")
  sub_sobj <- Quantif.unique.c(sobj = sub_sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i])

  ## Abundance analysis
  cat("\nAbundance analysis...\n")
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x),patchwork::wrap_elements(scRepertoire::clonalAbundance(sub_sobj, group.by = ident.name, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x, "_T"),patchwork::wrap_elements(scRepertoire::clonalAbundance(sub_sobj, group.by = ident.name, cloneCall = x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  ### Save
  png(paste0(clusters_output,'/clust_clonalAbundance_',samples.name[i],'.png'), width = 2000, height = 600)
  print(((plot_cluster_clonalAbundance_strict | plot_cluster_clonalAbundance_gene | plot_cluster_clonalAbundance_nt | plot_cluster_clonalAbundance_aa )/ (plot_cluster_clonalAbundance_strict_T | plot_cluster_clonalAbundance_gene_T | plot_cluster_clonalAbundance_nt_T | plot_cluster_clonalAbundance_aa_T )) +
    plot_annotation(title = samples.name[i], caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
  dev.off()

  ## Clonal Homeostasis analysis
  cat("\nClonal Homeostasis analysis...\n")
  sub_sobj <- Homeo.c(sobj = sub_sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i])

  ## Clonal Proportions analysis
  cat("\nClonal Proportions analysis...\n")
  sub_sobj <- Prop.c(sobj = sub_sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i])

  # ## Diversity analysis
  # cat("\nDiversity analysis...\n")
  # dir.create(paste0(sample_output,"/Diversity_scores/"), recursive = TRUE, showWarnings = FALSE)
  # sub_sobj <- Div.c(sobj = sub_sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = paste0(sample_output,"/Diversity_scores/"), caption = caption, sample.name = samples.name[i])

  ## Frequency analysis
  cat("\nFrequency analysis...\n")
  Freq.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], ident.name = ident.name, reduction=INT_GRP.reduction, freq_col = "Frequency_indiv")

  ## Clonal Overlap analysis (si plus de 1)
  cat("\nClonal Overlap analysis...\n")
  if(length(levels(Seurat::Idents(sub_sobj)))!=1 && length(unique(sub_sobj@meta.data[!is.na(sub_sobj@meta.data$CTstrict),ident.name]))!=1) Overlap.c(sobj = sub_sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i])

  ## Physico-chemical properties of the CDR3
  cat("\nPhysico-chemical properties of the CDR3 analysis...\n")
  Physicochemical_properties.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], ident.name = ident.name, type='TCR')
}

#All samples

## Quantification of unique contig analysis
cat("\nQuantification analysis...\n")
sobj <- Quantif.unique.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP)

## Abundance analysis
cat("\nAbundance analysis...\n")
### Plots
for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x),patchwork::wrap_elements(scRepertoire::clonalAbundance(sobj, group.by = ident.name, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x, "_T"),patchwork::wrap_elements(scRepertoire::clonalAbundance(sobj, group.by = ident.name, cloneCall = x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
### Save
png(paste0(clusters_output,'/clust_clonalAbundance.png'), width = 2000, height = 800)
print(((plot_cluster_clonalAbundance_strict | plot_cluster_clonalAbundance_gene | plot_cluster_clonalAbundance_nt | plot_cluster_clonalAbundance_aa )/ (plot_cluster_clonalAbundance_strict_T | plot_cluster_clonalAbundance_gene_T | plot_cluster_clonalAbundance_nt_T | plot_cluster_clonalAbundance_aa_T )) +
plot_annotation(title = sample.name.INT_GRP, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
dev.off()

## Clonal Homeostasis analysis
cat("\nClonal Homeostasis analysis...\n")
sobj <- Homeo.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP)

## Clonal Proportions analysis
cat("\nClonal Proportions analysis...\n")
sobj <- Prop.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP)

# ## Diversity analysis
# cat("\nDiversity analysis...\n")
# dir.create(paste0(clusters_output,"/Diversity_scores/"), recursive = TRUE, showWarnings = FALSE)
# sobj <- Div.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = paste0(clusters_output,"/Diversity_scores/"), caption = caption, sample.name = sample.name.INT_GRP)

## Frequency analysis
cat("\nFrequency analysis...\n")
Freq.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP, ident.name = ident.name, reduction=, freq_col="Frequency_all")

## Clonal Overlap analysis (si plus de 1)
cat("\nClonal Overlap analysis...\n")
if(length(levels(Seurat::Idents(sobj)))!=1 && length(unique(sobj@meta.data[!is.na(sobj@meta.data$CTstrict),ident.name]))!=1) Overlap.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP)

## Physico-chemical properties of the CDR3
cat("\nPhysico-chemical properties of the CDR3 analysis...\n")
Physicochemical_properties.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP, ident.name = ident.name, type='TCR')

## Renamme TCR columns with 'TCR_' prefix
toMatch <- c("^CTgene","^CTnt","^CTaa","^CTstrict","^Frequency","^clonalProportion","^clonalFrequency","^cloneSize","^TRAV","^TRAJ","^TRAC","^TRA_nt","^TRBV","^TRBD","^TRBJ","^TRBC","^TRB_nt","^TRA_nt_len","^TRB_nt_len","^highlight_aa")
matches <- grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data))
colnames(sobj@meta.data)[matches] <- paste0("TCR_", grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data), value=TRUE))

## Make the results tables (raw vdj data, sobj TCR data, nb cell/clust/clone and merge of this 3 tables)
cat("\nSaving tables...\n")
#get raw vdj data
cr_res <- lapply(seq_along(vdj.input.files.tcr), load.sc.tcr.bcr, sobj = sobj, vdj.input.file=vdj.input.files.tcr, sample.name=samples.name.ge)
for (i in 1:length(samples.name.ge)){
  cr_res[[i]]$barcode <- paste0(samples.name.ge[i],"_", cr_res[[i]]$barcode)
}
cr_res_unlist <- lapply(cr_res, function(x) {
  df <- x[[1]]
  df$sample_barcode <- x$barcode
  df
})
cr_res_unlist_merged <- do.call(rbind, cr_res_unlist)
cr_res_unlist_merged <- cr_res_unlist_merged[,c("barcode","is_cell","high_confidence","length","chain","v_gene","d_gene","j_gene","c_gene","full_length","productive","cdr3","cdr3_nt","reads","umis")]
# write.table(cr_res_unlist_merged, file = paste0(output_path_TCR,"/raw_vdj.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
#get sobj TCR data
TCR_col <- grep("^TCR", colnames(sobj@meta.data), value = TRUE)
TCR_col_to_remove <- grep("^TCR_highlight_aa_clust|^TCR_highlight_aa_top|^TCR_TRA_nt$|^TCR_TRB_nt$|^TCR_clonalProportion$|^TCR_clonalFrequency$|^TCR_cloneSize$|^TCR_Frequency_", colnames(sobj@meta.data), value = TRUE)
TCR_col <- TCR_col[! TCR_col %in% TCR_col_to_remove]
df_sobj <- data.frame(barcode=colnames(sobj), sobj@meta.data[,TCR_col], clusters = sobj@meta.data[[ident.name]]) #table with barcode, TCR results, and the cluster of each cell.
write.table(df_sobj, file=paste0(output_path_TCR, "/vdj_foreach_cell.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
#get table nb cell/cluster/clone
df_nb_cell <- data.frame(sobj@meta.data$TCR_highlight_aa_all, sobj@meta.data[[ident.name]])
df_nb_cell <- na.omit(df_nb_cell)
df_nb_cell <- table(df_nb_cell)
df_nb_cell <- data.frame(rbind(df_nb_cell))
colnames(df_nb_cell) <- sub("X", "nbCell_byClone_Clust", colnames(df_nb_cell))
df_nb_cell <- data.frame(TCR_highlight_aa_all = row.names(df_nb_cell), df_nb_cell)
write.table(df_nb_cell, file = paste0(output_path_TCR,"/nb_cell_byclone_bycluster.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
#merge of 3 previous table
df_merged <- merge(df_sobj, df_nb_cell, by = "TCR_highlight_aa_all", all.x = TRUE, all.y = TRUE)
df_merged$TCR_highlight_aa_all <- NULL
df_merged <- merge(cr_res_unlist_merged, df_merged, by = "barcode", all.x = TRUE, all.y = FALSE)
write.table(df_merged, file = paste0(output_path_TCR,"/vdj_merged.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

## Save packages versions
sobj@misc$technical_info$scRepertoire <- utils::packageVersion('scRepertoire')
sobj@misc$technical_info$alakazam <- utils::packageVersion('alakazam')

### Materials and Methods
if(all(file.exists(paste0(dirname(vdj.input.files.tcr), "/Materials_and_Methods.txt")))){
  tmp2=c()
  for (nb_file in 1:length(dirname(vdj.input.files.tcr))){
    tmp <- readr::read_tsv(paste0(dirname(vdj.input.files.tcr[nb_file]), "/Materials_and_Methods.txt"), col_names = FALSE, show_col_types = FALSE)$X1
    tmp2[nb_file]=""
    for (i in 1:length(tmp)) tmp2[nb_file] <- paste(tmp2[nb_file],tmp[i], sep="")
  }
  if(length(unique(tmp2)) == 1) sobj@misc$parameters$Materials_and_Methods$TCR <- tmp2[1]
} else sobj@misc$parameters$Materials_and_Methods$TCR <- NULL
sobj@misc$parameters$Materials_and_Methods$TCR <- paste0(sobj@misc$parameters$Materials_and_Methods$TCR, " The annotation was merged with corresponding cell barcode of gene expression. The scRepertoire package (version ",sobj@misc$technical_info$scRepertoire,") was used to process annotation to assign clone based on TCR chains. scRepertoire allows to study contig quantification, contig abundance, contig length, clonal space homeostasis, clonal proportion, clonal overlap beetween clusters and diversity. Physicochemical properties of the CDR3, based on amino-acid sequences, was determined by the alakazam R package (version ",sobj@misc$technical_info$alakazam,").")
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj = sobj, output.dir = output.dir)

### Saving GE_ADT_TCR object
cat("\nSaving object...\n")
GE_TCR_file <- paste0(output.dir, sub("\\.rda$|\\.RData$", "", basename(input.rda)), '_TCR')
save(sobj, file = paste0(GE_TCR_file, '.rda'), compress = "bzip2")
