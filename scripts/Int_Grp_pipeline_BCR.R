#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda", help="Input seurat object (in .rda format)."),
  make_option("--output.dir", help="Output path"),
  make_option("--vdj.input.files.bcr", help="File filtered_contig_annotations.csv from CellRanger aligment pipeline."),
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
vdj.input.files.bcr <- splitByComma.ifnotNULL(vdj.input.files.bcr)
list.author.name <- splitByComma.ifnotNULL(author.name)
list.author.mail <- splitByComma.ifnotNULL(author.mail)
#### Fixed parameters
output_path_BCR <- paste0(output.dir, "/BCR_results/")
list_type_clT <- c("strict", "gene", "nt", "aa")
list_type_contig <- c("nt","aa")
caption <- '"gene" - use the genes comprising the BCR
"nt" - use the nucleotide sequence of the CDR3 region
"aa" - use the amino acid sequence of the CDR3 region
"strict" - use the genes comprising the BCR + the nucleotide sequence of the CDR3 region for T cells. This is the proper definition of clone.'
### Clean
rm(args)

#### Check non-optional parameters ####
if (is.null(input.rda)) stop("input.rda parameter can't be empty!")
if (is.null(output.dir)) stop("output.dir parameter can't be empty!")
if (is.null(vdj.input.files.bcr)) stop("vdj.input.files.bcr parameter can't be empty!")

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
print(paste0("vdj.input.files.bcr : ",vdj.input.files.bcr))
print("###########################################")

## Load libraries
require(patchwork)
suppressMessages(require(Seurat))
library(dplyr)

## Set the seed
set.seed(sobj@misc$params$seed)

## Add authors
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)

## Remove scale.data to allow subset of seurat object
if(assay == "RNA") sobj <- reset_data_matrix(sobj, assay = assay, data = "scale.data")

## GLOBAL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cat("\nGlobal Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
global_output <- paste0(output_path_BCR, "Global_analysis")
dir.create(path = global_output, recursive = TRUE, showWarnings = FALSE)

## Loading input data and Combining contigs
cat("\nLoading input data and Combining contigs...\n")
cr_res <- lapply(seq_along(vdj.input.files.bcr), load.sc.tcr.bcr, sobj = sobj, vdj.input.file = vdj.input.files.bcr, sample.name = samples.name.ge)
bcr.combined <- scRepertoire::combineBCR(input.data = cr_res, samples = samples.name, ID = rep("GE", length(samples.name)), filterMulti = TRUE, filterNonproductive = TRUE)
bcr.combined <- lapply(bcr.combined, function(df) {
  df$CTstrict <- sub("^cluster", "clone", df$CTstrict)
  df
})

## Quantification of unique contig analysis
cat("\nQuantification analysis...\n")
Quantif.unique.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name.INT_GRP)

## Abundance analysis
cat("\nAbundance analysis...\n")
### Plots
for(x in list_type_clT) assign(paste0("plot_clonalAbundance_",x),patchwork::wrap_elements(scRepertoire::clonalAbundance(bcr.combined, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
for(x in list_type_clT) assign(paste0("plot_clonalAbundance_",x, "_T"),patchwork::wrap_elements(scRepertoire::clonalAbundance(bcr.combined, cloneCall = x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
### Save
png(paste0(global_output,'/clonalAbundance.png'), width = 1000, height = 600)
print(((plot_clonalAbundance_strict | plot_clonalAbundance_gene | plot_clonalAbundance_nt | plot_clonalAbundance_aa ) / (plot_clonalAbundance_strict_T | plot_clonalAbundance_gene_T | plot_clonalAbundance_nt_T | plot_clonalAbundance_aa_T ))+
plot_annotation(title = sample.name.INT_GRP, subtitle = paste0("(",dim(sobj)[2]," cells)"), caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
dev.off()

## Contigs Length analysis
cat("\nContigs Length analysis...\n")
for(x in list_type_contig){
    ### This should give multimodal plot
    assign(paste0("plot_clonalLength_",x,"_both"), scRepertoire::clonalLength(bcr.combined, cloneCall=x, chain = "both") + Seurat::NoLegend() + ggplot2::ggtitle('Both'))
    ### Plots the TRA chain distribution
    assign(paste0("plot_clonalLength_",x,"_IGH"), scRepertoire::clonalLength(bcr.combined, cloneCall=x, chain = "IGH") + Seurat::NoLegend() + ggplot2::ggtitle('TRA'))
    ### Plots the TRB chain distribution
    assign(paste0("plot_clonalLength_",x,"_IGL"), scRepertoire::clonalLength(bcr.combined, cloneCall=x, chain = "IGL") + Seurat::NoLegend() + ggplot2::ggtitle('TRB'))
}
### Save
png(paste0(global_output,'/clonalLength.png'), width = 1200, height = 800)
print((plot_clonalLength_nt_both | plot_clonalLength_nt_IGH | plot_clonalLength_nt_IGL) / (plot_clonalLength_aa_both | plot_clonalLength_aa_IGH | plot_clonalLength_aa_IGL ) +
plot_annotation(title = sample.name.INT_GRP, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
dev.off()

## Clonal Homeostasis analysis
cat("\nClonal Homeostasis analysis...\n")
Homeo.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name.INT_GRP)

## Clonal Proportions analysis
cat("\nClonal Proportions analysis...\n")
Prop.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name.INT_GRP)

# ## Diversity analysis
# cat("\nDiversity analysis...\n")
# dir.create(paste0(global_output,"/Diversity_scores/"), recursive = TRUE, showWarnings = FALSE)
# Div.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = paste0(global_output,"/Diversity_scores/"), caption = caption, sample.name = sample.name.INT_GRP)

## Combine BCR data with seurat object
cat("\nCombine BCR data with seurat object...\n")
### Corresponding barcode
for(i in 1:length(bcr.combined)) bcr.combined[[i]]$barcode <- gsub(pattern = "BCR", replacement = "GE", bcr.combined[[i]]$barcode)
### Combination
sobj <- scRepertoire::combineExpression(input.data = bcr.combined, sc = sobj, cloneCall="aa")
sobj@meta.data$Frequency_indiv <- sobj@meta.data$clonalFrequency
### Combination tout echantillon confondu
bcr.combined_unlist <- do.call("rbind", bcr.combined)
sobj <- scRepertoire::combineExpression(input.data = bcr.combined_unlist, sc = sobj, cloneCall="aa")
sobj@meta.data$Frequency_all <- sobj@meta.data$clonalFrequency
sobj@meta.data$clonalFrequency <- NULL
rm(bcr.combined,bcr.combined_unlist)
### Spliting CTnt and CTgenes (into separate columns for IG-V/L and corresponding sequences) and save as metadata
### and Adding length of IG sequence to meta.data
sobj <- split.bcr(sobj)
### for all samples
### Plots
for(x in c('IGHV','IGHJ','IGHD','Isotype','IGLV','IGLJ','IGLC')) assign(paste0("dimplot_",x), if (all(is.na(sobj@meta.data[[x]]))) patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x) else Seurat::DimPlot(sobj, group.by = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x))
for(x in c("cdr3_nt_IGH_len", "cdr3_nt_IGL_len"))  {
    assign(paste0("featureplot_",x),
          tryCatch({  print(Seurat::FeaturePlot(sobj, features = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))) },
              error = function(e) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))  } ))
}
### Save plots
png(paste0(global_output,'/cloneType_',sample.name.INT_GRP,'_uMAP.png'), width =700, height = 3000)
print(( wrap_elements(dimplot_IGHV / dimplot_IGHJ / dimplot_IGHD / dimplot_Isotype / dimplot_IGLV / dimplot_IGLJ / dimplot_IGLC / featureplot_cdr3_nt_IGH_len / featureplot_cdr3_nt_IGL_len) +
    plot_annotation(title = 'BCR', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ))
dev.off()
### by samples
for (i in seq(samples.name)){
  #### selection des data du sample
  cells_sample <- rownames(sobj@meta.data[sobj@meta.data$orig.ident == samples.name.ge[i],])
  sub_sobj <- subset(sobj, cells = cells_sample)
  #### Plots
  for(x in c('IGHV','IGHJ','IGHD','Isotype','IGLV','IGLJ','IGLC')) assign(paste0("dimplot_",x), tryCatch( {  print(Seurat::DimPlot(sub_sobj, group.by = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x)) },  error=function(err) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x)  } ))
  for(x in c("cdr3_nt_IGH_len", "cdr3_nt_IGL_len")) {
       assign(paste0("featureplot_",x),
             tryCatch({  print(Seurat::FeaturePlot(sub_sobj, features = x, reduction = INT_GRP.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))) },
                error = function(e) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))  } ))
  }
  #### Save plots
  png(paste0(global_output,'/cloneType_',samples.name[i],'_uMAP.png'), width =1400, height = 3000)
  print( ( wrap_elements(dimplot_IGHV / dimplot_IGHJ / dimplot_IGHD / dimplot_Isotype / dimplot_IGLV / dimplot_IGLJ / dimplot_IGLC / featureplot_cdr3_nt_IGH_len / featureplot_cdr3_nt_IGL_len) +
             plot_annotation(title = 'BCR', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ))
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
Physicochemical_properties.g(sobj = sobj, list_type_clT = list_type_clT, out.dir = global_output, sample.name = sample.name.INT_GRP, type = 'BCR')





## CLUSTERS LEVEL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
cat("\nClusters Level Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
clusters_output <- paste0(output_path_BCR, "Clusters_analysis")
dir.create(path = clusters_output, recursive = TRUE, showWarnings = FALSE)

## By sample
for (i in seq(samples.name)){
  #create directory
  sample_output <- paste0(clusters_output, "/", samples.name[i])
  dir.create(path = sample_output, recursive = TRUE, showWarnings = FALSE)
  
  #### selection des data du sample
  cells_sample <- rownames(sobj@meta.data[sobj@meta.data$orig.ident == samples.name.ge[i],])
  sub_sobj <- subset(sobj, cells = cells_sample)
  
  ## Quantification of unique contig analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors la quantif d'unique sera de 100% partout!!!!!!!!!!!!!!!!!!
  cat("\nQuantification analysis...\n")
  sub_sobj <- Quantif.unique.c(sobj = sub_sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i])
  
  ## Abundance analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors l'abondance sera de 1 partout!!!!!!!!!!!!!!!!!!!!!
  cat("\nAbundance analysis...\n")
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x),patchwork::wrap_elements(scRepertoire::clonalAbundance(sub_sobj, group.by = ident.name, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x, "_T"),patchwork::wrap_elements(scRepertoire::clonalAbundance(sub_sobj, group.by = ident.name, cloneCall = x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  ### Save
  png(paste0(clusters_output,'/clust_clonalAbundance_',samples.name[i],'.png'), width = 2000, height = 600)
  print((plot_cluster_clonalAbundance_strict | plot_cluster_clonalAbundance_gene | plot_cluster_clonalAbundance_nt | plot_cluster_clonalAbundance_aa ) +
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
  Freq.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], ident.name = ident.name, reduction = INT_GRP.reduction, freq_col = "Frequency_indiv")

  ## Clonal Overlap analysis (si plus de 1)
  cat("\nClonal Overlap analysis...\n")
  if(length(levels(Seurat::Idents(sub_sobj)))!=1 && length(unique(sub_sobj@meta.data[!is.na(sub_sobj@meta.data$CTstrict),ident.name]))!=1) Overlap.c(sobj = sub_sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i])
  
  ## Physico-chemical properties of the CDR3
  cat("\nPhysico-chemical properties of the CDR3 analysis...\n")
  Physicochemical_properties.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption = caption, sample.name = samples.name[i], ident.name = ident.name, type='BCR')
}

#All samples

## Quantification of unique contig analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors la quantif d'unique sera de 100% partout!!!!!!!!!!!!!!!!!!
cat("\nQuantification analysis...\n")
sobj <- Quantif.unique.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP)

## Abundance analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors l'abondance sera de 1 partout!!!!!!!!!!!!!!!!!!!!!
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

## Frequency analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors la fréquence vaut 1 partout!!!!!!!!!!!!!!!
cat("\nFrequency analysis...\n")
Freq.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP, ident.name = ident.name, reduction=INT_GRP.reduction, freq_col="Frequency_all")

## Clonal Overlap analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors pas d'overlap!!!!!!!!!!!!!!!
cat("\nClonal Overlap analysis...\n")
if(length(levels(Seurat::Idents(sobj)))!=1 && length(unique(sobj@meta.data[!is.na(sobj@meta.data$CTstrict),ident.name]))!=1) Overlap.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP)

## Physico-chemical properties of the CDR3
cat("\nPhysico-chemical properties of the CDR3 analysis...\n")
Physicochemical_properties.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name.INT_GRP, ident.name = ident.name, type='BCR')

## Renamme BCR columns with 'BCR_' prefix
toMatch <- c("^CTgene","^CTnt","^CTaa","^CTstrict","^Frequency","^clonalProportion","^clonalFrequency","^cloneSize","^IGHV","^IGHJ","^IGHD","^Isotype","^IGLV","^IGLJ","^IGLC","^cdr3_nt_IGH", "^cdr3_nt_IGL","^cdr3_nt_IGH_len","^cdr3_nt_IGL_len","^highlight_aa")
matches <- grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data))
colnames(sobj@meta.data)[matches]  <- paste0("BCR_", grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data), value=TRUE))

## Make the results tables (raw vdj data, sobj BCR data, nb cell/clust/clone and merge of this 3 tables)
cat("\nSaving tables...\n")
#get raw vdj data
cr_res <- lapply(seq_along(vdj.input.files.bcr), load.sc.tcr.bcr, sobj = sobj, vdj.input.file = vdj.input.files.bcr, sample.name = samples.name.ge)
for (i in 1:length(samples.name.ge)){
  cr_res[[i]]$barcode=paste0(samples.name.ge[i],"_", cr_res[[i]]$barcode)
}
cr_res_unlist <- lapply(cr_res, function(x) {
  df <- x[[1]]
  df$sample_barcode <- x$barcode
  df
})
cr_res_unlist_merged <- do.call(rbind, cr_res_unlist)
cr_res_unlist_merged <- cr_res_unlist_merged[,c("barcode","is_cell","high_confidence","length","chain","v_gene","d_gene","j_gene","c_gene","full_length","productive","cdr3","cdr3_nt","reads","umis")]
# write.table(cr_res_unlist_merged, file = paste0(output_path_BCR,"/raw_vdj.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
#get sobj BCR data
BCR_col = grep("^BCR", colnames(sobj@meta.data), value = TRUE)
BCR_col_to_remove <- grep("^BCR_highlight_aa_clust|^BCR_highlight_aa_top|^BCR_clonalProportion$|^BCR_clonalFrequency$|^BCR_cloneSize$|^BCR_Frequency_", colnames(sobj@meta.data), value = TRUE)
BCR_col = BCR_col[! BCR_col %in% BCR_col_to_remove]
df_sobj <- data.frame(barcode=colnames(sobj),sobj@meta.data[,BCR_col],clusters=sobj@meta.data[[ident.name]])
write.table(df_sobj, file=paste0(output_path_BCR,"/vdj_foreach_cell.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
#get table nb cell/cluster/clone
df_nb_cell <- data.frame(sobj@meta.data$BCR_highlight_aa_all,sobj@meta.data[[ident.name]])
df_nb_cell <- na.omit(df_nb_cell)
df_nb_cell <- table(df_nb_cell)
df_nb_cell <- data.frame(rbind(df_nb_cell))
colnames(df_nb_cell) <- sub("X", "nbCell_byClono_Clust",colnames(df_nb_cell))
df_nb_cell <- data.frame(BCR_highlight_aa_all=row.names(df_nb_cell),df_nb_cell)
write.table(df_nb_cell, file=paste0(output_path_BCR,"/nb_cell_byclone_bycluster.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
##merge of 3 previous table #error if NA into BCR_highlight_aa_all column
#df_merged <- merge(df_sobj, df_nb_cell, by = "BCR_highlight_aa_all",all.x=TRUE, all.y=TRUE)
#df_merged$BCR_highlight_aa_all <- NULL
#df_merged <- merge(cr_res_unlist_merged, df_merged, by = "barcode",all.x=TRUE, all.y=FALSE)
#write.table(df_merged, file=paste0(output_path_BCR,"/vdj_merged.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")

## Save packages versions
sobj@misc$technical_info$scRepertoire <- utils::packageVersion('scRepertoire')
sobj@misc$technical_info$alakazam <- utils::packageVersion('alakazam')

### Materials and Methods
if(all(file.exists(paste0(dirname(vdj.input.files.bcr), "/Materials_and_Methods.txt")))){
  tmp2=c()
  for (nb_file in 1:length(dirname(vdj.input.files.bcr))){
    tmp <- readr::read_tsv(paste0(dirname(vdj.input.files.bcr[nb_file]), "/Materials_and_Methods.txt"), col_names = FALSE, show_col_types = FALSE)$X1
    tmp2[nb_file]=""
    for (i in 1:length(tmp)) tmp2[nb_file] <- paste(tmp2[nb_file],tmp[i], sep="")
  }
  if(length(unique(tmp2)) == 1) sobj@misc$parameters$Materials_and_Methods$BCR <- tmp2[1]
} else sobj@misc$parameters$Materials_and_Methods$BCR <- NULL
sobj@misc$parameters$Materials_and_Methods$BCR <- paste0(sobj@misc$parameters$Materials_and_Methods$BCR, " The annotation was merged with corresponding cell barcode of gene expression. The scRepertoire package (version ",sobj@misc$technical_info$scRepertoire,") was used to process annotation to assign clone based on Ig chains. scRepertoire allows to study contig quantification, contig abundance, contig length, clonal space homeostasis, clonal proportion, clonal overlap beetween clusters. Physicochemical properties of the CDR3, based on amino-acid sequences, was determined by the alakazam R package (version ",sobj@misc$technical_info$alakazam,").")
if(!is.null(sobj@misc$parameters$Materials_and_Methods$TCR)){ #merge TCR/BCR is Materials and Methods are the same
  MandM_TCR <- gsub("TCR chains", "", sobj@misc$parameters$Materials_and_Methods$TCR)
  MandM_BCR <- gsub("Ig chains", "", sobj@misc$parameters$Materials_and_Methods$BCR)
  if(MandM_TCR == MandM_BCR){
      sobj@misc$parameters$Materials_and_Methods$Immune_profiling <- gsub("TCR chains", "TCR or Ig chains", sobj@misc$parameters$Materials_and_Methods$TCR)
      sobj@misc$parameters$Materials_and_Methods$TCR <- NULL
      sobj@misc$parameters$Materials_and_Methods$BCR <- NULL
  }
}
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj = sobj, output.dir=output.dir)

### Saving GE_ADT_BCR object
cat("\nSaving object...\n")
GE_BCR_file <- paste0(output.dir, sub("\\.rda$|\\.RData$", "", basename(input.rda)), '_BCR')
save(sobj, file = paste0(GE_BCR_file, '.rda'), compress = "bzip2")
