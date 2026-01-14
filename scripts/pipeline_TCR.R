#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda", help="Input seurat object (in .rda format)."),
  make_option("--output.dir", help="Output path"),
  make_option("--vdj.input.file.tcr", help="File filtered_contig_annotations.csv from CellRanger aligment pipeline."),
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
if (is.null(vdj.input.file.tcr)) stop("vdj.input.file.tcr parameter can't be empty!")

### Load data
load(input.rda)

#### Get Missing Paramaters ####
dimred.method <- sobj@misc$params$reductions$method
GE_file <- sub("\\.rda$|\\.RData$", "", input.rda)
dimred.method <- sobj@misc$params$reductions$method
ident.name <- sobj@misc$params$clustering$ident
RNA.reduction <- sobj@misc$params$clustering$umap
sample.name <- sub("_GE", "", sobj@misc$params$sample.name.ge)

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

### Add authors
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)

## Loading input data and Combining contigs
cat("\nLoading input data and Combining contigs...\n")
cr_res <- load.sc.tcr.bcr(sobj = sobj, vdj.input.file = vdj.input.file.tcr)
if(dim(cr_res[[1]])[1]!=0){ #if there are TCR
  tcr.combined <- scRepertoire::combineTCR(input.data = cr_res, samples = sample.name, ID = "TCR", filterMulti = FALSE, filterNonproductive = FALSE)

  ## GLOBAL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  cat("\nGlobal Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
  global_output <- paste0(output_path_TCR, "Global_analysis")
  dir.create(path = global_output, recursive = TRUE, showWarnings = FALSE)

  ## Quality analysis
  cat("\nQuantification analysis...\n")
  QC.tcr.bcr(cr_res = cr_res, out.dir = global_output, type = "TCR")
  
  ## Filtering data
  tcr.combined <- scRepertoire::combineTCR(input.data = cr_res, samples = sample.name, ID = "TCR", filterMulti = TRUE, filterNonproductive = TRUE)
  
  ## Quantification of unique contig analysis
  Quantif.unique.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name)

  ## Abundance analysis
  cat("\nAbundance analysis...\n")
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_clonalAbundance_",x),patchwork::wrap_elements(scRepertoire::clonalAbundance(tcr.combined, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
  for(x in list_type_clT) assign(paste0("plot_clonalAbundance_",x, "_T"),patchwork::wrap_elements(scRepertoire::clonalAbundance(tcr.combined, cloneCall = x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
  ### Save
  png(paste0(global_output,'/clonalAbundance.png'), width = 1000, height = 600)
  print(((plot_clonalAbundance_strict | plot_clonalAbundance_gene | plot_clonalAbundance_nt | plot_clonalAbundance_aa ) / (plot_clonalAbundance_strict_T | plot_clonalAbundance_gene_T | plot_clonalAbundance_nt_T | plot_clonalAbundance_aa_T ))+
    plot_annotation(title = sample.name, subtitle = paste0("(",dim(sobj)[2]," cells)"), caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
  dev.off()

  ## Contigs Length analysis
  cat("\nContigs Length analysis...\n")
  for(x in list_type_contig){
    ### This should give multimodal plot
    assign(paste0("plot_clonalLength_",x,"_both"), scRepertoire::clonalLength(tcr.combined, cloneCall=x, chain = "both") + Seurat::NoLegend() + ggplot2::ggtitle('Both'))
    ### Plots the TRA chain distribution
    assign(paste0("plot_clonalLength_",x,"_TRA"), scRepertoire::clonalLength(tcr.combined, cloneCall=x, chain = "TRA") + Seurat::NoLegend() + ggplot2::ggtitle('TRA'))
    ### Plots the TRG chain distribution
    assign(paste0("plot_clonalLength_",x,"_TRB"), scRepertoire::clonalLength(tcr.combined, cloneCall=x, chain = "TRB") + Seurat::NoLegend() + ggplot2::ggtitle('TRB'))
  }
  ### Save
  png(paste0(global_output,'/clonalLength.png'), width = 1200, height = 800)
  print((plot_clonalLength_nt_both | plot_clonalLength_nt_TRA | plot_clonalLength_nt_TRB) / (plot_clonalLength_aa_both | plot_clonalLength_aa_TRA | plot_clonalLength_aa_TRB ) +
    plot_annotation(title = sample.name, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
  dev.off()

  ## Clonal Homeostasis analysis
  cat("\nClonal Homeostasis analysis...\n")
  Homeo.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name)

  ## Clonal Proportions analysis
  cat("\nClonal Proportions analysis...\n")
  Prop.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name)

  # ## Diversity analysis
  # cat("\nDiversity analysis...\n")
  # dir.create(paste0(global_output,"/Diversity_scores/"), recursive = TRUE, showWarnings = FALSE)
  # Div.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = paste0(global_output,"/Diversity_scores/"), caption = caption, sample.name = sample.name)

  ## Combine TCR data with seurat object
  cat("\nCombine TCR data with seurat object...\n")
  ### Corresponding barcode
  tcr.combined[[1]]$barcode <- gsub(pattern = paste0(sample.name, "_TCR_"), replacement = '', tcr.combined[[1]]$barcode)
  ### Combination
  sobj <- scRepertoire::combineExpression(input.data = tcr.combined, sc = sobj, cloneCall="aa")

  ### Spliting CTstrict (into separate columns for TRA-V/J/C, TRB-V/J/C and corresponding sequences) and save as metadata
  ### and Adding length of TR sequence to meta.data
  sobj <- split.CTstrict.tcr(sobj)

  ### Plots
  for(x in c("TRAV","TRAJ","TRAD","TRAC","TRBV","TRBD","TRBJ","TRBC")) assign(paste0("dimplot_",x), if (all(is.na(sobj@meta.data[[x]]))) patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x) else Seurat::DimPlot(sobj, group.by = x, reduction = RNA.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x))
  for(x in c("TRA_nt_len","TRB_nt_len")) {
    assign(paste0("featureplot_",x),
            tryCatch({  print(Seurat::FeaturePlot(sobj, features = x, reduction = RNA.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))) },
                error = function(e) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))  } ))
  }
  ### Save plots
  png(paste0(global_output,'/cloneType_uMAP.png'), width = 1400, height = 2000)
  print( wrap_elements( (dimplot_TRAV / dimplot_TRAD / dimplot_TRAJ /dimplot_TRAC / featureplot_TRA_nt_len ) +
       plot_annotation(title = 'TRA', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) |
    wrap_elements( (dimplot_TRBV / dimplot_TRBD / dimplot_TRBJ / dimplot_TRBC / featureplot_TRB_nt_len) +
       plot_annotation(title = 'TRB', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) )
  dev.off()

  ## Frequency analysis
  cat("\nFrequency analysis...\n")
  sobj <- Freq.g(sobj = sobj, out.dir = global_output, sample.name = sample.name, reduction = RNA.reduction, freq_col = "clonalFrequency")

  ## Physicochemical properties of the CDR3
  cat("\nPhysicochemical properties of the CDR3 analysis...\n")
  Physicochemical_properties.g(sobj = sobj, list_type_clT = list_type_clT, out.dir = global_output, sample.name = sample.name, type = 'TCR')




  ## CLUSTERS LEVEL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  cat("\nClusters Level Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
  clusters_output <- paste0(output_path_TCR, "Clusters_analysis")
  dir.create(path = clusters_output, recursive = TRUE, showWarnings = FALSE)

  ## Quantification of unique contig analysis
  cat("\nQuantification analysis...\n")
  sobj <- Quantif.unique.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name)

  ## Abundance analysis
  cat("\nAbundance analysis...\n")
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x),patchwork::wrap_elements(scRepertoire::clonalAbundance(sobj, group.by = ident.name, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x, "_T"),patchwork::wrap_elements(scRepertoire::clonalAbundance(sobj, group.by = ident.name, cloneCall = x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  ### Save
  png(paste0(clusters_output,'/clust_clonalAbundance_',sample.name,'.png'), width = 2000, height = 800)
  print(((plot_cluster_clonalAbundance_strict | plot_cluster_clonalAbundance_gene | plot_cluster_clonalAbundance_nt | plot_cluster_clonalAbundance_aa )/ (plot_cluster_clonalAbundance_strict_T | plot_cluster_clonalAbundance_gene_T | plot_cluster_clonalAbundance_nt_T | plot_cluster_clonalAbundance_aa_T )) +
    plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
  dev.off()

  ## Clonal Homeostasis analysis
  cat("\nClonal Homeostasis analysis...\n")
  sobj <- Homeo.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name)

  ## Clonal Proportions analysis
  cat("\nClonal Proportions analysis...\n")
  sobj <- Prop.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name)

  # # Diversity analysis
  # cat("\nDiversity analysis...\n")
  # dir.create(paste0(clusters_output,"/Diversity_scores/"), recursive = TRUE, showWarnings = FALSE)
  # sobj <- Div.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = paste0(clusters_output,"/Diversity_scores/"), caption = caption, sample.name = sample.name)

  ## Frequency analysis
  cat("\nFrequency analysis...\n")
  Freq.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name, ident.name = ident.name, reduction=RNA.reduction, freq_col = "clonalFrequency")

  ## Clonal Overlap analysis (si plus de 1)
  cat("\nClonal Overlap analysis...\n")
  if(length(levels(Seurat::Idents(sobj)))!=1 && length(unique(sobj@meta.data[!is.na(sobj@meta.data$CTstrict),ident.name]))!=1) Overlap.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name)

  ## Physico-chemical properties of the CDR3
  cat("\nPhysico-chemical properties of the CDR3 analysis...\n")
  Physicochemical_properties.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name, ident.name = ident.name, type='TCR')

  #renamme TCR columns with 'TCR_' prefix
  toMatch <- c("^CTgene","^CTnt","^CTaa","^CTstrict","^clonalProportion","^clonalFrequency","^cloneSize","^TRAV","^TRAJ","^TRAC","^TRA_nt","^TRBV","^TRBD","^TRBJ","^TRBC","^TRB_nt","^TRA_nt_len","^TRB_nt_len","^highlight_aa")
  matches <- grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data))
  colnames(sobj@meta.data)[matches] <- paste0("TCR_", grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data), value=TRUE))

  ## Make the results tables (raw vdj data, sobj BCR data, nb cell/clust/clone and merge of this 3 tables)
  cat("\nSaving tables...\n")
  #get raw vdj data
  cr_res[[1]] <- cr_res[[1]][,c("barcode","is_cell","high_confidence","length","chain","v_gene","d_gene","j_gene","c_gene","full_length","productive","cdr3","cdr3_nt","reads","umis")]
  # write.table(cr_res[[1]], file = paste0(output_path_TCR,"/raw_vdj.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
  #get sobj TCR data
  TCR_col <- grep("^TCR_", colnames(sobj@meta.data), value = TRUE)
  TCR_col_to_remove <- grep("^TCR_highlight_aa_clust|^TCR_highlight_aa_top|^TCR_TRA_nt$|^TCR_TRB_nt$|^TCR_clonalProportion$|^TCR_clonalFrequency$|^TCR_cloneSize$", colnames(sobj@meta.data), value = TRUE)
  TCR_col <- TCR_col[! TCR_col %in% TCR_col_to_remove]
  df_sobj <- data.frame(barcode = colnames(sobj), sobj@meta.data[,TCR_col], cluster = sobj@meta.data[[ident.name]]) #table with barcode, TCR results, and the cluster of each cell.
  write.table(df_sobj, file = paste0(output_path_TCR, "/vdj_foreach_cell.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
  #get table nb cell/cluster/clone
  df_nb_cell <- data.frame(sobj@meta.data$TCR_highlight_aa_all, sobj@meta.data[[ident.name]])
  df_nb_cell <- na.omit(df_nb_cell)
  df_nb_cell <- table(df_nb_cell)
  df_nb_cell <- data.frame(rbind(df_nb_cell))
  colnames(df_nb_cell) <- sub("X", "nbCell_byClone_Clust", colnames(df_nb_cell))
  df_nb_cell <- data.frame(TCR_highlight_aa_all = row.names(df_nb_cell), df_nb_cell)
  write.table(df_nb_cell, file=paste0(output_path_TCR,"/nb_cell_byclone_bycluster.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
  ##merge of 3 previous table #error if NA into TCR_highlight_aa_all column
  #df_merged <- merge(df_sobj, df_nb_cell, by = "TCR_highlight_aa_all",all.x=TRUE, all.y=TRUE)
  #df_merged$TCR_highlight_aa_all <- NULL
  #df_merged <- merge(cr_res[[1]], df_merged, by = "barcode",all.x=TRUE, all.y=FALSE)
  #write.table(df_merged, file=paste0(output_path_TCR,"/vdj_merged.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")

  ## Save packages versions
  sobj@misc$technical_info$scRepertoire <- utils::packageVersion('scRepertoire')
  sobj@misc$technical_info$alakazam <- utils::packageVersion('alakazam')

  ### Materials and Methods
  if(file.exists(paste0(dirname(vdj.input.file.tcr), "/Materials_and_Methods.txt"))){
    tmp <- readr::read_tsv(paste0(dirname(vdj.input.file.tcr), "/Materials_and_Methods.txt"), col_names = FALSE, show_col_types = FALSE)$X1
    tmp2 <- ""
    for (i in 1:length(tmp)) tmp2 <- paste(tmp2,tmp[i], sep="")
    sobj@misc$parameters$Materials_and_Methods$TCR <- tmp2
  } else sobj@misc$parameters$Materials_and_Methods$TCR <- NULL
  sobj@misc$parameters$Materials_and_Methods$TCR <- paste0(sobj@misc$parameters$Materials_and_Methods$TCR, " The annotation was merged with corresponding cell barcode of gene expression. The scRepertoire package (version ",sobj@misc$technical_info$scRepertoire,") was used to process annotation to assign clone based on TCR chains. scRepertoire allows to study contig quantification, contig abundance, contig length, clonal space homeostasis, clonal proportion, clonal overlap beetween clusters. Physicochemical properties of the CDR3, based on amino-acid sequences, was determined by the alakazam R package (version ",sobj@misc$technical_info$alakazam,").")
  sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
  write_MandM(sobj=sobj, output.dir=output.dir)

}else message("No TCR found!")

### Saving GE_ADT_TCR object
cat("\nSaving object...\n")
GE_TCR_file <- paste0(output.dir, basename(GE_file), '_TCR')
save(sobj, file = paste0(GE_TCR_file, '.rda'), compress = "bzip2")
