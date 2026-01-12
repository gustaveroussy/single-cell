#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda", help="Input seurat object (in .rda format)."),
  make_option("--output.dir", help="Output path"),
  make_option("--vdj.input.file.bcr", help="File filtered_contig_annotations.csv from CellRanger aligment pipeline."),
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
output_path_BCR <- paste0(output.dir, "/BCR_results/")
list_type_clT <- c("strict", "gene", "nt", "aa")
list_type_contig <- c("nt","aa")
caption <- '"gene" - use the genes comprising the BCR
"nt" - use the nucleotide sequence of the CDR3 region
"aa" - use the amino acid sequence of the CDR3 region
"strict" - use the genes comprising the BCR + the nucleotide sequence of the CDR3 region for B cells.'
### Clean
rm(args)

#### Check non-optional parameters ####
if (is.null(input.rda)) stop("input.rda parameter can't be empty!")
if (is.null(output.dir)) stop("output.dir parameter can't be empty!")
if (is.null(vdj.input.file.bcr)) stop("vdj.input.file.bcr parameter can't be empty!")

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
print(paste0("vdj.input.file.bcr : ",vdj.input.file.bcr))
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
cr_res <- load.sc.tcr.bcr(sobj = sobj, vdj.input.file = vdj.input.file.bcr)
if(dim(cr_res[[1]])[1]!=0){ #if there are BCR
  bcr.combined <- scRepertoire::combineBCR(input.data = cr_res, samples = sample.name, ID = "BCR", filterMulti = FALSE, filterNonproductive = FALSE)
  bcr.combined[[1]]$CTstrict <- sub("^cluster", "clone", bcr.combined[[1]]$CTstrict)

  ## GLOBAL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  cat("\nGlobal Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
  global_output <- paste0(output_path_BCR, "Global_analysis")
  dir.create(path = global_output, recursive = TRUE, showWarnings = FALSE)

  ## Quality analysis
  cat("\nQuantification analysis...\n")
  QC.tcr.bcr(cr_res = cr_res, out.dir = global_output, type = "BCR")

  ## Filtering data
  bcr.combined <- scRepertoire::combineBCR(input.data = cr_res, samples = sample.name, ID = "BCR", filterMulti = TRUE, filterNonproductive = TRUE)

  ## Quantification of unique contig analysis
  Quantif.unique.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name)

  ## Abundance analysis
  cat("\nAbundance analysis...\n")
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_clonalAbundance_",x),patchwork::wrap_elements(scRepertoire::clonalAbundance(bcr.combined, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
  for(x in list_type_clT) assign(paste0("plot_clonalAbundance_",x, "_T"),patchwork::wrap_elements(scRepertoire::clonalAbundance(bcr.combined, cloneCall = x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
  ### Save
  png(paste0(global_output,'/clonalAbundance.png'), width = 1000, height = 600)
  print(((plot_clonalAbundance_strict | plot_clonalAbundance_gene | plot_clonalAbundance_nt | plot_clonalAbundance_aa ) / (plot_clonalAbundance_strict_T | plot_clonalAbundance_gene_T | plot_clonalAbundance_nt_T | plot_clonalAbundance_aa_T ))+
    plot_annotation(title = sample.name, subtitle = paste0("(",dim(sobj)[2]," cells)"), caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
  dev.off()

  ## Contigs Length analysis
  cat("\nContigs Length analysis...\n")
  for(x in list_type_contig){
    ### This should give multimodal plot
    assign(paste0("plot_clonalLength_",x,"_both"), scRepertoire::clonalLength(bcr.combined, cloneCall=x, chain = "both") + Seurat::NoLegend() + ggplot2::ggtitle('Both'))
    ### Plots the IGH chain distribution
    assign(paste0("plot_clonalLength_",x,"_IGH"), scRepertoire::clonalLength(bcr.combined, cloneCall=x, chain = "IGH") + Seurat::NoLegend() + ggplot2::ggtitle('IGH'))
    ### Plots the IGL chain distribution
    assign(paste0("plot_clonalLength_",x,"_IGL"), scRepertoire::clonalLength(bcr.combined, cloneCall=x, chain = "IGL") + Seurat::NoLegend() + ggplot2::ggtitle('IGL'))
  }
  ### Save
  png(paste0(global_output,'/clonalLength.png'), width = 1200, height = 800)
  print((plot_clonalLength_nt_both | plot_clonalLength_nt_IGH | plot_clonalLength_nt_IGL) / (plot_clonalLength_aa_both | plot_clonalLength_aa_IGH | plot_clonalLength_aa_IGL ) +
    plot_annotation(title = sample.name, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
  dev.off()

  ## Clonal Homeostasis analysis
  cat("\nClonal Homeostasis analysis...\n")
  Homeo.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name)

  ## Clonal Proportions analysis
  cat("\nClonal Proportions analysis...\n")
  Prop.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name)

  # ## Diversity analysis
  # cat("\nDiversity analysis...\n")
  # dir.create(paste0(global_output,"/Diversity_scores/"), recursive = TRUE, showWarnings = FALSE)
  # Div.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = paste0(global_output,"/Diversity_scores/"), caption = caption, sample.name = sample.name)

  ## Combine BCR data with seurat object
  cat("\nCombine BCR data with seurat object...\n")
  ### Corresponding barcode
  bcr.combined[[1]]$barcode <- gsub(pattern = paste0(sample.name, "_BCR_"), replacement = '', bcr.combined[[1]]$barcode)
  ### Combination
  sobj <- scRepertoire::combineExpression(input.data = bcr.combined, sc = sobj, cloneCall="aa")

  ### Spliting CTnt and CTgenes (into separate columns for IG-V/L and corresponding sequences) and save as metadata
  ### and Adding length of IG sequence to meta.data
  sobj <- split.bcr(sobj)

  ### Plots
  for(x in c('IGHV','IGHJ','IGHD','Isotype','IGLV','IGLJ','IGLC')) assign(paste0("dimplot_",x), if (all(is.na(sobj@meta.data[[x]]))) patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x) else Seurat::DimPlot(sobj, group.by = x, reduction = RNA.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x))
  for(x in c("cdr3_nt_IGH_len", "cdr3_nt_IGL_len"))  {
    assign(paste0("featureplot_",x),
           tryCatch({ Seurat::FeaturePlot(sobj, features = x, reduction = RNA.reduction) },
              error=function(e) { return(patchwork::plot_spacer()) }
           ) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x)))
  }
  ### Save plots
  png(paste0(global_output,'/cloneType_uMAP.png'), width =700, height = 3000)
  print(( wrap_elements(dimplot_IGHV / dimplot_IGHJ / dimplot_IGHD / dimplot_Isotype / dimplot_IGLV / dimplot_IGLJ / dimplot_IGLC / featureplot_cdr3_nt_IGH_len / featureplot_cdr3_nt_IGL_len) +
      plot_annotation(title = 'BCR', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ))
  dev.off()

  ## Frequency analysis
  cat("\nFrequency analysis...\n")
  sobj <- Freq.g(sobj = sobj, out.dir = global_output, sample.name = sample.name, reduction = RNA.reduction, freq_col = "clonalFrequency")

  ## Physicochemical properties of the CDR3
  cat("\nPhysicochemical properties of the CDR3 analysis...\n")
  Physicochemical_properties.g(sobj = sobj, list_type_clT = list_type_clT, out.dir = global_output, sample.name = sample.name, type='BCR')





  ## CLUSTERS LEVEL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  cat("\nClusters Level Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
  clusters_output <- paste0(output_path_BCR, "Clusters_analysis")
  dir.create(path = clusters_output, recursive = TRUE, showWarnings = FALSE)

  ## Quantification of unique contig analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors la quantif d'unique sera de 100% partout!!!!!!!!!!!!!!!!!!
  cat("\nQuantification analysis...\n")
  sobj <- Quantif.unique.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name)

  ## Abundance analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors l'abondance sera de 1 partout!!!!!!!!!!!!!!!!!!!!!
  cat("\nAbundance analysis...\n")
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x),patchwork::wrap_elements(scRepertoire::clonalAbundance(sobj, group.by = ident.name, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  for(x in list_type_clT) assign(paste0("plot_cluster_clonalAbundance_",x, "_T"),patchwork::wrap_elements(scRepertoire::clonalAbundance(sobj, group.by = ident.name, cloneCall = x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  ### Save
  png(paste0(clusters_output,'/clust_clonalAbundance_',sample.name,'.png'), width = 2000, height = 600)
  print(((plot_cluster_clonalAbundance_strict | plot_cluster_clonalAbundance_gene | plot_cluster_clonalAbundance_nt | plot_cluster_clonalAbundance_aa )/ (plot_cluster_clonalAbundance_strict_T | plot_cluster_clonalAbundance_gene_T | plot_cluster_clonalAbundance_nt_T | plot_cluster_clonalAbundance_aa_T )) +
    plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
  dev.off()

  ## Clonal Homeostasis analysis
  cat("\nClonal Homeostasis analysis...\n")
  sobj <- Homeo.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name)

  ## Clonal Proportions analysis
  cat("\nClonal Proportions analysis...\n")
  sobj <- Prop.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name)

  # ## Diversity analysis
  # cat("\nDiversity analysis...\n")
  # dir.create(paste0(clusters_output,"/Diversity_scores/"), recursive = TRUE, showWarnings = FALSE)
  # sobj <- Div.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = paste0(clusters_output,"/Diversity_scores/"), caption = caption, sample.name = sample.name)

  ## Frequency analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors la fréquence vaut 1 partout!
  cat("\nFrequency analysis...\n")
  Freq.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name, ident.name = ident.name, reduction=RNA.reduction, freq_col = "clonalFrequency")

  ## Clonal Overlap analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors pas d'overlap!
  cat("\nClonal Overlap analysis...\n")
  if(length(levels(Seurat::Idents(sobj)))!=1 && length(unique(sobj@meta.data[!is.na(sobj@meta.data$CTstrict),ident.name]))!=1) Overlap.c(sobj = sobj, ident.name = ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name)

  ## Physico-chemical properties of the CDR3
  cat("\nPhysico-chemical properties of the CDR3 analysis...\n")
  Physicochemical_properties.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption = caption, sample.name = sample.name, ident.name = ident.name, type='TCR')

  ## Renamme BCR columns with 'BCR_' prefix
  toMatch <- c("^CTgene","^CTnt","^CTaa","^CTstrict","^clonalProportion","^clonalFrequency","^cloneSize","^IGHV","^IGHJ","^IGHD","^Isotype","^IGLV","^IGLJ","^IGLC","^cdr3_nt_IGH", "^cdr3_nt_IGL","^cdr3_nt_IGH_len","^cdr3_nt_IGL_len","^highlight_aa")
  matches <- grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data))
  colnames(sobj@meta.data)[matches]  <- paste0("BCR_", grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data), value=TRUE))

  ## Make the results tables (raw vdj data, sobj BCR data, nb cell/clust/clone and merge of this 3 tables)
  cat("\nSaving tables...\n")
  #get raw vdj data
  cr_res[[1]] <- cr_res[[1]][,c("barcode","is_cell","high_confidence","length","chain","v_gene","d_gene","j_gene","c_gene","full_length","productive","cdr3","cdr3_nt","reads","umis")]
  # write.table(cr_res[[1]], file=paste0(output_path_BCR,"/raw_vdj.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
  #get sobj BCR data
  BCR_col = grep("^BCR_", colnames(sobj@meta.data), value = TRUE)
  BCR_col_to_remove <- grep("^BCR_highlight_aa_clust|^BCR_highlight_aa_top|^BCR_clonalProportion$|^BCR_clonalFrequency$|^BCR_cloneSize$", colnames(sobj@meta.data), value = TRUE)
  BCR_col = BCR_col[! BCR_col %in% BCR_col_to_remove]
  df_sobj <- data.frame(barcode=colnames(sobj),sobj@meta.data[,BCR_col],cluster=sobj@meta.data[[ident.name]]) #table with barcode, BCR results, and the cluster of each cell.
  write.table(df_sobj, file=paste0(output_path_BCR,"/vdj_foreach_cell.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
  #get table nb cell/cluster/clone
  df_nb_cell <- data.frame(sobj@meta.data$BCR_highlight_aa_all,sobj@meta.data[[ident.name]])
  df_nb_cell <- na.omit(df_nb_cell)
  df_nb_cell <- table(df_nb_cell)
  df_nb_cell <- data.frame(rbind(df_nb_cell))
  colnames(df_nb_cell) <- sub("X", "nbCell_byClone_Clust",colnames(df_nb_cell))
  df_nb_cell <- data.frame(BCR_highlight_aa_all=row.names(df_nb_cell),df_nb_cell)
  write.table(df_nb_cell, file=paste0(output_path_BCR,"/nb_cell_byclone_bycluster.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")
  #merge of 3 previous table
  df_merged <- merge(df_sobj, df_nb_cell, by = "BCR_highlight_aa_all",all.x=TRUE, all.y=TRUE)
  df_merged$BCR_highlight_aa_all <- NULL
  df_merged <- merge(cr_res[[1]], df_merged, by = "barcode",all.x=TRUE, all.y=FALSE)
  write.table(df_merged, file=paste0(output_path_BCR,"/vdj_merged.tsv"), quote=FALSE, row.names=FALSE, sep = "\t")

  ## Save packages versions
  sobj@misc$technical_info$scRepertoire <- utils::packageVersion('scRepertoire')
  sobj@misc$technical_info$alakazam <- utils::packageVersion('alakazam')

  ### Materials and Methods
  if(file.exists(paste0(dirname(vdj.input.file.bcr), "/Materials_and_Methods.txt"))){
    tmp <- readr::read_tsv(paste0(dirname(vdj.input.file.bcr), "/Materials_and_Methods.txt"), col_names = FALSE, show_col_types = FALSE)$X1
    tmp2 <- ""
    for (i in 1:length(tmp)) tmp2=paste(tmp2,tmp[i], sep="")
    sobj@misc$parameters$Materials_and_Methods$BCR <- tmp2
  } else sobj@misc$parameters$Materials_and_Methods$BCR <- NULL
  sobj@misc$parameters$Materials_and_Methods$BCR <- paste0(sobj@misc$parameters$Materials_and_Methods$BCR, " The annotation was merged with corresponding cell barcode of gene expression. The scRepertoire package (version ",sobj@misc$technical_info$scRepertoire,") was used to process annotation to assign clone based on Ig chains. scRepertoire allows to study contig quantification, contig abundance, contig length, clonal space homeostasis, clonal proportion, clonal overlap beetween clusters and diversity. Physicochemical properties of the CDR3, based on amino-acid sequences, was determined by the alakazam R package (version ",sobj@misc$technical_info$alakazam,").")
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
  write_MandM(sobj=sobj, output.dir=output.dir)

}else message("No BCR found!")

### Saving GE_ADT_BCR object
cat("\nSaving object...\n")
GE_BCR_file <- paste0(output.dir, basename(GE_file), '_BCR')
save(sobj, file = paste0(GE_BCR_file, '.rda'), compress = "bzip2")
