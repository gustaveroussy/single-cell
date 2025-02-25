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
vdj.input.file.bcr <- args$options$vdj.input.file.bcr
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
if (is.null(vdj.input.file.bcr)) stop("vdj.input.file.bcr parameter can't be empty!")

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
GE_file <- sub("\\.rda$|\\.RData$", "", input.rda)
dimred.method <- sobj@misc$params$reductions$method
ident.name <- sobj@misc$params$clustering$ident
RNA.reduction <- sobj@misc$params$clustering$umap
sample.name <- sub("_GE", "", sobj@misc$params$sample.name.GE)

#### Fixed parameters ####
output_path_BCR <- paste0(output.dir, "/BCR_results/")
list_type_clT <- c("gene+nt", "gene", "nt", "aa")
list_type_contig <- c("nt","aa")
caption <- '“gene” - use the genes comprising the BCR
“nt” - use the nucleotide sequence of the CDR3 region
“aa” - use the amino acid sequence of the CDR3 region
“gene+nt” - use the genes comprising the BCR + the nucleotide sequence of the CDR3 region for T cells. This is the proper definition of clonotype.'

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

## Loading input data and Combining contigs
cat("\nLoading input data and Combining contigs...\n")
cr_res <- load.sc.tcr.bcr(sobj=sobj, vdj.input.file = vdj.input.file.bcr)
if(dim(cr_res)[1]!=0){ #if there are BCR
  bcr.combined <- scRepertoire::combineBCR(df = list(cr_res), samples = sample.name, ID = "BCR")

  ## GLOBAL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  cat("\nGlobal Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
  global_output <- paste0(output_path_BCR, "Global_analysis")
  dir.create(path = global_output, recursive = TRUE, showWarnings = FALSE)

  ## Quantification analysis
  cat("\nQuantification analysis...\n")
  QC.tcr.bcr(cr_res=cr_res, out.dir = global_output, type="BCR")

  ## Quantification of unique contig analysis
  Quantif.unique.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name)

  ## Abundance analysis
  cat("\nAbundance analysis...\n")
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_abundanceContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::abundanceContig(bcr.combined, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
  ### Save
  png(paste0(global_output,'/abundanceContig.png'), width = 800, height = 300)
  (plot_abundanceContig_gene_nt | plot_abundanceContig_gene | plot_abundanceContig_nt | plot_abundanceContig_aa ) +
    plot_annotation(title = sample.name, subtitle = paste0("(",dim(sobj)[2]," cells)"), caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
  dev.off()

  ## Contigs Length analysis
  cat("\nContigs Length analysis...\n")
  for(x in list_type_contig){
    ### This should give multimodal plot
    assign(paste0("plot_lengthContig_",x,"_comb"), scRepertoire::lengthContig(bcr.combined, cloneCall=x, chains = "combined") + Seurat::NoLegend())
    ### Plots the A and B chains distribution separately
    assign(paste0("plot_lengthContig_",x,"_sin"), scRepertoire::lengthContig(bcr.combined, cloneCall=x, chains = "single") + Seurat::NoLegend())
  }
  ### Save
  png(paste0(global_output,'/lengthContig.png'), width = 800, height = 800)
  (plot_lengthContig_nt_comb | plot_lengthContig_nt_sin) / (plot_lengthContig_aa_comb | plot_lengthContig_aa_sin ) +
    plot_annotation(title = sample.name, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
  dev.off()

  ## Clonal Homeostasis analysis
  cat("\nClonal Homeostasis analysis...\n")
  Homeo.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name)

  ## Clonal Proportions analysis
  cat("\nClonal Proportions analysis...\n")
  Prop.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name)

  ## Diversity analysis
  cat("\nDiversity analysis...\n")
  Div.g(combined = bcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name)

  ## Combine BCR data with seurat object
  cat("\nCombine BCR data with seurat object...\n")
  ### Corresponding barcode
  bcr.combined[[1]]$barcode <- gsub(pattern = paste0(sample.name, "_BCR_"), replacement = '', bcr.combined[[1]]$barcode)
  ### Combination
  sobj <- scRepertoire::combineExpression(df = bcr.combined, sc = sobj, cloneCall="aa")

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
  png(paste0(global_output,'/cloneType.png'), width =700, height = 3000)
  ( wrap_elements(dimplot_IGHV / dimplot_IGHJ / dimplot_IGHD / dimplot_Isotype / dimplot_IGLV / dimplot_IGLJ / dimplot_IGLC / featureplot_cdr3_nt_IGH_len / featureplot_cdr3_nt_IGL_len) +
      plot_annotation(title = 'BCR', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) )
  dev.off()

  ## Frequency analysis
  cat("\nFrequency analysis...\n")
  sobj <- Freq.g(sobj=sobj, out.dir = global_output, sample.name=sample.name, reduction=RNA.reduction, freq_col="Frequency")

  ## Physicochemical properties of the CDR3
  cat("\nPhysicochemical properties of the CDR3 analysis...\n")
  Physicochemical_properties.g(sobj=sobj, list_type_clT = list_type_clT, out.dir = global_output, sample.name=sample.name, type='BCR')





  ## CLUSTERS LEVEL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  cat("\nClusters Level Analysis >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n")
  clusters_output <- paste0(output_path_BCR, "Clusters_analysis")
  dir.create(path = clusters_output, recursive = TRUE, showWarnings = TRUE)

  ## Filter cells that have no value for the x concerned + conversion to a list by clusters
  ## Need to filter, otherwise the functions count the 'NA' as a sequence.
  for(x in list_type_clT){
    if(x=="gene+nt") y="strict" else y=x
    filtred_sobj = sobj[,!is.na(sobj@meta.data[paste0("CT", y)])]
    assign(paste0("filtred_metadata_", sub("\\+","_",x)), scRepertoire::expression2List(sc=filtred_sobj, group=ident.name))
  }
  rm(filtred_sobj)

  ## Quantification of unique contig analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors la quantif d'unique sera de 100% partout!!!!!!!!!!!!!!!!!!
  cat("\nQuantification analysis...\n")
  sobj <- Quantif.unique.c(sobj = sobj, ident.name=ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

  ## Abundance analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors l'abondance sera de 1 partout!!!!!!!!!!!!!!!!!!!!!
  cat("\nAbundance analysis...\n")
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_cluster_abundanceContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::abundanceContig(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
  ### Save
  png(paste0(clusters_output,'/abundanceContig.png'), width = 2000, height = 600)
  (plot_cluster_abundanceContig_gene_nt | plot_cluster_abundanceContig_gene | plot_cluster_abundanceContig_nt | plot_cluster_abundanceContig_aa ) +
    plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
  dev.off()

  ## Clonal Homeostasis analysis
  cat("\nClonal Homeostasis analysis...\n")
  sobj <- Homeo.c(sobj = sobj, ident.name=ident.name, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

  ## Clonal Proportions analysis
  cat("\nClonal Proportions analysis...\n")
  sobj <- Prop.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

  ## Diversity analysis
  cat("\nDiversity analysis...\n")
  sobj <- Div.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

  ## Frequency analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors la fréquence vaut 1 partout!
  cat("\nFrequency analysis...\n")
  Freq.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, ident.name=ident.name, reduction=RNA.reduction, freq_col="Frequency", filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

  ## Clonal Overlap analysis ##UTILE??? NB: il n'y a qu'un BCR par cellule alors pas d'overlap!
  cat("\nClonal Overlap analysis...\n")
  if(length(levels(Seurat::Idents(sobj)))!=1 && length(unique(sobj@meta.data[!is.na(sobj@meta.data$CTstrict),ident.name]))!=1) Overlap.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)

  ## Physico-chemical properties of the CDR3
  cat("\nPhysico-chemical properties of the CDR3 analysis...\n")
  Physicochemical_properties.c(sobj = sobj, list_type_clT = list_type_clT, out.dir = clusters_output, caption=caption, sample.name=sample.name, ident.name=ident.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt, type='BCR')

  #renamme BCR columns with 'BCR_' prefix
  toMatch <- c("^CTgene","^CTnt","^CTaa","^CTstrict","^Frequency","^cloneType","^IGHV","^IGHJ","^IGHD","^Isotype","^IGLV","^IGLJ","^IGLC","^cdr3_nt_IGH", "^cdr3_nt_IGL","^cdr3_nt_IGH_len","^cdr3_nt_IGL_len","^highlight_aa")
  matches <- grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data))
  colnames(sobj@meta.data)[matches]  <- paste0("BCR_", grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data), value=TRUE))

  ## Save packages versions
  sobj@misc$technical_info$scRepertoire <- utils::packageVersion('scRepertoire')
  sobj@misc$technical_info$alakazam <- utils::packageVersion('alakazam')

  ### Materials and Methods
  if(file.exists(paste0(dirname(vdj.input.file.bcr), "/../../Materials_and_Methods.txt"))){
    tmp <- readr::read_tsv(paste0(dirname(vdj.input.file.bcr), "/../../Materials_and_Methods.txt"), col_names = FALSE)$X1
    tmp2 <- ""
    for (i in 1:length(tmp)) tmp2=paste(tmp2,tmp[i], sep="")
    sobj@misc$parameters$Materials_and_Methods$BCR <- tmp2
  } else sobj@misc$parameters$Materials_and_Methods$BCR <- NULL
  sobj@misc$parameters$Materials_and_Methods$BCR <- paste0(sobj@misc$parameters$Materials_and_Methods$BCR, " The annotation was merged with corresponding cell barcode of gene expression. The scRepertoire package (version ",sobj@misc$technical_info$scRepertoire,") was used to process annotation to assign clonotype based on Ig chains. scRepertoire allows to study contig quantification, contig abundance, contig length, clonal space homeostasis, clonal proportion, clonal overlap beetween clusters and diversity. Physicochemical properties of the CDR3, based on amino-acid sequences, was determined by the alakazam R package (version ",sobj@misc$technical_info$alakazam,").")
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

}else message("No TCR found!")

### Saving GE_ADT_BCR object
cat("\nSaving object...\n")
GE_BCR_file <- paste0(output.dir, basename(GE_file), '_BCR')
save(sobj, file = paste0(GE_BCR_file, '.rda'), compress = "bzip2")
