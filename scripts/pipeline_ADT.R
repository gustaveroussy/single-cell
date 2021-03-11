#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.ge", help="Input seurat object (in .rda format)."),
  make_option("--output.dir", help="Output path"),
  make_option("--input.dir.adt", help="Input path to the KALLISTOBUS result."),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  make_option("--gene.names", help="List of gene names wich correspond to the ADT proteins."),
  make_option("--ADT_min.cutoff", help="List of quantile min to cutoff protein expression for plot."),
  make_option("--ADT_max.cutoff", help="List of quantile max to cutoff protein expression for plot."),
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
input.rda.ge <- args$options$input.rda.ge
output.dir <- args$options$output.dir
input.dir.adt <- args$options$input.dir.adt
author.name <- args$options$author.name
author.mail <- args$options$author.mail
### Computational Parameters
nthreads <- as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
gene.names <- unlist(stringr::str_split(args$options$gene.names, ",")) #gene list correponding to the protein ADT
ADT_min.cutoff <- if (!is.null(args$options$ADT_min.cutoff)) unlist(stringr::str_split(args$options$ADT_min.cutoff, ","))
ADT_max.cutoff <- if (!is.null(args$options$ADT_max.cutoff)) unlist(stringr::str_split(args$options$ADT_max.cutoff, ","))
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
if (is.null(input.rda.ge)) stop("input.rda.ge parameter can't be empty!")
if (is.null(output.dir)) stop("output.dir parameter can't be empty!")
if (is.null(input.dir.adt)) stop("input.dir.adt parameter can't be empty!")
if (is.null(gene.names)) stop("gene.names parameter can't be empty!")

### Load data
load(input.rda.ge)

### Save project parameters
if (!is.null(author.name) && !tolower(author.name) %in% tolower(sobj@misc$params$author.name)) sobj@misc$params$author.name <- c(sobj@misc$params$author.name, author.name)
if (!is.null(author.mail) && !tolower(author.mail) %in% tolower(sobj@misc$params$author.mail)) sobj@misc$params$author.mail <- c(sobj@misc$params$author.mail, author.mail)

#### Get Missing Paramaters ####
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
## Normalization and dimension reduction
dimred.method <- sobj@misc$params$reductions$method
norm.method_GE <- sobj@misc$params$normalization$normalization.method
assay <- if(norm.method_GE == "SCTransform") 'SCT' else 'RNA'
## Clustering
GE_file <- sub('.rda', '', input.rda.ge)
RNA.reduction <- sobj@misc$params$clustering$umap
sample.name.ADT <- sub("_GE", "_ADT", sobj@misc$params$sample.name.GE)
## ADT
if (is.null(ADT_min.cutoff))  ADT_min.cutoff <- rep("q30", length(gene.names))
if (is.null(ADT_max.cutoff))  ADT_max.cutoff <- rep("q95", length(gene.names))

#### Fixed parameters ####
output_path_ADT <- paste0(output.dir, "/ADT_results/")
cor.method <- 'spearman'
norm.method_ADT <- 'LogNormalize' #gave good VISUAL results. Avoid sctransform on such small dataset
slot <- 'data' #for correlation and umap


#### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

#########
## MAIN
#########

### printing parameters:
print("###########################################")
print(paste0("input.rda.ge : ",input.rda.ge))
print(paste0("output.dir : ",output.dir))
print(paste0("input.dir.adt : ",input.dir.adt))
print(paste0("gene.names : ",paste0(gene.names,collapse = ", ")))
print("###########################################")

## Load libraries
require(patchwork)

## Set the seed
set.seed(sobj@misc$params$seed)

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Loading raw count matrix
dir.create(path = output_path_ADT, recursive = TRUE, showWarnings = TRUE)
sobjADT <- load.sc.data(data.path = input.dir.adt, sample.name = sample.name.ADT, assay = 'ADT', droplets.limit = NULL, emptydrops.fdr = NULL, BPPARAM = cl, my.seed = sobj@misc$params$seed, out.dir = output_path_ADT, draw_plots = FALSE)

### Check number of protein names and gene.names and quantiles cutoff
if(length(rownames(sobjADT)) != length(gene.names)) stop(paste0("The number of gene.names is not the same as the proteins in the ADT count table: ", length(gene.names), " genes (", paste0(gene.names, collapse=","),") and ",length(rownames(sobjADT))," proteins (",paste0(rownames(sobjADT), collapse=","),")."))
if(length(rownames(sobjADT)) != length(ADT_min.cutoff)) stop(paste0("The number of ADT_min.cutoff is not the same as the proteins in the ADT count table: ", length(ADT_min.cutoff), " quantiles (", paste0(ADT_min.cutoff, collapse=","),") and ",length(rownames(sobjADT))," proteins (",paste0(rownames(sobjADT), collapse=","),")."))
if(length(rownames(sobjADT)) != length(ADT_max.cutoff)) stop(paste0("The number of ADT_max.cutoff is not the same as the proteins in the ADT count table: ", length(ADT_max.cutoff), " quantiles (", paste0(ADT_max.cutoff, collapse=","),") and ",length(rownames(sobjADT))," proteins (",paste0(rownames(sobjADT), collapse=","),")."))

### Synching ADT to GE cells
sobjADT <- sobjADT[, colnames(sobjADT) %in% colnames(sobj)]
#sobj <- sobj[, colnames(sobj) %in% colnames(sobjADT)]
if (!all(sort(colnames(sobj@assays$RNA@counts)) == sort(colnames(sobjADT@assays$ADT@counts)))){
  #identification of cells which are in sobj but not in ADT
  colnames_ADT <- colnames(sobjADT)
  cells_without_ADT <- colnames(sobj[, !(colnames(sobj) %in% colnames_ADT)])
  #addition of 0 in the expression of these cells in ADT to be able to CreateAssayObject + addition of a new meta.data
  sobjADT@assays$ADT@counts <- cbind(sobjADT@assays$ADT@counts, rep(rep(0, dim(sobjADT@assays$ADT@counts)[1]),length(cells_without_ADT)))
  colnames(sobjADT@assays$ADT@counts) <- c(colnames_ADT, cells_without_ADT)
  meta.data_ADT <- data.frame(nCount_ADT = c(sobjADT@meta.data$nCount_ADT, rep(0,length(cells_without_ADT))), nFeature_ADT = c(sobjADT@meta.data$nFeature_ADT, rep(0,length(cells_without_ADT))), log_nCount_ADT = c(sobjADT@meta.data$log_nCount_ADT, rep(0,length(cells_without_ADT))), row.names=c(colnames_ADT, cells_without_ADT))
  sobj <- Seurat::AddMetaData(sobj, meta.data_ADT, col.name = c("nCount_ADT","nFeature_ADT","log_nCount_ADT"))
}else{
  sobj <- Seurat::AddMetaData(sobj, sobjADT@meta.data[,c("nCount_ADT","nFeature_ADT","log_nCount_ADT")], col.name = c("nCount_ADT","nFeature_ADT","log_nCount_ADT"))
}

### Merging
sobj[['ADT']] <- Seurat::CreateAssayObject(sobjADT@assays[['ADT']]@counts)
sobj@assays[['ADT']]@misc <- sobjADT@misc
sobj@misc$pipeline_commands <- c(sobj@misc$pipeline_commands, sobjADT@misc$pipeline_commands)
sobj@misc$params$ADT <- sobjADT@misc$params$sobj_creation
if(sobj@misc$technical_info$BUSpaRse != sobjADT@misc$technical_info$BUSpaRse) sobj@misc$technical_info$BUSpaRse <- paste0(sobj@misc$technical_info$BUSpaRse, ", ", sobjADT@misc$technical_info$BUSpaRse)
if(sobj@misc$technical_info$DropletUtils != sobjADT@misc$technical_info$DropletUtils) sobj@misc$technical_info$DropletUtils <- paste0(sobj@misc$technical_info$DropletUtils, ", ", sobjADT@misc$technical_info$DropletUtils)
if(sobj@misc$technical_info$Seurat != sobjADT@misc$technical_info$Seurat) sobj@misc$technical_info$Seurat <- paste0(sobj@misc$technical_info$Seurat, ", ", sobjADT@misc$technical_info$Seurat)
sobj@misc$parameters$Materials_and_Methods$ADT <- sobjADT@misc$parameters$Materials_and_Methods$part0_Alignment
rm(sobjADT)

### Normalization
sobj <- Seurat::NormalizeData(sobj, assay = 'ADT', normalization.method = norm.method_ADT)

### Computing correlations
cor.df <- data.frame(RNA_feature = gene.names, ADT_feature = rownames(sobj@assays[['ADT']]@counts), stringsAsFactors = FALSE)
suppressWarnings(cor.unfiltered <- feature.cor(sobj = sobj, assay1 = assay, assay2 = 'ADT', assay1.features = gene.names, assay2.features = rownames(sobj@assays[['ADT']]@counts), slot = slot, cor.method = cor.method, zero.filter = FALSE, gene.names = gene.names))
suppressWarnings(cor.filtered <- feature.cor(sobj = sobj, assay1 = assay, assay2 = 'ADT', assay1.features = gene.names, assay2.features = rownames(sobj@assays[['ADT']]@counts), slot = slot, cor.method = cor.method, zero.filter = TRUE, gene.names = gene.names))
cor.df <- cbind(cor.df, cor.unfiltered, cor.filtered)
sobj@assays[['ADT']]@misc$cor <- cor.df
rm(cor.df,cor.unfiltered,cor.filtered)

### Co-plot gene expression AND ADT protein level
#### withtout customized cutoff
RNA_data_plot <- feature_plots(sobj, assay = assay, features = gene.names, slot = slot, reduction = RNA.reduction, min.cutoff = rep(0,length(gene.names)), max.cutoff = rep("q95",length(gene.names)))
ADT_data_plot <- feature_plots(sobj, assay = 'ADT', features = rownames(sobj@assays[['ADT']]@counts), slot = slot, reduction = RNA.reduction, min.cutoff = rep(0,length(rownames(sobj@assays[['ADT']]@counts))), max.cutoff = rep("q95",length(rownames(sobj@assays[['ADT']]@counts))))
png(paste0(output_path_ADT,'/ADT_dimplot.png'), width = 1200, height = 600 * length(gene.names))
wrap_elements(RNA_data_plot) + wrap_elements(ADT_data_plot)
dev.off()
#### with customized cutoff
RNA_data_plot <- feature_plots(sobj, assay = assay, features = gene.names, slot = slot, reduction = RNA.reduction, min.cutoff = rep(0,length(gene.names)), max.cutoff = rep("q95",length(gene.names)))
ADT_data_plot <- feature_plots(sobj, assay = 'ADT', features = rownames(sobj@assays[['ADT']]@counts), slot = slot, reduction = RNA.reduction, min.cutoff = ADT_min.cutoff, max.cutoff = ADT_max.cutoff)
png(paste0(output_path_ADT,'/ADT_dimplot_legend_cutoff.png'), width = 1200, height = 600 * length(gene.names))
wrap_elements(RNA_data_plot) + wrap_elements(ADT_data_plot)
dev.off()
rm(RNA_data_plot,ADT_data_plot)

### Adding ADT expression (from @data slot) as metadata
adt.values <- t(as.matrix(sobj@assays$ADT@data))
colnames(adt.values) <- paste0("ADT_", colnames(adt.values))
sobj <-  Seurat::AddMetaData(sobj, adt.values, col.name = colnames(adt.values))
rm(adt.values)

### Save parameters
#gene.names:
sobj@assays[['ADT']]@misc$paramters$gene.names <- gene.names
sobj@misc$params$ADT$gene.names <- gene.names
#normalization:
sobj@assays[['ADT']]@misc$paramters$normalization$method <- norm.method_ADT
sobj@misc$params$ADT$normalization$method <- norm.method_ADT
#correlation:
sobj@assays[['ADT']]@misc$paramters$cor$method <- cor.method
sobj@assays[['ADT']]@misc$paramters$cor$slot <- slot
sobj@misc$params$ADT$cor$method <- cor.method
sobj@misc$params$ADT$cor$slot <- slot
#cutoff:
sobj@assays[['ADT']]@misc$paramters$cutoff_min = ADT_min.cutoff
sobj@assays[['ADT']]@misc$paramters$cutoff_max = ADT_max.cutoff
sobj@misc$params$ADT$cutoff_min = ADT_min.cutoff
sobj@misc$params$ADT$cutoff_max = ADT_max.cutoff

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$ADT <- paste0(sobj@misc$parameters$Materials_and_Methods$ADT," Only cell barcodes corresponding to the cell barcodes of gene expression were kept. Counting table was log-normalize (NormalizeData() function from Seurat with normalization.method parameters setting to '", norm.method_ADT,"') and ", cor.method," correlation scores beetween protein levels and gene expression levels was computed and ploted on UMAP.")
sobj@misc$parameters$Materials_and_Methods$packages_references <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)

### Saving GE_ADT_TCR object
GE_ADT_file <- paste0(output.dir, basename(GE_file), '_ADT')
save(sobj, file = paste0(GE_ADT_file, '.rda'), compress = "bzip2")
