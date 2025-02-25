#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--sample.name.adt", help="Name of CITE-seq sample (ADT)"),
  make_option("--input.rda.ge", help="Input seurat object (in .rda format)."),
  make_option("--output.dir", help="Output path"),
  make_option("--input.dir.adt", help="Input path to the KALLISTOBUS result."),
  make_option("--author.name", help="Name of author of the analysis"),
  make_option("--author.mail", help="Email of author of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  make_option("--gene.names", help="List of gene names wich correspond to the ADT proteins."),
  make_option("--ADT.min.cutoff", help="List of quantile min to cutoff protein expression for plot (default: q30)."),
  make_option("--ADT.max.cutoff", help="List of quantile max to cutoff protein expression for plot (default: q95)."),
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
sample.name.ADT <- args$options$sample.name.adt
input.rda.ge <- args$options$input.rda.ge
output.dir <- args$options$output.dir
input.dir.adt <- args$options$input.dir.adt
list.author.name <- if (!is.null(args$options$author.name)) unlist(stringr::str_split(args$options$author.name, ","))
list.author.mail <- if (!is.null(args$options$author.mail)) unlist(stringr::str_split(args$options$author.mail, ","))
### Computational Parameters
nthreads <- as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
gene.names <- unlist(stringr::str_split(args$options$gene.names, ",")) #gene list correponding to the protein ADT
ADT.min.cutoff <- if (!is.null(args$options$ADT.min.cutoff)) unlist(stringr::str_split(args$options$ADT.min.cutoff, ","))
ADT.max.cutoff <- if (!is.null(args$options$ADT.max.cutoff)) unlist(stringr::str_split(args$options$ADT.max.cutoff, ","))
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
    if(i %in% c("nthreads")) assign(i, as.numeric(yaml_options[[i]])) else assign(i, yaml_options[[i]])
  }
  rm(yaml_options, i)
}
### Clean
rm(args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(sample.name.ADT)) stop("sample.name.ADT parameter can't be empty!")
if (is.null(input.rda.ge)) stop("input.rda.ge parameter can't be empty!")
if (is.null(output.dir)) stop("output.dir parameter can't be empty!")
if (is.null(input.dir.adt)) stop("input.dir.adt parameter can't be empty!")
if (is.null(gene.names)) stop("gene.names parameter can't be empty!")

### Load data
load(input.rda.ge)

### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

### Save project parameters
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)

#### Get Missing Paramaters ####
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
## Normalization and dimension reduction
dimred.method <- sobj@misc$params$reductions$method
norm.method_GE <- sobj@misc$params$normalization$normalization.method
assay <- if(norm.method_GE == "SCTransform") 'SCT' else 'RNA'
## Clustering
RNA.reduction <- sobj@misc$params$clustering$umap
## ADT
if (is.null(ADT.min.cutoff))  ADT.min.cutoff <- rep("q30", length(gene.names))
if (is.null(ADT.max.cutoff))  ADT.max.cutoff <- rep("q95", length(gene.names))

#### Fixed parameters ####
output_path_ADT <- paste0(output.dir, "/ADT_results/")
cor.method <- 'spearman'
norm.method_ADT <- 'LogNormalize' #gave good VISUAL results. Avoid sctransform on such small dataset
slot <- 'data' #for correlation and umap

#########
## MAIN
#########

### printing parameters:
print("###########################################")
print(paste0("sample.name.adt : ",sample.name.ADT))
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
cat("\nLoading raw count matrix of ADT...\n")
dir.create(path = output_path_ADT, recursive = TRUE, showWarnings = TRUE)
sobjADT <- load.sc.data(data.path = input.dir.adt, sample.name = sample.name.ADT, assay = 'ADT', droplets.limit = NULL, emptydrops.fdr = NULL, BPPARAM = cl, my.seed = sobj@misc$params$seed, out.dir = output_path_ADT, draw_plots = FALSE)

### Check number of protein names and gene.names and quantiles cutoff
if(length(rownames(sobjADT)) != length(gene.names)) stop(paste0("The number of gene.names is not the same as the proteins in the ADT count table: ", length(gene.names), " genes (", paste0(gene.names, collapse=","),") and ",length(rownames(sobjADT))," proteins (",paste0(rownames(sobjADT), collapse=","),")."))
if(length(rownames(sobjADT)) != length(ADT.min.cutoff)) stop(paste0("The number of ADT.min.cutoff is not the same as the proteins in the ADT count table: ", length(ADT.min.cutoff), " quantiles (", paste0(ADT.min.cutoff, collapse=","),") and ",length(rownames(sobjADT))," proteins (",paste0(rownames(sobjADT), collapse=","),")."))
if(length(rownames(sobjADT)) != length(ADT.max.cutoff)) stop(paste0("The number of ADT.max.cutoff is not the same as the proteins in the ADT count table: ", length(ADT.max.cutoff), " quantiles (", paste0(ADT.max.cutoff, collapse=","),") and ",length(rownames(sobjADT))," proteins (",paste0(rownames(sobjADT), collapse=","),")."))

### Synching ADT to GE cells
cat("\nSynching ADT to GE cells...\n")
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
cat("\nMerging ADT to GE...\n")
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
cat("\nNormalization ADT expressions...\n")
sobj <- Seurat::NormalizeData(sobj, assay = 'ADT', normalization.method = norm.method_ADT)

### Computing correlations
cat("\nComputing correlations...\n")
cor.df <- data.frame(RNA_feature = gene.names, ADT_feature = rownames(sobj@assays[['ADT']]@counts), stringsAsFactors = FALSE)
cor.unfiltered <- feature.cor(sobj = sobj, assay1 = assay, assay2 = 'ADT', assay1.features = gene.names, assay2.features = rownames(sobj@assays[['ADT']]@counts), slot = slot, cor.method = cor.method, zero.filter = FALSE, gene.names = gene.names, min.cutoff = NULL, max.cutoff = NULL)
cor.filtered <- feature.cor(sobj = sobj, assay1 = assay, assay2 = 'ADT', assay1.features = gene.names, assay2.features = rownames(sobj@assays[['ADT']]@counts), slot = slot, cor.method = cor.method, zero.filter = TRUE, gene.names = gene.names, min.cutoff = NULL, max.cutoff = NULL)
cor.quantile <- feature.cor(sobj = sobj, assay1 = assay, assay2 = 'ADT', assay1.features = gene.names, assay2.features = rownames(sobj@assays[['ADT']]@counts), slot = slot, cor.method = cor.method, zero.filter = FALSE, gene.names = gene.names, min.cutoff = ADT.min.cutoff, max.cutoff = ADT.max.cutoff)
cor.df <- cbind(cor.df, cor.unfiltered, cor.filtered,cor.quantile)
sobj@assays[['ADT']]@misc$cor <- cor.df
write.table(cor.df, file = paste0(output_path_ADT,'/ADT_correlations.csv'),sep = ";", row.names = FALSE, quote = FALSE)
rm(cor.df,cor.unfiltered,cor.filtered)

### Co-plot gene expression and ADT protein level
cat("\nCo-plot gene expression and ADT protein level...\n")
#### withtout customized cutoff
RNA_data_plot <- feature_plots(sobj, assay = assay, features = gene.names, slot = slot, reduction = RNA.reduction, min.cutoff = rep(0,length(gene.names)), max.cutoff = rep("q100",length(gene.names)))
ADT_data_plot <- feature_plots(sobj, assay = 'ADT', features = rownames(sobj@assays[['ADT']]@counts), slot = slot, reduction = RNA.reduction, min.cutoff = rep(0,length(rownames(sobj@assays[['ADT']]@counts))), max.cutoff = rep("q100",length(rownames(sobj@assays[['ADT']]@counts))))
png(paste0(output_path_ADT,'/ADT_dimplot.png'), width = 1200, height = 600 * length(gene.names))
wrap_elements(RNA_data_plot) + wrap_elements(ADT_data_plot)
dev.off()
#### with customized cutoff
cat("\nCo-plot gene expression and ADT protein level with customized cutoff...\n")
RNA_data_plot <- feature_plots(sobj, assay = assay, features = gene.names, slot = slot, reduction = RNA.reduction, min.cutoff = rep(0,length(gene.names)), max.cutoff = rep("q100",length(gene.names)))
ADT_data_plot <- feature_plots(sobj, assay = 'ADT', features = rownames(sobj@assays[['ADT']]@counts), slot = slot, reduction = RNA.reduction, min.cutoff = ADT.min.cutoff, max.cutoff = ADT.max.cutoff)
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
sobj@misc$params$ADT$gene.names <- paste0(gene.names, collapse=",")
#normalization:
sobj@assays[['ADT']]@misc$paramters$normalization$method <- norm.method_ADT
sobj@misc$params$ADT$normalization$method <- norm.method_ADT
#correlation:
sobj@assays[['ADT']]@misc$paramters$cor$method <- cor.method
sobj@assays[['ADT']]@misc$paramters$cor$slot <- slot
sobj@misc$params$ADT$cor$method <- cor.method
sobj@misc$params$ADT$cor$slot <- slot
#cutoff:
sobj@assays[['ADT']]@misc$paramters$cutoff_min = ADT.min.cutoff
sobj@assays[['ADT']]@misc$paramters$cutoff_max = ADT.max.cutoff
sobj@misc$params$ADT$cutoff_min = paste0(ADT.min.cutoff, collapse=",")
sobj@misc$params$ADT$cutoff_max = paste0(ADT.max.cutoff, collapse=",")

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$ADT <- paste0(sobj@misc$parameters$Materials_and_Methods$ADT," Only cell barcodes corresponding to the cell barcodes of gene expression were kept. Counting table was log-normalize (NormalizeData() function from Seurat with normalization.method parameters setting to '", norm.method_ADT,"') and ploted on UMAP near corresponding gene expression. ", stringr::str_to_title(cor.method)," correlation scores beetween proteins levels and genes expression levels was computed.")
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj=sobj, output.dir=output.dir)

### Saving GE_ADT object
cat("\nSaving object...\n")
GE_ADT_file <- paste0(output.dir, sub("\\.rda$|\\.RData$", "", basename(input.rda.ge)), '_ADT')
save(sobj, file = paste0(GE_ADT_file, '.rda'), compress = "bzip2")
