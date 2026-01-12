#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--samples.name.adt", help="Names of CITE-seq samples (ADT)"),
  make_option("--input.rda", help="Input seurat object (in .rda format)."),
  make_option("--output.dir", help="Output path"),
  make_option("--input.dirs.adt", help="Input paths to the KALLISTOBUS results."),
  make_option("--author.name", help="Name of author of the analysis"),
  make_option("--author.mail", help="Email of author of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder"),
  ### Analysis Parameters
  make_option("--gene.names", help="List of gene names wich correspond to the ADT proteins.")
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
samples.name.ADT <- splitByComma.ifnotNULL(samples.name.adt)
input.dirs.adt <- splitByComma.ifnotNULL(input.dirs.adt)
list.author.name <- splitByComma.ifnotNULL(author.name)
list.author.mail <- splitByComma.ifnotNULL(author.mail)
### Computational Parameters
nthreads <- asNumeric.ifnotNULLelse(nthreads, 1)
### Analysis Parameters
gene.names <- splitByComma.ifnotNULL(gene.names)
#### Fixed parameters
output_path_ADT <- paste0(output.dir, "/ADT_results/")
cor.method <- 'spearman'
norm.method_ADT <- 'CLR'
slot <- 'data' #for correlation and umap
### Clean
rm(args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(samples.name.ADT)) stop("samples.name.adt parameter can't be empty!")
if (is.null(input.rda)) stop("input.rda parameter can't be empty!")
if (is.null(output.dir)) stop("output.dir parameter can't be empty!")
if (is.null(input.dirs.adt)) stop("input.dirs.adt parameter can't be empty!")
if (is.null(gene.names)) stop("gene.names parameter can't be empty!")

### Load data
load(input.rda)

#### Get Missing Paramaters ####
### Analysis Parameters
dimred.method <- sobj@misc$params$reductions$method
norm.method_GE <- sobj@misc$params$normalization$normalization.method
assay <- if(norm.method_GE == "SCTransform") 'SCT' else 'RNA'
INT_GRP.reduction <- sobj@misc$params$clustering$umap
name.int_grp <- Seurat::Project(sobj)

###########################
## MAIN: ADD ADT DATA
###########################

### printing parameters:
print("###########################################")
print(paste0("samples.name.adt : ",samples.name.ADT))
print(paste0("input.rda : ",input.rda))
print(paste0("output.dir : ",output.dir))
print(paste0("input.dirs.adt : ",input.dirs.adt))
print(paste0("gene.names : ",paste0(gene.names,collapse = ", ")))
print("###########################################")

## Load libraries
require(patchwork)

## Set the seed
set.seed(sobj@misc$params$seed)

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Add authors
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)

### Loading raw count matrix
cat("\nLoading raw count matrix of ADT...\n")
dir.create(path = output_path_ADT, recursive = TRUE, showWarnings = FALSE)
sobjADT.list <- sapply(seq_along(input.dirs.adt), function(x) {
  sobjADT <- load.sc.data(data.path = input.dirs.adt[x], sample.name = samples.name.ADT[x], assay = 'ADT', droplets.limit = NULL, emptydrops.fdr = NULL, BPPARAM = cl, my.seed = sobj@misc$params$seed, out.dir = output_path_ADT, draw_plots = FALSE)
  return(sobjADT)
})
names(sobjADT.list) <- sobj@misc$params$names.ge

### Merge ADT of all samples
sobjADT <- merge(x = sobjADT.list[[1]], y = sobjADT.list[-1], add.cell.ids = names(sobjADT.list), project = name.int_grp, merge.data = TRUE)
sobjADT[["ADT"]] <- SeuratObject::JoinLayers(sobjADT[["ADT"]])
for (x in seq_along(sobjADT.list)) sobjADT@misc$pipeline_commands <- c(sobjADT@misc$pipeline_commands, sobjADT.list[[x]]@misc$pipeline_commands)
sobjADT@misc$params$ADT <- sobjADT.list[[1]]@misc$params$sobj_creation
sobjADT@misc$technical_info <- sobjADT.list[[1]]@misc$technical_info
tmp=c()
for (nb_file in 1:length(sobjADT.list)){
  tmp[nb_file] <-  if(length(sobjADT.list[[nb_file]]@misc$parameters$Materials_and_Methods$part0_Alignment) == 0) NA else sobjADT.list[[nb_file]]@misc$parameters$Materials_and_Methods$part0_Alignment
}
if(length(unique(tmp)) == 1 && !any(is.na(unique(tmp)))){
  sobjADT@misc$parameters$Materials_and_Methods$part0_Alignment <- tmp[1]
} else {
  sobjADT@misc$parameters$Materials_and_Methods$part0_Alignment <- NULL
}

### Check number of protein names and gene.names and quantiles cutoff
if(length(rownames(sobjADT)) != length(gene.names)) stop(paste0("The number of gene.names is not the same as the proteins in the ADT count table: ", length(gene.names), " genes (", paste0(gene.names, collapse=","),") and ",length(rownames(sobjADT))," proteins (",paste0(rownames(sobjADT), collapse=","),")."))

### Synching ADT to GE cells
cat("\nSynching ADT to GE cells...\n")
sobjADT <- sobjADT[, colnames(sobjADT) %in% colnames(sobj)]
if (!all(sort(colnames(sobj@assays$RNA)) == sort(colnames(sobjADT@assays$ADT)))){
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

### Merging GE/ADT
cat("\nMerging ADT to GE...\n")
sobj[['ADT']] <- Seurat::CreateAssayObject(Seurat::GetAssayData(sobjADT, assay = 'ADT', layer = "counts"))
sobj@assays[['ADT']]@misc <- sobjADT@misc
sobj@misc$pipeline_commands <- c(sobj@misc$pipeline_commands, sobjADT@misc$pipeline_commands)
sobj@misc$params$ADT <- sobjADT@misc$params$sobj_creation
sobj@misc$technical_info$BUSpaRse <- sobjADT@misc$technical_info$BUSpaRse
sobj@misc$technical_info$DropletUtils <- sobjADT@misc$technical_info$DropletUtils
if(sobj@misc$technical_info$Seurat != sobjADT@misc$technical_info$Seurat) sobj@misc$technical_info$Seurat <- paste0(sobj@misc$technical_info$Seurat, ", ", sobjADT@misc$technical_info$Seurat)
sobj@misc$parameters$Materials_and_Methods$ADT <- sobjADT@misc$parameters$Materials_and_Methods$part0_Alignment
rm(sobjADT)

### Normalization
cat("\nNormalization ADT expressions...\n")
sobj <- Seurat::NormalizeData(sobj, assay = 'ADT', normalization.method = norm.method_ADT)

### Computing correlations
cat("\nComputing correlations...\n")
cor.df <- data.frame(RNA_feature = gene.names, ADT_feature = rownames(sobj@assays$ADT), stringsAsFactors = FALSE)
cor.filtered <- feature.cor(sobj = sobj, assay1 = assay, assay2 = 'ADT', assay1.features = gene.names, assay2.features = rownames(sobj@assays$ADT), slot = slot, cor.method = cor.method, zero.filter = TRUE, gene.names = gene.names, min.cutoff = NULL, max.cutoff = NULL)
cor.df <- cbind(cor.df, cor.filtered)
for (i_min in seq(0, 100, 5)){
    for(j_max in seq(5, 100, 5)){
        if (i_min < j_max) {
            cor.quantile <- feature.cor(sobj = sobj, assay1 = assay, assay2 = 'ADT', assay1.features = gene.names, assay2.features = rownames(sobj@assays$ADT), slot = slot, cor.method = cor.method, zero.filter = FALSE, gene.names = gene.names, min.cutoff = rep(i_min,length(gene.names)), max.cutoff = rep(j_max,length(gene.names)))
            cor.df <- cbind(cor.df,cor.quantile)
        }
    }
}
sobj@assays[['ADT']]@misc$cor <- cor.df
write.table(cor.df, file = paste0(output_path_ADT,'/ADT_correlations.csv'),sep = ";", row.names = FALSE, quote = FALSE)
rm(cor.df,cor.filtered,cor.quantile)

### Co-plot gene expression and ADT protein level
cat("\nCo-plot gene expression and ADT protein level...\n")
for (k_gene in 1:length(gene.names)){
    #plot gene expression
    RNA_data_plot <- feature_plots(sobj, assay = assay, features = gene.names[k_gene], slot = slot, reduction = INT_GRP.reduction, min.cutoff = "q0", max.cutoff = "q100")
    png(paste0(output_path_ADT,'/',gene.names[k_gene], "__", rownames(sobj@assays$ADT)[k_gene],'_association_RNA_dimplot.png'), width = 600, height = 600)
    print(RNA_data_plot)
    dev.off()
    #plot protein expression
    ADT_data_plot <- list()
    for (i_min in c(0,5,seq(10, 90, 10),95)){
        for(j_max in c(seq(10, 90, 10),95, 100)){
            if (i_min < j_max) {
                #plot expression
                ADT_data_plot[[paste0("q",i_min,"_q",j_max)]] <- feature_plots(sobj, assay = 'ADT', features = rownames(sobj@assays$ADT)[k_gene], slot = slot, reduction = INT_GRP.reduction, min.cutoff = paste0("q",i_min), max.cutoff = paste0("q",j_max))
                ADT_data_plot[[paste0("q",i_min,"_q",j_max)]][[1]] <- ADT_data_plot[[paste0("q",i_min,"_q",j_max)]][[1]] + ggplot2::labs(title =paste0("q",i_min," to q",j_max))
            }else
                #plot vide
                ADT_data_plot[[paste0(i_min,"_",j_max)]] <- plot_spacer()
        }
    }
    png(paste0(output_path_ADT,'/',gene.names[k_gene], "__", rownames(sobj@assays$ADT)[k_gene],'_association_ADT_dimplot.png'), width = 6000, height = 6000)
    print(wrap_plots(ADT_data_plot) + plot_layout(ncol = 11))
    dev.off()
}

### Adding ADT expression (from @data slot) as metadata
adt.values <- t(as.matrix(sobj@assays$ADT@data))
colnames(adt.values) <- paste0("ADT_", colnames(adt.values))
sobj <-  Seurat::AddMetaData(sobj, adt.values, col.name = colnames(adt.values))
rm(adt.values)

### Save parameters
#samples names
sobj@misc$params$names.adt <- samples.name.ADT
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

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$ADT <- paste0(sobj@misc$parameters$Materials_and_Methods$ADT," Only cell barcodes corresponding to the cell barcodes of gene expression were kept. Counting table was log-normalize (NormalizeData() function from Seurat with normalization.method parameters setting to '", norm.method_ADT,"') and ", cor.method," correlation scores beetween proteins levels and genes expression levels was computed and ploted on UMAP.")
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj=sobj, output.dir=output.dir)

### Saving GE_ADT object
cat("\nSaving object...\n")
GE_ADT_file <- paste0(output.dir, sub("\\.rda$|\\.RData$", "", basename(input.rda)), '_ADT')
save(sobj, file = paste0(GE_ADT_file, '.rda'), compress = "bzip2")
