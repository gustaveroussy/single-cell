#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.ge", help="Input filtred, normalized and reducted seurat object (in .rda format)."),
  make_option("--output.dir.ge", help="Output path"),
  make_option("--markfile", help="Genes to plot on umap (# )format: 2 columns named Genes and Signature"),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  # Clustering
  make_option("--keep.dims", help="Number of dimension to keep for clustering (from 0 to keep.dims)"),
  make_option("--keep.res", help="Resolution value for clustering"),
  # Annotation
  make_option("--cfr.minscore", help="Minimum correlation score for clustifyr to consider"),
  make_option("--sr.minscore", help="Minimum correlation score for SingleR to consider"),
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
output.dir.ge <- args$options$output.dir.ge
markfile <- if (!is.null(args$options$markfile)) unlist(stringr::str_split(args$options$markfile, ","))
author.name <- args$options$author.name
author.mail <- args$options$author.mail
### Computational Parameters
nthreads <-  if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Clustering
keep.dims <- if (!is.null(args$options$keep.dims)) as.numeric(args$options$keep.dims)
keep.res <- if (!is.null(args$options$keep.res)) as.numeric(args$options$keep.res)
# Annotation
cfr.minscore <- if (!is.null(args$options$cfr.minscore)) as.numeric(args$options$cfr.minscore)
sr.minscore <- if (!is.null(args$options$sr.minscore)) as.numeric(args$options$sr.minscore)
### Yaml parameters file to remplace all parameters before (usefull to use R script without snakemake)
if (!is.null(args$options$yaml)){
  yaml_options <- yaml::yaml.load_file(args$options$yaml)
  for(i in names(yaml_options)) assign(i, yaml_options[[i]])
  rm(yaml_options)
}
### Clean
rm(option_list,parser,args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(input.rda.ge)) stop("input.rda.ge parameter can't be empty!")
if (is.null(output.dir.ge)) stop("output.dir.ge parameter can't be empty!")
if (is.null(keep.dims)) stop("keep.dims parameter can't be empty!")
if (is.null(keep.res)) stop("keep.res parameter can't be empty!")

### Load data
load(input.rda.ge)

### Save project parameters
if (!is.null(author.name) && !tolower(author.name) %in% tolower(sobj@misc$params$author.name)) sobj@misc$params$author.name <- c(sobj@misc$params$author.name, author.name)
if (!is.null(author.mail) && !tolower(author.mail) %in% tolower(sobj@misc$params$author.mail)) sobj@misc$params$author.mail <- c(sobj@misc$params$author.mail, author.mail)

#### Get Missing Paramaters ####
### Project
sample.name.GE <- sobj@misc$params$sample.name.GE
species <- sobj@misc$params$species
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
# Normalization and dimension reduction
assay <- Seurat::DefaultAssay(sobj)
norm.method <- sobj@assays[[assay]]@misc$params$normalization$normalization.method
dimred.method <- sobj@assays[[assay]]@misc$params$reductions$method
norm_vtr <- paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr[1])) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL, collapse = '_')
dimred_vtr <- paste0(c(dimred.method, if(!is.na(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr[1])) paste(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
# Annotation
if (is.null(cfr.minscore)) cfr.minscore <- 0.35
if (is.null(sr.minscore)) sr.minscore <- 0.25

#### Fixed parameters ####
# Annotation
if (species == "homo_sapiens") {
  singler.setnames <- c("HumanPrimaryCellAtlasData", "NovershternHematopoieticData", "DatabaseImmuneCellExpressionData", "MonacoImmuneData")
  clustifyr.setnames <- c("pbmc_avg", "hema_microarray_matrix", "gtex_bulk_matrix")
}
if (species == "mus_musculus") {
  singler.setnames <- c("MouseRNAseqData", "ImmGenData")
  clustifyr.setnames <- c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq")
}
if (species == "rattus_norvegicus") {
  singler.setnames <- c("MouseRNAseqData", "ImmGenData")
  clustifyr.setnames <- c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq")
}

#### Fixed parameters ####
# Plots
solo.pt.size <- 3
multi.pt.size <- 2
gradient.cols <- c("gold", "blue")

#### Get genes markers ####
if (is.null(markfile)){
  markers <- NULL
}else{
  markers <- c()
  for (i in markfile){
    mark.xl <- openxlsx::read.xlsx(i, sheet = 1, startRow = 1, fillMergedCells = TRUE, colNames = TRUE)
    mark.xl <- mark.xl[order(mark.xl[,2]),]
    markers_tmp <- setNames(mark.xl[,1], mark.xl[,2])
    markers <- c(markers,markers_tmp)
  }
}

#### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))


#########
## MAIN
#########

### printing parameters:
print("###########################################")
print(paste0("sample : ",sample.name.GE))
print(paste0("input.rda.ge : ",input.rda.ge))
print(paste0("output.dir.ge : ",output.dir.ge))
print(paste0("Dimension: ", keep.dims))
print(paste0("Resolution: ", keep.res))
print("###########################################")

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Building clustered output directory
clust.dir <- paste(output.dir.ge, paste0("dims", keep.dims, "_res", keep.res), sep = '/')
dir.create(path = clust.dir, recursive = TRUE, showWarnings = FALSE)

### Replotting final clusters
cat("\nClustering...\n")
sobj <- louvain.cluster(sobj = sobj, reduction = paste0(assay, "_", dimred.method), max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size)

## Setting ident name and RNA.reduction
ident.name <- paste0(paste0(assay, "_", dimred.method, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))
RNA.reduction <- paste(c(assay, dimred.method, keep.dims, 'umap'), collapse = '_')

### Technical plots
cat("\nSaving technical plots...\n")
technical.plot(sobj = sobj, ident = ident.name, out.dir = clust.dir, multi.pt.size = multi.pt.size)

### Finding markers
cat("\nFinding markers...\n")
sobj <- find.markers.quick(sobj = sobj, ident = ident.name, test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, adjp.p.max = 5E-02, topn = 10, heatmap.cols = gradient.cols, out.dir = clust.dir)

### Automatic cell type annotation
cat("\nAutomatic cell type annotation...\n")
sobj <- cells.annot(sobj = sobj, ident = ident.name, singler.setnames = singler.setnames, clustifyr.setnames = clustifyr.setnames, sr.minscore = sr.minscore, cfr.minscore = cfr.minscore, out.dir = clust.dir, solo.pt.size = solo.pt.size)

### Assessing clusters : Plotting provided marker genes
cat("\nPlotting provided marker genes...\n")
if(!is.null(markers)) sobj <- markers.umap.plot(sobj = sobj, markers = markers, ident = ident.name, out.dir = clust.dir, dimplot.cols = gradient.cols, multi.pt.size = 2)

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$part4_Clust_Markers_Annot <- paste0(
"An automatic annotation of cell types was perfom by SingleR (",sobj@misc$technical_info$SingleR,") (with fine-tuning step) and ClustifyR (",sobj@misc$technical_info$clustifyr,"), using packages built-in references. It labels clusters (or cells) from a dataset based on similarity (Spearman correlation score) to a reference dataset with known labels. The labels with a correlation score greater than ",sr.minscore," for SingleR or greater than ",cfr.minscore," for ClustifyR were kept.",
"Marker genes for Louvain clusters were identified through a «one versus others» differential anaylisis using the Wilcoxon test through the FindAllMarkers() function from Seurat, considering only genes with a minimum log fold-change of 0.5 in at least 75% of cells from one of the groups compared, and FDR-adjusted p-values <0.05 (Benjaminin-Hochberg method)."
)
sobj@misc$parameters$Materials_and_Methods$packages_references <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj=sobj, output.dir=output.dir)

### Saving final object
cat("\nSaving object...\n")
GE_file = paste0(clust.dir, '/', paste(c(sample.name.GE, norm_vtr, dimred_vtr, keep.dims, keep.res), collapse = '_'))
save(sobj, file = paste0(GE_file, '.rda'), compress = "bzip2")
