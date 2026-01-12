#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.ge", help="Input filtred seurat object (in .rda format)."),
  make_option("--output.dir.ge", help="Output path"),
  make_option("--eval.markers", help="List of genes to evaluate normalization and dimension reduction"),
  make_option("--author.name", help="Name of author of the analysis"),
  make_option("--author.mail", help="Email of author of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder"),
  ### Analysis Parameters
  # Normalization and dimension reduction
  make_option("--features.n", help="Number of High Variable Genes to consider"),
  make_option("--norm.method", help="Name of normalization method (LogNormalize or SCTransform)"),
  make_option("--HVG.FindVariableFeaturesMix", help="TRUE to user FindVariableFeaturesMix method to select HVG (after LogNormalize)"),
  make_option("--regex.genes.to.remove.from.HVG", help="Regular expression to select genes to remove of the HVG identification"),
  make_option("--dimred.method", help="Name of dimension reduction method (scbfa or bpca or pca or ica or mds)"),
  make_option("--vtr.biases.norm", help="List of biases to regress (percent_mt, percent_rb, nFeature_RNA, percent_st, Cyclone.Phase, and all other column name in metadata) into normalization (for LogNormalization or SCTransform)."),
  make_option("--vtr.biases.dimred", help="List of biases to regress (percent_mt, percent_rb, nFeature_RNA, percent_st, Cyclone.Phase, and all other column name in metadata) into dimension reduction (for scbfa or bpca)."),
  make_option("--vtr.scale", help="TRUE to center biaises to regress (for scbfa and bpca only)"),
  make_option("--dims.max", help="Number max of dimensions to compute for dimension reduction (depends on sample complexity and number of cells)"),
  make_option("--skip.eval_dims_res", help="Allow to skip the step of the evaluation of dimensions and resolutions"),
  make_option("--eval.dims.max", help="Number max of dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--eval.dims.min", help="Number min of dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--eval.dims.steps", help="Steps for dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--eval.res.max", help="Number max of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--eval.res.min", help="Number min of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--eval.res.steps", help="Steps for resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--eval.pt.size", help="Adjust point size on umap for evaluation"),
  ### Metadata
  make_option("--metadata.file", help="csv file with the metadata to add in the seurat object")
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
eval.markers <-  splitByComma.ifnotNULL(eval.markers)
list.author.name <- splitByComma.ifnotNULL(author.name)
list.author.mail <- splitByComma.ifnotNULL(author.mail)
### Computational Parameters
nthreads <- asNumeric.ifnotNULLelse(nthreads, 4)
### Analysis Parameters
# Normalization and dimension reduction
features.n <- asNumeric.ifnotNULLelse(features.n, 3000)
if (is.null(norm.method)) norm.method <- 'SCTransform'
if (is.null(HVG.FindVariableFeaturesMix)) HVG.FindVariableFeaturesMix <- FALSE
if (is.null(dimred.method)) dimred.method <- 'pca'
vtr.biases.norm <- sort(splitByComma.ifnotNULL(vtr.biases.norm))
vtr.biases.dimred <- sort(splitByComma.ifnotNULL(vtr.biases.dimred))
if (is.null(vtr.scale)) vtr.scale <- TRUE
if (vtr.scale && !(dimred.method %in% c('scbfa', 'bpca', 'mds'))) vtr.scale <- FALSE
dims.max <- asNumeric.ifnotNULLelse(dims.max, 49)
if (is.null(skip.eval_dims_res)) skip.eval_dims_res <- FALSE
eval.dims.max <- asNumeric.ifnotNULLelse(eval.dims.max, 49)
eval.dims.min <- asNumeric.ifnotNULLelse(eval.dims.min, 9)
eval.dims.steps <- asNumeric.ifnotNULLelse(eval.dims.steps, 2)
eval.res.max <- asNumeric.ifnotNULLelse(eval.res.max, 1.2)
eval.res.min <- asNumeric.ifnotNULLelse(eval.res.min, 0.1)
eval.res.steps <- asNumeric.ifnotNULLelse(eval.res.steps, 0.1)
eval.pt.size <- asNumeric.ifnotNULLelse(eval.pt.size, 2L)
### Metadata
metadata.file <- splitByComma.ifnotNULL(metadata.file)
### Fixed parameters
assay <- 'RNA'
raw.methods <- c('scbfa', 'bpca')
all.methods <- c(raw.methods, 'pca', 'mds', 'ica')
### Clean
rm(option_list,parser,args)

#### Check non-optional parameters ####
if (is.null(input.rda.ge)) stop("input.rda.ge parameter can't be empty!")
if (is.null(output.dir.ge)) stop("output.dir.ge parameter can't be empty!")

### Load data
load(input.rda.ge)

#### Get Missing Paramaters ####
sample.name.ge <- sobj@misc$params$sample.name.ge

#### Check optional parameters ####
if (!(norm.method %in% c('SCTransform','LogNormalize'))) stop("Normalization method unknown! ('LogNormalize' or 'SCTransform')")
if (norm.method == 'SCTransform' && HVG.FindVariableFeaturesMix) {
  warning("HVG.FindVariableFeaturesMix can't be used with SCTransform normalization! HVG.FindVariableFeaturesMix set to FALSE.")
  HVG.FindVariableFeaturesMix <- FALSE
}
if (!(dimred.method %in% c('pca','scbfa','bpca','mds'))) stop("Dimension Reduction method unknown! ('pca', 'scbfa', 'bpca' or 'mds')")
normalization.vtr <- vtr.biases.norm
reduction.vtr <- vtr.biases.dimred
if((dimred.method %in% c('pca', 'mds')) && !(is.null(reduction.vtr))) stop("vtr.biases.dimred can be used only with 'scbfa' or 'bpca' methods!")
if(!is.null(normalization.vtr) && !is.null(reduction.vtr)) warning("vtr.biases.norm et vtr.biases.dimred are both set!")

################################################################################
## MAIN : NORMALIZATION, DIMENSION REDUCTION AND PARAMETERS EVALUATION
################################################################################

### printing parameters
print("###########################################")
print(paste0("sample : ",sample.name.ge))
print(paste0("input.rda.ge : ",input.rda.ge))
print(paste0("output.dir.ge : ",output.dir.ge))
print("###########################################")

### Seurat multithreading
RcppParallel::setThreadOptions(numThreads = 1) #PCA, scaling
future::plan("multicore", workers = 1)  #FindMarkers, SCTransform, integration
options(future.globals.maxSize = 1000 * 1024^3)  #for a Seurat object of 1To (avoid error like: "Error: The total size of the global variables exported to the workers; for future expression (‘globals’) exceeds the maximum allowed size.")
options(future.rng.onMisuse = "ignore")       #remove useless warnings
options(SeuratCommand.umap.threads = 1) #for RunUMAP() (replace the n.thread option of the function)

### Add authors and metadata
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)
if(!is.null(metadata.file)) sobj <- add_metadata_sobj(sobj=sobj, metadata.file = metadata.file)

### Normalization and dimension reduction
cat("\nNormalization...\n")
sobj <- sc.normalization(sobj = sobj, assay = assay, normalization.method = norm.method, features.n = features.n, vtr.biases = normalization.vtr, regex.genes.to.remove.from.HVG = regex.genes.to.remove.from.HVG, HVG.FindVariableFeaturesMix = HVG.FindVariableFeaturesMix)
if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
cat("\nDimensions reduction...\n")
sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = dims.max, vtr.biases = reduction.vtr, vtr.scale = vtr.scale)

### Building reduced normalized output dir
norm_vtr = paste0(c(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr.biases[1])) paste(sobj@assays[[assay]]@misc$scaling$vtr.biases, collapse = '_') else NULL), collapse = '_')
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr.biases[1])) paste(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr.biases, collapse = '_') else NULL), collapse = '_')
norm.dim.red.dir = paste0(output.dir.ge, norm_vtr, '/', dimred_vtr)
dir.create(path = norm.dim.red.dir, recursive = TRUE, showWarnings = FALSE)

### Save packages versions
sobj@misc$technical_info$clustree <- utils::packageVersion('clustree')
sobj@misc$technical_info$patchwork <- utils::packageVersion('patchwork')

### Parameters for Materials and Methods
MM_tmp <- if(!is.null(regex.genes.to.remove.from.HVG)) paste0(" (some genes were exluded from HVG with the regex '", regex.genes.to.remove.from.HVG,"') ") else NULL
if(!is.null(normalization.vtr)){
  normalization.vtr <- stringr::str_replace(normalization.vtr, "nCount_RNA", "the number of detected transcripts")
  normalization.vtr <- stringr::str_replace(normalization.vtr, "sizeFactor", "the number of detected transcripts")
  normalization.vtr <- stringr::str_replace(normalization.vtr, "nFeature_RNA", "the number of detected genes")
  normalization.vtr <- stringr::str_replace(normalization.vtr, "percent_mt", "the proportion of mitochondrial transcripts")
  normalization.vtr <- stringr::str_replace(normalization.vtr, "percent_rb", "the proportion of ribosomal transcripts")
  normalization.vtr <- stringr::str_replace(normalization.vtr, "percent_st", "the proportion of mechanical stress response transcripts")
  normalization.vtr <- stringr::str_replace(normalization.vtr, "Cyclone.Phase", "the cell cycle phase determined by Cyclone")
  normalization.vtr <- stringr::str_replace(normalization.vtr, "Seurat.Phase", "the cell cycle phase determined by Seurat")
  MM_tmp2 <- paste0("And bias factors were regress out (",paste0(normalization.vtr, collapse = ", "),"). ")
}else {
  MM_tmp2 <- NULL
}
if(!is.null(reduction.vtr)){
  reduction.vtr <- stringr::str_replace(reduction.vtr, "nCount_RNA", "the number of detected transcripts")
  reduction.vtr <- stringr::str_replace(reduction.vtr, "sizeFactor", "the number of detected transcripts")
  reduction.vtr <- stringr::str_replace(reduction.vtr, "nFeature_RNA", "the number of detected genes")
  reduction.vtr <- stringr::str_replace(reduction.vtr, "percent_mt", "the proportion of mitochondrial transcripts")
  reduction.vtr <- stringr::str_replace(reduction.vtr, "percent_rb", "the proportion of ribosomal transcripts")
  reduction.vtr <- stringr::str_replace(reduction.vtr, "percent_st", "the proportion of mechanical stress response transcripts")
  reduction.vtr <- stringr::str_replace(reduction.vtr, "Cyclone.Phase", "the cell cycle phase determined by Cyclone")
  reduction.vtr <- stringr::str_replace(reduction.vtr, "Seurat.Phase", "the cell cycle phase determined by Seurat")
  MM_tmp3 <- paste0("Per-cell bias factors (including ", paste0(reduction.vtr, collapse = ", "),") were regressed out during the dimension reduction. ")
}else {
  MM_tmp3 <- NULL
}
if(dimred.method == 'scbfa') dimred.method.name <- 'scBFA'
if(dimred.method == 'bpca') dimred.method.name <- 'BinaryPCA'

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$part3_Norm_DimRed_Eval <- paste0("Seurat (version ",sobj@misc$technical_info$Seurat,") was applied for further data processing. ",
if(norm.method == 'SCTransform') paste0("The SCTransform normalization method (Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol. 2019;20 10.1186/s13059-019-1874-1.) was used to normalize, select ",features.n," Highly Variable Genes ", MM_tmp,"and scale them. ", MM_tmp2),
if(norm.method == 'LogNormalize') paste0("The LogNormalize normalization method from Seurat was used to transform data, then ",features.n," Highly Variable Genes, selected by ", if (HVG.FindVariableFeaturesMix) paste("FindVariableFeaturesMix() function from mixhvg R packages (version ", sobj@misc$technical_info$mixhvg,")") else "FindVariableFeatures() function ", MM_tmp,", were scaled. ", MM_tmp2),
if(dimred.method == 'pca') paste0("Person residuals were used for dimension reduction by Principal Component Analysis (PCA). "),
if(dimred.method == 'ica') paste0("Person residuals were used for dimension reduction by Independent Component Analysis (ICA). "),
if(dimred.method == 'mds') paste0("Person residuals were used for dimension reduction by Multidimensional Scaling (MDS). "),
if(dimred.method %in% raw.methods) paste0("As the ", dimred.method.name," dimension reduction method (version ",sobj@misc$technical_info$scBFA,") is meant to be applied on a subset of the count matrix, we followed the authors recommendation and applied it on the HVG. ", MM_tmp3),
if (!(skip.eval_dims_res)) "The number of dimensions to keep for further analysis was evaluated by assessing a range of reduced spaces using ",eval.dims.min," to ",eval.dims.max," dimensions, with a step of ",eval.dims.steps,". For each generated space, Louvain clustering of cells was performed using a range of values for the resolution parameter from ",eval.res.min," to ",eval.res.max," with a step of ",eval.res.steps,". The optimal space was manually evaluated as the one combination of kept dimensions and clustering resolution resolving the best structure (clusters homogeneity and compacity) in a Uniform Manifold Approximation and Projection space (UMAP). Additionaly, we used the clustree method (version ",sobj@misc$technical_info$clustree,") to assess if the selected optimal space corresponded to a relatively stable position in the clustering results tested for these dimensions / resolution combinations; and the standard deviation of each dimension was evaluated on an elbowplot."
)
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)

### Saving reduced normalized object
cat("\nSaving object...\n")
save(sobj, file = paste0(norm.dim.red.dir, '/', paste(c(sample.name.ge, norm_vtr, dimred_vtr), collapse = '_'), '.rda'), compress = "bzip2")

### Correlating reduction dimensions with biases and markers expression
cat("\nCorrelation of dimensions...\n")
dimensions.eval(sobj = sobj, reduction = paste0(assay, "_", dimred.method), eval.markers = eval.markers, slot = 'data', out.dir = norm.dim.red.dir, nthreads = floor(nthreads/2))

### Elbowplot
cat("\nElbowPlot...\n")
elb <- Seurat::ElbowPlot(sobj, ndims = dims.max, reduction = paste0(assay, "_", dimred.method)) + ggplot2::theme_classic()
ggplot2::ggsave(filename = paste0(norm.dim.red.dir,"/", sample.name.ge, '_', assay, "_", dimred.method, '_elbowplot.png'), plot = elb, width = 7, height = 4)

### Testing multiple clustering parameters (nb dims kept + Louvain resolution)
if (!(skip.eval_dims_res)){
  cat("\nEvaluation of multiple clustering parameters...\n")
  clustering.eval.mt(sobj = sobj, reduction = paste0(assay, "_", dimred.method), dimsvec = seq.int(eval.dims.min, eval.dims.max, eval.dims.steps), resvec = seq(eval.res.min,eval.res.max,eval.res.steps), out.dir = norm.dim.red.dir, solo.pt.size = eval.pt.size, nthreads = nthreads)
}
