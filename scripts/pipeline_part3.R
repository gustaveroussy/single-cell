#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.ge", help="Input filtred seurat object (in .rda format)."),
  make_option("--output.dir.ge", help="Output path"),
  make_option("--eval.markers", help="Genes to evaluate to check normalization and dimension reduction"),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  # Normalization and dimension reduction
  make_option("--features.n", help="Number of High Variable Genes to consider"),
  make_option("--norm.method", help="Name of normalization method (LogNormalize or SCTransform)"),
  make_option("--dimred.method", help="Name of dimension reduction method (scbfa or bpca or pca or ica or mds)"),
  make_option("--vtr", help="Liste of biases to regress (percent_mt, percent_rb, nFeature_RNA, percent_st, Cyclone.Phase, and all other column name in metadata)"),
  make_option("--vtr.scale", help="TRUE to center biaises to regress (for scbfa and bpca only)"),
  make_option("--dims.max", help="Number max of dimensions to compute (depends on sample complexity and number of cells)"),
  make_option("--dims.min", help="Number min of dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--dims.steps", help="Steps for dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.max", help="Number max of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.min", help="Number min of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.steps", help="Steps for resolution to compute for evaluation (depends on sample complexity and number of cells)"),
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
eval.markers <- unlist(stringr::str_split(args$options$eval.markers, ","))
author.name <- args$options$author.name
author.mail <- args$options$author.mail
### Computational Parameters
nthreads <-  if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Normalization and dimension reduction
features.n <- if (!is.null(args$options$features.n)) as.numeric(args$options$features.n)
norm.method <- args$options$norm.method
dimred.method <- args$options$dimred.method
vtr <- sort(unlist(stringr::str_split(args$options$vtr, ",")))
vtr.scale <- args$options$vtr.scale
dims.max <- if (!is.null(args$options$dims.max)) as.numeric(args$options$dims.max)
dims.min <- if (!is.null(args$options$dims.min)) as.numeric(args$options$dims.min)
dims.steps <- if (!is.null(args$options$dims.steps)) as.numeric(args$options$dims.steps)
res.max <- if (!is.null(args$options$res.max)) as.numeric(args$options$res.max)
res.min <- if (!is.null(args$options$res.min)) as.numeric(args$options$res.min)
res.steps <- if (!is.null(args$options$res.steps)) as.numeric(args$options$res.steps)
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

### Load data
load(input.rda.ge)

### Save project parameters
if (!is.null(author.name) && !tolower(author.name) %in% tolower(sobj@misc$params$author.name)) sobj@misc$params$author.name <- c(sobj@misc$params$author.name, author.name)
if (!is.null(author.mail) && !tolower(author.mail) %in% tolower(sobj@misc$params$author.mail)) sobj@misc$params$author.mail <- c(sobj@misc$params$author.mail, author.mail)

#### Get Missing Paramaters ####
### Project
sample.name.GE <- sobj@misc$params$sample.name.GE
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
# Normalization and dimension reduction
if (is.null(features.n)) features.n <- 3000
if (is.null(norm.method)) norm.method <- 'SCTransform'
if (is.null(dimred.method)) dimred.method <- 'pca'
if (is.null(vtr)) vtr <- NULL
if (is.null(vtr.scale)) vtr.scale <- TRUE
if (vtr.scale && !(dimred.method %in% c('scbfa', 'bpca', 'mds'))) vtr.scale <- FALSE
if (is.null(dims.max)) dims.max <- 49
if (is.null(dims.min)) dims.min <- 3
if (is.null(dims.steps)) dims.steps <- 2
if (is.null(res.max)) res.max <- 1.2
if (is.null(res.min)) res.min <- 0.1
if (is.null(res.steps)) res.steps <- 0.1
#### Fixed parameters ####
### Analysis Parameters
assay <- 'RNA'

#### Check optional parameters ####
if (!(norm.method %in% c('SCTransform','LogNormalize'))) stop("Normalization method unknown! (LogNormalize or SCTransform)")
if (!(dimred.method %in% c('pca','scbfa','bpca','mds'))) stop("Dimension Reduction method unknown! (pca, scbfa, bpca or mds)")
normalization.vtr <- if (norm.method == 'SCTransform') vtr else NULL
reduction.vtr <- if (dimred.method %in% c('scbfa','bpca','mds')) vtr else NULL
if (!(norm.method == 'SCTransform' || dimred.method %in% c('scbfa', 'bpca', 'mds')) && !is.null(vtr)) stop("vtr can be used only with SCtransform, scbfa, bpca or mds methods!")
if (!is.null(normalization.vtr) && !is.null(reduction.vtr)) message(paste0("Warning: vtr were set in Normalisation (", norm.method, ") and Dimension reduction (", dimred.method,")!"))

#### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

#########
## MAIN
#########

### E. NORMALIZATION, DIMENSION REDUCTION AND PARAMETERS EVALUATION
##-------------------

### printing parameters:
print("###########################################")
print(paste0("sample : ",sample.name.GE))
print(paste0("input.rda.ge : ",input.rda.ge))
print(paste0("output.dir.ge : ",output.dir.ge))
print("###########################################")

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Normalization and dimension reduction
cat("\nNormalization...\n")
sobj <- sc.normalization(sobj = sobj, assay = assay, normalization.method = norm.method, features.n = features.n, vtr = normalization.vtr)
if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
cat("\nDimensions reduction...\n")
sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = dims.max, vtr = reduction.vtr, vtr.scale = vtr.scale)

### Building reduced normalized output dir
norm_vtr = paste0(c(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr[1])) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL), collapse = '_')
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr[1])) paste(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
norm.dim.red.dir = paste0(output.dir.ge, norm_vtr, '/', dimred_vtr)
dir.create(path = norm.dim.red.dir, recursive = TRUE, showWarnings = FALSE)

### Materials and Methods
MM_tmp <- if(dimred.method == 'pca') 'PCA' else if(dimred.method == 'scbfa') 'scbfa'
if(!is.null(vtr)){
  vtr <- stringr::str_replace(vtr, "sizeFactor", "the number of detected transcripts")
  vtr <- stringr::str_replace(vtr, "nFeature_RNA", "the number of detected genes")
  vtr <- stringr::str_replace(vtr, "percent_mt", "the proportion of mitochondrial transcripts")
  vtr <- stringr::str_replace(vtr, "percent_rb", "the proportion of ribosomal transcripts")
  vtr <- stringr::str_replace(vtr, "percent_st", "the proportion of mechanical stress response transcripts")
  vtr <- stringr::str_replace(vtr, "Cyclone.Phase", "the cell cycle phase determined by Cyclone")
  vtr <- stringr::str_replace(vtr, "Seurat.Phase", "the cell cycle phase determined by Seurat")
  MM_tmp2 <- if(norm.method == 'SCTransform' && dimred.method == 'pca') paste0(" and regress out bias factors (",paste0(vtr, collapse = ", "),")") else if(norm.method == 'LogNormalize' && dimred.method == 'scbfa') paste0("Per-cell bias factors (including ", paste0(vtr, collapse = ", "),") were regressed out during the scBFA dimension reduction.") else NULL
}
sobj@misc$parameters$Materials_and_Methods$part3_Norm_DimRed_Eval <- paste0("Seurat (",sobj@misc$technical_info$Seurat,") was applied for further data processing. ",
if(norm.method == 'SCTransform' && dimred.method == 'pca') { paste0("The SCTransform normalization method (Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol. 2019;20 10.1186/s13059-019-1874-1.) was used to normalize, scale, select ",features.n," Highly Variable Genes", if(!is.null(vtr)) MM_tmp2,". Person residuals from this regression were used for dimension reduction by Principal Component Analysis (PCA).") },
if(norm.method == 'LogNormalize' && dimred.method == 'scbfa') { paste0("As the scBFA dimension reduction method (",sobj@misc$technical_info$scBFA,") is meant to be applied on a subset of the count matrix, we followed the authors recommendation and identified ",features.n," HVG (highly variable genes) using the FindVariableFeatures() method from Seurat applied on data transformed by its LogNormalize method. ", MM_tmp2)},
"The number of ",MM_tmp," dimensions to keep for further analysis was evaluated by assessing a range of reduced ",MM_tmp," spaces using ",dims.min," to ",dims.max," dimensions, with a step of ",dims.steps,". For each generated ",MM_tmp," space, Louvain clustering of cells was performed using a range of values for the resolution parameter from ",res.min," to ",res.max," with a step of ",res.steps,". The optimal space was manually evaluated as the one combination of kept dimensions and clustering resolution resolving the best structure (clusters homogeneity and compacity) in a Uniform Manifold Approximation and Projection space (UMAP). Additionaly, we used the clustree method (",sobj@misc$technical_info$clustree,") to assess if the selected optimal space corresponded to a relatively stable position in the clustering results tested for these dimensions / resolution combinations."
)
sobj@misc$parameters$Materials_and_Methods$packages_references <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)

### Saving reduced normalized object
cat("\nSaving object...\n")
save(sobj, file = paste0(norm.dim.red.dir, '/', paste(c(sample.name.GE, norm_vtr, dimred_vtr), collapse = '_'), '.rda'), compress = "bzip2")

### Correlating reduction dimensions with biases and markers expression
cat("\nCorrelation of dimensions...\n")
dimensions.eval(sobj = sobj, reduction = paste0(assay, "_", dimred.method), eval.markers = eval.markers, slot = 'data', out.dir = norm.dim.red.dir, nthreads = floor(nthreads/2))

### Testing multiple clustering parameters (nb dims kept + Louvain resolution)
cat("\nEvaluation of multiple clustering parameters...\n")
clustering.eval.mt(sobj = sobj, reduction = paste0(assay, "_", dimred.method), dimsvec = seq.int(dims.min, dims.max, 2), resvec = seq(res.min,res.max,.1), out.dir = norm.dim.red.dir, solo.pt.size = 3L, BPPARAM = cl)
