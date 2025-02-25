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
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  # Normalization and dimension reduction
  make_option("--features.n", help="Number of High Variable Genes to consider"),
  make_option("--norm.method", help="Name of normalization method (LogNormalize or SCTransform)"),
  make_option("--dimred.method", help="Name of dimension reduction method (scbfa or bpca or pca or ica or mds)"),
  make_option("--vtr.biases", help="List of biases to regress (percent_mt, percent_rb, nFeature_RNA, percent_st, Cyclone.Phase, and all other column name in metadata)"),
  make_option("--vtr.scale", help="TRUE to center biaises to regress (for scbfa and bpca only)"),
  make_option("--dims.max", help="Number max of dimensions to compute (depends on sample complexity and number of cells)"),
  make_option("--dims.min", help="Number min of dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--dims.steps", help="Steps for dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.max", help="Number max of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.min", help="Number min of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.steps", help="Steps for resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  ### Metadata
  make_option("--metadata.file", help="csv file with the metadata to add in the seurat object"),
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
input.rda.ge <- args$options$input.rda.ge
output.dir.ge <- args$options$output.dir.ge
eval.markers <- unlist(stringr::str_split(args$options$eval.markers, ","))
list.author.name <- if (!is.null(args$options$author.name)) unlist(stringr::str_split(args$options$author.name, ","))
list.author.mail <- if (!is.null(args$options$author.mail)) unlist(stringr::str_split(args$options$author.mail, ","))
### Computational Parameters
nthreads <-  if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Normalization and dimension reduction
features.n <- if (!is.null(args$options$features.n)) as.numeric(args$options$features.n)
norm.method <- args$options$norm.method
dimred.method <- args$options$dimred.method
vtr.biases <- sort(unlist(stringr::str_split(args$options$vtr.biases, ",")))
vtr.scale <- args$options$vtr.scale
dims.max <- if (!is.null(args$options$dims.max)) as.numeric(args$options$dims.max)
dims.min <- if (!is.null(args$options$dims.min)) as.numeric(args$options$dims.min)
dims.steps <- if (!is.null(args$options$dims.steps)) as.numeric(args$options$dims.steps)
res.max <- if (!is.null(args$options$res.max)) as.numeric(args$options$res.max)
res.min <- if (!is.null(args$options$res.min)) as.numeric(args$options$res.min)
res.steps <- if (!is.null(args$options$res.steps)) as.numeric(args$options$res.steps)
### Metadata
metadata.file <-  if (!is.null(args$options$metadata.file)) unlist(stringr::str_split(args$options$metadata.file, ","))
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
    if(i %in% c("nthreads","features.n","dims.max","dims.min","dims.steps","res.max", "res.min", "res.steps")) assign(i, as.numeric(yaml_options[[i]]))else assign(i, yaml_options[[i]])
    
  }
  rm(yaml_options, i)
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

### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

### Save project parameters
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)

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
if (is.null(vtr.biases)) vtr.biases <- NULL
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
normalization.vtr <- if (norm.method == 'SCTransform') vtr.biases else NULL
reduction.vtr <- if (dimred.method %in% c('scbfa','bpca','mds')) vtr.biases else NULL
if (all(!any((!is.null(norm.method) && norm.method == 'SCTransform') || (!is.null(dimred.method) && dimred.method %in% c('scbfa', 'bpca'))) && !is.null(vtr.biases))) stop("vtr.biases can be used only with SCTransform, scbfa or bpca methods!")
if (!is.null(normalization.vtr) && !is.null(reduction.vtr)) message(paste0("Warning: vtr.biases were set in Normalisation (", norm.method, ") and Dimension reduction (", dimred.method,")!"))

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

### Add metadata
if(!is.null(metadata.file)) sobj <- add_metadata_sobj(sobj=sobj, metadata.file = metadata.file)

### Normalization and dimension reduction
cat("\nNormalization...\n")
sobj <- sc.normalization(sobj = sobj, assay = assay, normalization.method = norm.method, features.n = features.n, vtr.biases = normalization.vtr)
if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
cat("\nDimensions reduction...\n")
sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = dims.max, vtr.biases = reduction.vtr, vtr.scale = vtr.scale)

### Building reduced normalized output dir
norm_vtr = paste0(c(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr.biases[1])) paste(sobj@assays[[assay]]@misc$scaling$vtr.biases, collapse = '_') else NULL), collapse = '_')
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr.biases[1])) paste(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr.biases, collapse = '_') else NULL), collapse = '_')
norm.dim.red.dir = paste0(output.dir.ge, norm_vtr, '/', dimred_vtr)
dir.create(path = norm.dim.red.dir, recursive = TRUE, showWarnings = FALSE)

## Save packages versions
sobj@misc$technical_info$clustree <- utils::packageVersion('clustree')
sobj@misc$technical_info$patchwork <- utils::packageVersion('patchwork')

### Materials and Methods
MM_tmp <- if(dimred.method == 'pca') 'PCA' else dimred.method
if(!is.null(vtr.biases)){
  vtr.biases <- stringr::str_replace(vtr.biases, "nCount_RNA", "the number of detected transcripts")
  vtr.biases <- stringr::str_replace(vtr.biases, "sizeFactor", "the number of detected transcripts")
  vtr.biases <- stringr::str_replace(vtr.biases, "nFeature_RNA", "the number of detected genes")
  vtr.biases <- stringr::str_replace(vtr.biases, "percent_mt", "the proportion of mitochondrial transcripts")
  vtr.biases <- stringr::str_replace(vtr.biases, "percent_rb", "the proportion of ribosomal transcripts")
  vtr.biases <- stringr::str_replace(vtr.biases, "percent_st", "the proportion of mechanical stress response transcripts")
  vtr.biases <- stringr::str_replace(vtr.biases, "Cyclone.Phase", "the cell cycle phase determined by Cyclone")
  vtr.biases <- stringr::str_replace(vtr.biases, "Seurat.Phase", "the cell cycle phase determined by Seurat")
  
  MM_tmp2 <- if(norm.method == 'SCTransform') paste0(" and regress out bias factors (",paste0(vtr.biases, collapse = ", "),")") else NULL
  MM_tmp3 <- if(dimred.method == 'scbfa') paste0("Per-cell bias factors (including ", paste0(vtr.biases, collapse = ", "),") were regressed out during the scBFA dimension reduction.") else NULL
}else {
  MM_tmp2 <- NULL
  MM_tmp3 <- NULL
}
sobj@misc$parameters$Materials_and_Methods$part3_Norm_DimRed_Eval <- paste0("Seurat (version ",sobj@misc$technical_info$Seurat,") was applied for further data processing. ",
if(norm.method == 'SCTransform')  paste0("The SCTransform normalization method (Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol. 2019;20 10.1186/s13059-019-1874-1.) was used to normalize, scale, select ",features.n," Highly Variable Genes", MM_tmp2,"."),
if(norm.method == 'LogNormalize')  paste0(features.n," Highly Variable Genes (HVG) were identified using the FindVariableFeatures() method from Seurat applied on data transformed by its LogNormalize method."),
if(dimred.method == 'pca'){
  if(norm.method == 'SCTransform') paste0(" Person residuals from this regression were used for dimension reduction by Principal Component Analysis (PCA).")
  if(norm.method == 'LogNormalize') paste0(" HVG were scaled and and centered, providing person residuals used for dimension reduction by Principal Component Analysis (PCA).")
},
if(dimred.method == 'ica'){
  if(norm.method == 'SCTransform') paste0(" Person residuals from this regression were used for dimension reduction by Independent Component Analysis (ICA).")
  if(norm.method == 'LogNormalize') paste0(" HVG were scaled and and centered, providing person residuals used for dimension reduction by Independent Component Analysis (ICA).")
},
if(dimred.method == 'mds'){
  if(norm.method == 'SCTransform') paste0(" Person residuals from this regression were used for dimension reduction by Multidimensional Scaling (MDS).")
  if(norm.method == 'LogNormalize') paste0(" HVG were scaled and and centered, providing person residuals used for dimension reduction by Multidimensional Scaling (MDS).")
},
if(dimred.method == 'scbfa') paste0(" As the scBFA dimension reduction method (version ",sobj@misc$technical_info$scBFA,") is meant to be applied on a subset of the count matrix, we followed the authors recommendation and applied it on the HVG. ", MM_tmp3),
"The number of ",MM_tmp," dimensions to keep for further analysis was evaluated by assessing a range of reduced ",MM_tmp," spaces using ",dims.min," to ",dims.max," dimensions, with a step of ",dims.steps,". For each generated ",MM_tmp," space, Louvain clustering of cells was performed using a range of values for the resolution parameter from ",res.min," to ",res.max," with a step of ",res.steps,". The optimal space was manually evaluated as the one combination of kept dimensions and clustering resolution resolving the best structure (clusters homogeneity and compacity) in a Uniform Manifold Approximation and Projection space (UMAP). Additionaly, we used the clustree method (version ",sobj@misc$technical_info$clustree,") to assess if the selected optimal space corresponded to a relatively stable position in the clustering results tested for these dimensions / resolution combinations."
)
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)

### Saving reduced normalized object
cat("\nSaving object...\n")
save(sobj, file = paste0(norm.dim.red.dir, '/', paste(c(sample.name.GE, norm_vtr, dimred_vtr), collapse = '_'), '.rda'), compress = "bzip2")

### Correlating reduction dimensions with biases and markers expression
cat("\nCorrelation of dimensions...\n")
dimensions.eval(sobj = sobj, reduction = paste0(assay, "_", dimred.method), eval.markers = eval.markers, slot = 'data', out.dir = norm.dim.red.dir, nthreads = floor(nthreads/2))

### Testing multiple clustering parameters (nb dims kept + Louvain resolution)
cat("\nEvaluation of multiple clustering parameters...\n")
clustering.eval.mt(sobj = sobj, reduction = paste0(assay, "_", dimred.method), dimsvec = seq.int(dims.min, dims.max, dims.steps), resvec = seq(res.min,res.max,res.steps), out.dir = norm.dim.red.dir, solo.pt.size = 3L, BPPARAM = cl)
