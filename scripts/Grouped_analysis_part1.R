## GROUPED PROTOCOL (no integration)
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.list.rda", help="List of .rda file with individual analysis"),
  make_option("--output.dir.grp", help="Output path"),
  make_option("--name.grp", help="Name of the group (advice: include integration method name)"),
  make_option("--eval.markers", help="Genes to evaluate to check normalization and dimension reduction"),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  # Load data
  make_option("--min.cells", help="Number minimum of cells to keep a dataset"),
  # Metadata
  make_option("--metadata.file", help="csv file with the metadata to add in the seurat objects"),
  # Normalization and dimension reduction
  make_option("--keep.norm", help="If individual normalization must be kept or not (TRUE/FALSE)."),
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
input.list.rda <- unlist(stringr::str_split(args$options$input.list.rda, ","))
output.dir.grp <- args$options$output.dir.grp
name.grp <- args$options$name.grp
eval.markers <- unlist(stringr::str_split(args$options$eval.markers, ","))
list.author.name <- if (!is.null(args$options$author.name)) unlist(stringr::str_split(args$options$author.name, ","))
list.author.mail <- if (!is.null(args$options$author.mail)) unlist(stringr::str_split(args$options$author.mail, ","))
### Computational Parameters
nthreads <-  if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Load data
min.cells <- args$options$min.cells
# Metadata
metadata.file <- unlist(stringr::str_split(args$options$metadata.file, ","))
# Normalization and dimension reduction
keep.norm <- args$options$keep.norm
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
if(is.null(pipeline.path)) stop("pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(input.list.rda[1])) stop("input.list.rda parameter can't be empty!")
if (is.null(output.dir.grp)) stop("output.dir.grp parameter can't be empty!")

#### Get Missing Paramaters ####
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
# Load data
if (is.null(min.cells)) min.cells <- 0
# Normalization and dimension reduction
if (is.null(keep.norm)) keep.norm <- FALSE
if (is.null(features.n)) features.n <- 3000
if (is.null(norm.method) && !keep.norm) norm.method <- 'SCTransform'
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
solo.pt.size <- 2
raw.methods <- c('scbfa', 'bpca')
all.methods <- c(raw.methods, 'pca', 'mds', 'ica')

## Cheking parameters
if (keep.norm && !is.null(norm.method)){
  warning(paste0("To keep normalization, the norm.method parameter will not be used. Set norm.method to NULL."))
  norm.method <- NULL
}
if (!is.null(norm.method) && !(norm.method %in% c('SCTransform','LogNormalize'))) stop("Normalization method unknown! (LogNormalize or SCTransform)")
if (!is.null(dimred.method) && !(dimred.method %in% c('pca','scbfa','bpca','mds', 'Liger'))) stop("Dimension Reduction method unknown! (pca, scbfa, bpca, mds or Liger)")
normalization.vtr <- if (!is.null(norm.method) && norm.method == 'SCTransform') vtr.biases else NULL
reduction.vtr <- if (!is.null(dimred.method) && dimred.method %in% c('scbfa','bpca','mds')) vtr.biases else NULL
if (all(!any((!is.null(norm.method) && norm.method == 'SCTransform') || (!is.null(dimred.method) && dimred.method %in% c('scbfa', 'bpca'))) && !is.null(vtr.biases))) stop("vtr.biases can be used only with SCTransform, scbfa or bpca methods!")
if (!is.null(normalization.vtr) && !is.null(reduction.vtr)) warning(paste0("vtr.biases were set in Normalisation (", norm.method, ") and Dimension reduction (", dimred.method,")!"))

## Sourcing functions
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

## RUN
######

print("#########################")
print("Grouped_analysis_part1")
print(paste0("keep.norm: ",keep.norm))
print(paste0("name.grp: ",name.grp))
print("#########################")

data.path <- paste0(output.dir.grp,'/GROUPED_ANALYSIS/NO_INTEGRATED/',name.grp ,'/')
dir.create(data.path, recursive = TRUE, showWarnings = FALSE)

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Building the list of Seurat objects.
sobj.list <- sapply(seq_along(input.list.rda), function(x) {
  message(paste0("Loading '", input.list.rda[x], "' ..."))
  load(input.list.rda[x])
  ## Cleaning assays, reductions and graphs
  sobj@reductions <- list()
  sobj@graphs <- list()
  if(!keep.norm) {
    Seurat::DefaultAssay(sobj) <- "RNA"
    sobj@assays$SCT <- NULL
  }
  ### Add metadata
  if(!is.null(metadata.file)) sobj <- add_metadata_sobj(sobj=sobj, metadata.file = metadata.file)
  return(sobj)
})
names(sobj.list) <- vapply(sobj.list, Seurat::Project, 'a')
message(paste0("There are ", length(sobj.list), " samples."))
if(length(sobj.list) == 1) stop("We can't mix only one sample!") 

### Filtering low cells datasets # pour seurat surtout!
sobj.cells <- vapply(sobj.list, ncol, 1L)
if(any(sobj.cells < min.cells)) warning(paste0('Some datasets had less than ', min.cells, ' cells, thus were removed !'))
sobj.list <- sobj.list[sobj.cells >= min.cells]

## Save sample_GE names
names.GE <- names(sobj.list)

## Get species parameter
species <- sapply(seq_along(sobj.list), function(x) { sobj.list[[x]]@misc$params$species })
if(length(unique(species)) != 1) stop(paste0("We can't mix several species: ", paste0(species, collapse = ", ")))
species = species[1]

## Get assay parameter
if(keep.norm) { #keep normalisation
  n.meth <- sapply(seq_along(sobj.list), function(x) { sobj.list[[x]]@misc$params$normalization$normalization.method })
  if(length(unique(n.meth)) != 1) stop(paste0("We can't mix several normalisation method if normalization is kept: ",paste0(n.meth, collapse = ", ")))
  assay <- sobj.list[[1]]@misc$params$normalization$assay.out
  norm.method <- unique(n.meth)
}else assay <- 'RNA' # redo normalisation

### Add prefix for colnames of sample clustering and clean ADT/TCR/BCR
for (i in names(sobj.list)){
  # add prefix for colnames of sample clustering
  to_rename=grep("_res\\.",colnames(sobj.list[[i]]@meta.data), value = TRUE)
  for (j in to_rename){
    sobj.list[[i]]@meta.data[[paste0(i,'_',j)]]=sobj.list[[i]]@meta.data[[j]]
    sobj.list[[i]]@meta.data[[j]]=NULL
  }
  # cleaning sobj for ADT, TCR and BCR part
  TCR_BCR_col=grep("^ADT|^TCR|^BCR", colnames(sobj.list[[i]]@meta.data), value = TRUE)
  if(length(TCR_BCR_col) > 0) sobj.list[[i]]@meta.data[TCR_BCR_col] <- NULL
}

### If Keeping normlization, get scale.data for each sample (with SCT is Very long: average 1h per sample!)
if(keep.norm){
  cat("\nKeep Normalization...\n")
  ## Get scale.data for each sample
  names.sobj.list <- names(sobj.list)
  message("Get scale.data for each sample:")
  for (x in seq_along(sobj.list)){
    ## Scaling if necessary
    if (sum(dim(sobj.list[[x]]@assays[[assay]]@scale.data)) < 3) {
      #Check vtr.biases
      scale.vtr.all <- NULL
      if(!any(is.na(sobj.list[[x]]@assays[[assay]]@misc$scaling$vtr.biases))) {
        scale.vtr.all <- c(sobj.list[[x]]@assays[[assay]]@misc$scaling$vtr.biases)
        message(paste0("Found scaling coveriate(s) '", paste(scale.vtr.all, collapse = "', '"), "' to regress from normalization ..."))
      }
      #Scaling
      if(assay == 'SCT') {
        sobj.list[[x]] <- Seurat::ScaleData(object = sobj.list[[x]],
                                  vars.to.regress = scale.vtr.all, do.scale = FALSE, scale.max = Inf, block.size = 750)
      }else{
        sobj.list[[x]] <- Seurat::ScaleData(object = sobj.list[[x]],
                                  vars.to.regress = scale.vtr.all, do.scale = TRUE, scale.max = 10, block.size = 1000)
      }
    }
  }
  names(sobj.list) <- names.sobj.list
}

### Merge
sobj <- merge(x = sobj.list[[1]], y = sobj.list[-1], add.cell.ids = names(sobj.list), project = name.grp, merge.data = TRUE)
sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed

### Clean
rm(sobj.list)
gc()

### Complete the normalization
if(keep.norm){
  ## Get HVG (because merge delete HVG slot)
  Seurat::VariableFeatures(sobj[[assay]]) <- rownames(sobj[[assay]]@scale.data)
  sobj@assays[[assay]]@misc$params$normalization <- list(normalization.method = norm.method, assay.ori = "RNA", assay.out = assay, features.used = NA)
  sobj@misc$params$normalization$normalization.method <- norm.method
  sobj@assays[[assay]]@misc$scaling$vtr.biases <- NA
  norm_vtr = "NORMKEPT"
}else{
  ## Normalisation
  cat("\nNormalization...\n")
  sobj <- sc.normalization(sobj = sobj, assay = assay, normalization.method = norm.method, features.n = features.n, vtr.biases = normalization.vtr)
  if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
  norm_vtr = paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr.biases[1])) paste(sobj@assays[[assay]]@misc$scaling$vtr.biases, collapse = '_') else NULL)
}

### Reduction dimension
cat("\nDimensions reduction...\n")
sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = dims.max, vtr.biases = reduction.vtr, vtr.scale = vtr.scale)

### Building reduced normalized output dir
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr.biases[1])) paste(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr.biases, collapse = '_') else NULL), collapse = '_')
norm.dim.red.dir = paste0(data.path, norm_vtr, '/', dimred_vtr)
dir.create(path = norm.dim.red.dir, recursive = TRUE, showWarnings = FALSE)

## Save packages versions
sobj@misc$technical_info$clustree <- utils::packageVersion('clustree')
sobj@misc$technical_info$patchwork <- utils::packageVersion('patchwork')
sobj@misc$technical_info$Seurat <- utils::packageVersion('Seurat')

### Materials and Methods
MM_tmp <- if(dimred.method == 'pca') 'PCA' else dimred.method
if(!is.null(vtr.biases)){
  vtr.biases <- stringr::str_replace(vtr.biases, "sizeFactor", "the number of detected transcripts")
  vtr.biases <- stringr::str_replace(vtr.biases, "nFeature_RNA", "the number of detected genes")
  vtr.biases <- stringr::str_replace(vtr.biases, "percent_mt", "the proportion of mitochondrial transcripts")
  vtr.biases <- stringr::str_replace(vtr.biases, "percent_rb", "the proportion of ribosomal transcripts")
  vtr.biases <- stringr::str_replace(vtr.biases, "percent_st", "the proportion of mechanical stress response transcripts")
  vtr.biases <- stringr::str_replace(vtr.biases, "Cyclone.Phase", "the cell cycle phase determined by Cyclone")
  vtr.biases <- stringr::str_replace(vtr.biases, "Seurat.Phase", "the cell cycle phase determined by Seurat")

  MM_tmp2 <- if(norm.method == 'SCTransform') paste0(" and regress out bias factors (",paste0(vtr.biases, collapse = ", "),")") else NULL
  MM_tmp3 <- if(dimred.method == 'scbfa') paste0("Per-cell bias factors (including ", paste0(vtr.biases, collapse = ", "),") were regressed out during the scBFA dimension reduction.") else NULL
}else{
  MM_tmp2 <- NULL
  MM_tmp3 <- NULL
}
MM_tmp4 <- paste0("Seurat (version ",sobj@misc$technical_info$Seurat,") was applied for further data processing. ")
if(keep.norm){
  MM_tmp4 <- paste0(MM_tmp4, paste0("Each dataset was normalized independently by ",norm.method,", as described in the Individual Analysis section, then data were merged using the merge() function from Seurat and a common dimension reduction was performed by ",MM_tmp,"."))
}else{
 paste0("Datasets were merged using the merge() function from Seurat and a common normalization and dimension reduction were performed. ")
 if(norm.method == 'SCTransform') MM_tmp4 <- paste0(MM_tmp4, paste0("The SCTransform normalization method (Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol. 2019;20 10.1186/s13059-019-1874-1.) was used to normalize, scale, select ",features.n," Highly Variable Genes", MM_tmp2,"."))
 if(norm.method == 'LogNormalize') MM_tmp4 <- paste0(MM_tmp4, paste0(features.n," Highly Variable Genes (HVG) were identified using the FindVariableFeatures() method from Seurat applied on data transformed by its LogNormalize method."))
 if(dimred.method == 'pca'){
   if(norm.method == 'SCTransform') MM_tmp4 <- paste0(MM_tmp4, paste0(" Person residuals from this regression were used for dimension reduction by Principal Component Analysis (PCA)."))
   if(norm.method == 'LogNormalize') MM_tmp4 <- paste0(MM_tmp4, paste0(" HVG were scaled and and centered, providing person residuals used for dimension reduction by Principal Component Analysis (PCA)."))
 }
}
if(dimred.method == 'scbfa') MM_tmp4 <- paste0(MM_tmp4, paste0(" As the scBFA dimension reduction method (version ",sobj@misc$technical_info$scBFA,") is meant to be applied on a subset of the count matrix, we followed the authors recommendation and applied it on the HVG. ", MM_tmp3))
MM_tmp4 <- paste0(MM_tmp4, paste0("The number of ",MM_tmp," dimensions to keep for further analysis was evaluated by assessing a range of reduced ",MM_tmp," spaces using ",dims.min," to ",dims.max," dimensions, with a step of ",dims.steps,". For each generated ",MM_tmp," space, Louvain clustering of cells was performed using a range of values for the resolution parameter from ",res.min," to ",res.max," with a step of ",res.steps,". The optimal space was manually evaluated as the one combination of kept dimensions and clustering resolution resolving the best structure (clusters homogeneity and compacity) in a Uniform Manifold Approximation and Projection space (UMAP). Additionaly, we used the clustree method (version ",sobj@misc$technical_info$clustree,") to assess if the selected optimal space corresponded to a relatively stable position in the clustering results tested for these dimensions / resolution combinations."))
sobj@misc$parameters$Materials_and_Methods$Grouped_analysis_Norm_DimRed_Eval <- MM_tmp4
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)

### Saving reduced normalized object
cat("\nSaving object...\n")
sobj@misc$params$analysis_type <- paste0("Grouped analysis; Keep individual normalization: ", keep.norm)
sobj@misc$params$sobj_creation$Rsession <- utils::capture.output(devtools::session_info())
sobj@misc$params$species <- species
sobj@misc$params$group$keep.norm <- keep.norm
sobj@misc$params$name.grp <- name.grp
sobj@misc$params$names.ge <- names.GE
Seurat::Project(sobj) <- name.grp
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)
save(sobj, file = paste0(norm.dim.red.dir, '/', paste(c(name.grp, norm_vtr, dimred_vtr), collapse = '_'), '.rda'), compress = "bzip2")

### Correlating reduction dimensions with biases and markers expression
cat("\nCorrelation of dimensions...\n")
dimensions.eval(sobj = sobj, reduction = paste0(assay, "_", dimred.method), eval.markers = eval.markers, slot = 'data', out.dir = norm.dim.red.dir, nthreads = floor(nthreads/2))
gc()

### Testing multiple clustering parameters (nb dims kept + Louvain resolution)
cat("\nEvaluation of multiple clustering parameters...\n")
clustering.eval.mt(sobj = sobj, reduction = paste0(assay, "_", dimred.method), dimsvec = seq.int(dims.min, dims.max, dims.steps), resvec = seq(res.min,res.max,res.steps), out.dir = norm.dim.red.dir, solo.pt.size = solo.pt.size, BPPARAM = cl)
