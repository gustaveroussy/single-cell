## INTEGRATION PROTOCOL (Seurat or scbfa)
## To Do: test Harmony and LIGER integration
#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.list.rda", help="List of .rda file with individual analysis"),
  make_option("--output.dir.int", help="Output path"),
  make_option("--name.int", help="Name of the integration (advice: include integration method name)"),
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
  # Integration
  make_option("--integration.method", help="Name of integration method (Seurat, scbfa, Harmony or Liger)"),
  make_option("--vtr.batch", help="List of batch effect to regress into Harmony, Liger, scbfa or bpca correction ('orig.ident')"),
  # Normalization and dimension reduction
  make_option("--features.n", help="Number of High Variable Genes to consider"),
  make_option("--norm.method", help="Name of normalization method (LogNormalize or SCTransform)"),
  make_option("--dimred.method", help="Name of dimension reduction method (scbfa or bpca or pca or ica or mds or Liger)"),
  make_option("--vtr.biases", help="List of biases to regress into normalisation or dimension reduction (percent_mt, percent_rb, nFeature_RNA, percent_st, Cyclone.Phase, and all other column name in metadata)"),
  make_option("--vtr.scale", help="TRUE to center biaises to regress (for scbfa and bpca only)"),
  make_option("--dims.max", help="Number max of dimensions to compute (depends on sample complexity and number of cells)"),
  make_option("--dims.min", help="Number min of dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--dims.steps", help="Steps for dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.max", help="Number max of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.min", help="Number min of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.steps", help="Steps for resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  #### Yaml parameters file to remplace all parameters before (to use R script without snakemake)
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
output.dir.int <- args$options$output.dir.int
name.int <- args$options$name.int
eval.markers <- unlist(stringr::str_split(args$options$eval.markers, ","))
list.author.name <- unlist(stringr::str_split(args$options$author.name, ","))
list.author.mail <- unlist(stringr::str_split(args$options$author.mail, ","))
### Computational Parameters
nthreads <- if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads) else NULL
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Load data
min.cells <- if (!is.null(args$options$min.cells)) as.numeric(args$options$min.cells) else NULL
# Metadata
metadata.file <- unlist(stringr::str_split(args$options$metadata.file, ","))
# Integration
integration.method <- args$options$integration.method
vtr.batch <- unlist(stringr::str_split(args$options$vtr.batch, ","))
# Normalization and dimension reduction
features.n <- if (!is.null(args$options$features.n)) as.numeric(args$options$features.n) else NULL
norm.method <- args$options$norm.method
dimred.method <- args$options$dimred.method
vtr.biases <- sort(unlist(stringr::str_split(args$options$vtr.biases, ",")))
vtr.scale <- args$options$vtr.scale
dims.max <- if (!is.null(args$options$dims.max)) as.numeric(args$options$dims.max) else NULL
dims.min <- if (!is.null(args$options$dims.min)) as.numeric(args$options$dims.min) else NULL
dims.steps <- if (!is.null(args$options$dims.steps)) as.numeric(args$options$dims.steps) else NULL
res.max <- if (!is.null(args$options$res.max)) as.numeric(args$options$res.max) else NULL
res.min <- if (!is.null(args$options$res.min)) as.numeric(args$options$res.min) else NULL
res.steps <- if (!is.null(args$options$res.steps)) as.numeric(args$options$res.steps) else NULL

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
    if(i %in% c("nthreads","min.cells", "features.n", "dims.max", "dims.min", "dims.steps", "res.max", "res.min", "res.steps")) assign(i, as.numeric(yaml_options[[i]])) else assign(i, yaml_options[[i]])
  }
  name.int <- name.int
  rm(yaml_options, i)
}
### Clean
rm(option_list,parser,args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Get Missing Paramaters ####
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
# Load data
if (is.null(min.cells)) min.cells <- 0
# Normalization and dimension reduction
if (is.null(features.n)) features.n <- 3000
if (is.null(vtr.scale)) vtr.scale <- TRUE
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

#### Check non-optional parameters ####
if (is.null(input.list.rda)) stop("input.list.rda parameter can't be empty!")
if (is.null(output.dir.int)) stop("output.dir.int parameter can't be empty!")
if (is.null(integration.method)) stop("integration.method parameter can't be empty!")

## Cheking parameters
if (!(integration.method %in% c("Seurat", "scbfa", "Harmony", "Liger"))) stop("Integration method unknown! (Seurat, scbfa, Harmony or Liger)")
if (((integration.method == "Seurat" )||(integration.method == "Liger")) && !is.null(norm.method)){
  warning(paste0("With ",integration.method," integration, normalization of each sample is kept. Set norm.method to NULL."))
  norm.method <- NULL
}
if (integration.method == "Liger"){
  if(!(dimred.method == "Liger")){
    warning(paste0("With ",integration.method," integration, dimensions reduction is included. Set dimred.method to Liger"))
    dimred.method <- "Liger"
  }
}

if (is.null(norm.method) && ((integration.method == "Harmony") || (integration.method == "scbfa"))) norm.method <- 'SCTransform'
if (is.null(dimred.method) && ((integration.method == "Seurat") || (integration.method == "Harmony"))) dimred.method <- 'pca'
if(is.null(dimred.method) && (integration.method %in% raw.methods)) dimred.method <- integration.method
if((integration.method != dimred.method) && (integration.method %in% raw.methods)){
  warning(paste("For ", integration.method, " integration, dimensions reduction is included. Set dimred.method to ", integration.method, "."))
  dimred.method <- integration.method
}
if (is.null(vtr.batch) && (integration.method %in% c("Harmony","Liger","scbfa","bpca"))) stop(paste0("vtr.batch parameter can't be empty for ",integration.method,"! (example: orig.ident)"))
if ((integration.method == "Seurat") && !is.null(vtr.batch)){
  warning(paste0("With ",integration.method," integration, vtr.batch is not used. Set vtr.batch to NULL."))
  vtr.batch <- NULL
}

if (!is.null(norm.method) && !(norm.method %in% c('SCTransform','LogNormalize'))) stop("Normalization method unknown! (LogNormalize or SCTransform)")
if (!is.null(dimred.method) && !(dimred.method %in% c('pca','scbfa','bpca','mds', 'Liger'))) stop("Dimension Reduction method unknown! (pca, scbfa, bpca, mds or Liger)")
normalization.vtr <- if (!is.null(norm.method) && norm.method == 'SCTransform') vtr.biases else NULL
reduction.vtr <- if (!is.null(dimred.method) && dimred.method %in% c('scbfa','bpca')) vtr.biases else NULL
if (all(!any((!is.null(norm.method) && norm.method == 'SCTransform') || (!is.null(dimred.method) && dimred.method %in% c('scbfa', 'bpca'))) && !is.null(vtr.biases))) stop("vtr.biases can be used only with SCTransform, scbfa or bpca methods!")
if (!is.null(normalization.vtr) && !is.null(reduction.vtr)) warning(paste0("vtr.biases were set in Normalisation (", norm.method, ") and Dimension reduction (", dimred.method,")!"))
if (vtr.scale && !(dimred.method %in% c('scbfa', 'bpca'))){
  message(paste0("vtr.scale is used only for scbfa or bpca. Set vtr.scale to FALSE."))
  vtr.scale <- FALSE
}

## Sourcing functions
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

## RUN
######

print("#########################")
print(paste0("integration.method: ",integration.method))
print(paste0("name.int: ",name.int))
if (!is.null(norm.method)) print(paste0("norm.method: ",norm.method))
if (!is.null(dimred.method)) print(paste0("dimred.method: ",dimred.method))
print("#########################")

data.path <- paste0(output.dir.int,'/GROUPED_ANALYSIS/INTEGRATED/',name.int ,'/')
dir.create(data.path, recursive = TRUE, showWarnings = FALSE)

## Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

## Building the list of Seurat objects.
sobj.list <- sapply(seq_along(input.list.rda), function(x) {
  message(paste0("Loading '", input.list.rda[x], "' ..."))
  load(input.list.rda[x])
  ## Cleaning reductions and graphs
  sobj@reductions <- list()
  sobj@graphs <- list()
  ## Add metadata
  if(!is.null(metadata.file)) sobj <- add_metadata_sobj(sobj=sobj, metadata.file = metadata.file)
  return(sobj)
})
names(sobj.list) <- vapply(sobj.list, Seurat::Project, 'a')
message(paste0("There are ", length(sobj.list), " samples."))
if(length(sobj.list) == 1) stop("We can't mix only one sample!") 

## Filtering low cells datasets # pour seurat surtout!
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
if(integration.method  %in% c("scbfa", "Harmony")){ # redo normalisation
  assay <- 'RNA'
}
if(integration.method %in% c("Seurat","Liger")) { #keep normalisation
  n.meth <- sapply(seq_along(sobj.list), function(x) { sobj.list[[x]]@misc$params$normalization$normalization.method })
  if(length(unique(n.meth)) != 1) stop(paste0("We can't mix several normalisation method with ",integration.method,": ",paste0(n.meth, collapse = ", ")))
  assay <- sobj.list[[1]]@misc$params$normalization$assay.out
  norm.method <- unique(n.meth)
}

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

## Seurat Integration
if((integration.method == "Seurat") || (integration.method == 'Liger')){
  ## Get scale.data for each sample
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
      Seurat::DefaultAssay(sobj.list[[x]]) <- assay
      if(assay == 'SCT') {
        sobj.list[[x]] <- Seurat::ScaleData(object = sobj.list[[x]],
                                  vars.to.regress = scale.vtr.all,
                                  do.scale = FALSE, scale.max = Inf, block.size = 750)
      } else {
        sobj.list[[x]] <- Seurat::ScaleData(object = sobj.list[[x]],
                                  vars.to.regress = scale.vtr.all,
                                  do.scale = if(integration.method == 'Liger') FALSE else TRUE, scale.max = 10, block.size = 1000)
      }
    }
  }
  ## Integration
  if(integration.method == "Seurat"){
    message(paste0(integration.method," integration..."))
    sobj.list <- sapply(seq_along(sobj.list), function(x) { #Rename cells
      new_cells_name <- paste0(names(sobj.list)[x],"_",colnames(sobj.list[[x]]))
      sobj.list[[x]] <- Seurat::RenameCells(sobj.list[[x]], new.names = new_cells_name)
      return(sobj.list[[x]])
    })
    if(assay == 'SCT') int.norm.method <- 'SCT' else  int.norm.method <- 'LogNormalize'
    sobj.features <- Seurat::SelectIntegrationFeatures(object.list = sobj.list, nfeatures = features.n) #Sélection des marqueurs biologiques partagés
    sobj.list <- suppressWarnings(Seurat::PrepSCTIntegration(object.list = sobj.list, anchor.features = sobj.features)) #Verifie que les résidus de Pearson ont bien été calculés
    sobj.anchors <- Seurat::FindIntegrationAnchors(object.list = sobj.list, normalization.method = int.norm.method, anchor.features = sobj.features) #CCA + L2normalisation; puis KNN; puis MNNs : identification des paires de cellules; filtrage des anchors; calcul des scores
    sobj <- Seurat::IntegrateData(anchorset = sobj.anchors, normalization.method = int.norm.method) #Calcul des poids; application des poids sur la matrice d'expression: intégration

    # Params
    
    Seurat::Project(sobj) <- name.int
    sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed
    DefaultAssay(sobj) <- "integrated"
    sobj@assays[["integrated"]]@misc$scaling$vtr.biases = NA
    sobj@assays[["integrated"]]@misc$params$normalization <- list(normalization.method = norm.method, assay.ori = NA, assay.out = NA, features.used = NA)
    sobj@misc$params$normalization$normalization.method <- norm.method
    sobj@misc$params$integration$method <- integration.method
    sobj@misc$params$integration$orig.assay <- assay
    sobj@misc$params$integration$out.assay <- "integrated"
    sobj@misc$params$integration$vtr.batch <- NA
    assay <- "integrated"
    # Cleaning
    rm(sobj.list, sobj.features, sobj.anchors)
    gc()
    ## Dimensions reduction
    cat("\nDimensions reduction...\n")
    red.name <- paste(c("integrated", dimred.method, integration.method), collapse = '_')
    sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = "integrated", max.dims = dims.max, vtr.biases = reduction.vtr, vtr.scale = vtr.scale, red.name = red.name)

  } else if(integration.method == 'Liger'){
    ## Merge data
    cat("\nMerge data...\n")
    sobj <- merge(x = sobj.list[[1]], y = sobj.list[-1], add.cell.ids = names(sobj.list), project = name.int, merge.data = TRUE)
    Seurat::Project(sobj) <- name.int
    sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed
    ## Cleaning
    rm(sobj.list)
    gc()
    ## Integration LIGER
    message(paste0(integration.method," integration..."))
    red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
    sobj <- SeuratWrappers::RunOptimizeALS(sobj, k = dims.max, lambda = 5, split.by = vtr.batch, rand.seed = sobj@misc$params$seed) #calcul les matrices (fait la red dim)
    sobj <- SeuratWrappers::RunQuantileNorm(sobj, resolution = 1, split.by = vtr.batch, reduction.name = red.name, rand.seed = sobj@misc$params$seed) #fait le SFN graph + clustering + correction de la red dim
    # Params
    sobj@assays[[assay]]@misc$scaling$vtr.biases[1] <- NA
    sobj@reductions[[red.name]]@misc$vtr.biases <- vtr.batch
    sobj@reductions[[red.name]]@misc$from.assay <- assay
    sobj@assays[[assay]]@misc$params$normalization <- list(normalization.method = norm.method, assay.ori = NA, assay.out = NA, features.used = NA)
    sobj@assays[[assay]]@misc$params$reductions <- list(method = integration.method, assay = assay, max.dims = dims.max, vtr.biases = vtr.batch, scale_regressed_variables = NA)
    sobj@misc$params$integration$method <- integration.method
    sobj@misc$params$integration$orig.assay <- assay
    sobj@misc$params$integration$out.assay <- assay
    sobj@misc$params$integration$vtr.batch <- vtr.batch
    sobj@misc$technical_info$SeuratWrappers <- utils::packageVersion('SeuratWrappers')
    # Cleaning
    sobj@assays[[assay]]@scale.data <- matrix(nrow = 0, ncol = 0)
    gc()
  }
}

## scbfa/bpca (or Harmony integration beging)
if((integration.method %in% raw.methods) || (integration.method == 'Harmony')){
  ## Merge data
  cat("\nMerge data...\n")
  sobj <- merge(x = sobj.list[[1]], y = sobj.list[-1], add.cell.ids = names(sobj.list), project = name.int, merge.data = TRUE)
  Seurat::Project(sobj) <- name.int
  sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed
  ## Cleaning
  rm(sobj.list)
  gc()
  #Normalization
  cat("\nNormalization...\n")
  sobj <- sc.normalization(sobj = sobj, assay = assay, normalization.method = norm.method, features.n = features.n, vtr.biases = normalization.vtr)
  if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
  ## Dimensions reduction
  cat("\nDimensions reduction...\n")
  if(integration.method %in% raw.methods) message(paste0(integration.method," integration..."))
  red.name <- if(integration.method %in% raw.methods) paste(c(assay, dimred.method, integration.method), collapse = '_') else paste(c(assay, dimred.method), collapse = '_')
  sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = dims.max, vtr.biases = if(integration.method == 'Harmony') reduction.vtr else c(reduction.vtr, vtr.batch), vtr.scale = vtr.scale, red.name = red.name)
  # Params
  sobj@misc$params$integration$method <- integration.method
  sobj@misc$params$integration$orig.assay <- "RNA"
  sobj@misc$params$integration$out.assay <- assay
  sobj@misc$params$integration$vtr.batch <- vtr.batch
}

### Building reduced normalized output dir
if (integration.method %in% c("Seurat", "Liger")){
  norm_vtr = "NORMKEPT"
}else{
  norm_vtr = paste0(c(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr.biases[1])) paste(sobj@assays[[assay]]@misc$scaling$vtr.biases, collapse = '_') else NULL), collapse = '_')
}
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[red.name]]@misc$vtr.biases[1])) paste(sobj@reductions[[red.name]]@misc$vtr.biases, collapse = '_') else NULL), collapse = '_')
norm.dim.red.dir = paste0(data.path, norm_vtr, '/', dimred_vtr)
dir.create(path = norm.dim.red.dir, recursive = TRUE, showWarnings = FALSE)

## Harmony integration ending
if (integration.method == 'Harmony'){
  message(paste0(integration.method," integration..."))
  ## Scaling if necessary #NB: harmony needs scale.data
  if (sum(dim(sobj@assays[[assay]]@scale.data)) < 3) {
    #Check vtr.biases
    scale.vtr.all <- NULL
    if(!any(is.na(sobj@assays[[assay]]@misc$scaling$vtr.biases))) {
      scale.vtr.all <- c(sobj@assays[[assay]]@misc$scaling$vtr.biases)
      message(paste0("Found scaling coveriate(s) '", paste(scale.vtr.all, collapse = "', '"), "' to regress from normalization ..."))
    }
    #Scaling
    Seurat::DefaultAssay(sobj) <- assay
    if(assay == 'SCT') {
      sobj <- Seurat::ScaleData(object = sobj,
                                vars.to.regress = scale.vtr.all,
                                do.scale = FALSE, scale.max = Inf, block.size = 750)
    } else {
      sobj <- Seurat::ScaleData(object = sobj,
                                vars.to.regress = scale.vtr.all,
                                do.scale = TRUE, scale.max = 10, block.size = 1000)
    }
  }
  ## Integration
  red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
  png(paste0(norm.dim.red.dir, '/harmony_convergence_plot.png'), width = 1000, height = 1000)
  sobj <- harmony::RunHarmony(sobj, vtr.batch, reduction = paste(c(assay, dimred.method), collapse = '_'), assay.use = assay, plot_convergence = TRUE, reduction.save = red.name) #, do_pca=FALSE ??
  dev.off()
  ##Clean
  sobj@assays[[assay]]@scale.data <- matrix(nrow = 0, ncol = 0)
  gc()
  #Params
  sobj@misc$technical_info$Harmony <- utils::packageVersion('harmony')
  sobj@reductions[[red.name]]@misc$vtr.biases = NA
  sobj@reductions[[red.name]]@misc$from.assay <- assay
}

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
}else {
  MM_tmp2 <- NULL
  MM_tmp3 <- NULL
}
if(integration.method == "Seurat"){
  sobj@misc$parameters$Materials_and_Methods$Integration_Norm_DimRed_Eval <- paste0(
    "Datasets were integrated using canonical correlation analysis (CCA) to identify pairwise anchors between datasets and using the anchors to harmonize the datasets, as implemented in Seurat (version ",sobj@misc$technical_info$Seurat,").",
    "In practice, datasets were normalized independently using ",norm.method,", as described before (see Individual Analysis section). ",
    "The top ",features.n," Highly Variable Genes across all samples were selected by the SelectIntegrationFeatures() function and used for integration with the PrepSCTIntegration(), FindIntegrationAnchors() and IntegrateData() functions (with default parameters). ",
    if(dimred.method == 'pca') paste("After integration, a ",MM_tmp," dimension reduction was performed on integrated dataset."),
    if(dimred.method == 'scbfa') paste("After integration, a ",dimred.method," dimension reduction was performed on HVG of integrated dataset.", MM_tmp3)
  )
}else if(integration.method == "scbfa"){
  sobj@misc$parameters$Materials_and_Methods$Integration_Norm_DimRed_Eval <- paste0("Datasets were integrated using the scBFA dimension reduction method. Datasets were merged by the merge() function from Seurat (version ",sobj@misc$technical_info$Seurat,"), and ",
    if(norm.method == 'SCTransform')  paste0(" the SCTransform normalization method (Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol. 2019;20 10.1186/s13059-019-1874-1.) was used to normalize, scale, select ",features.n," Highly Variable Genes", MM_tmp2,"."),
    if(norm.method == 'LogNormalize')  paste0(features.n," Highly Variable Genes (HVG) were identified using the FindVariableFeatures() method from Seurat applied on data transformed by its LogNormalize method."),
    if(dimred.method == 'scbfa') paste0(" As the scBFA dimension reduction method (version ",sobj@misc$technical_info$scBFA,") is meant to be applied on a subset of the count matrix, we followed the authors recommendation and applied it on the HVG. ", MM_tmp3, "The batch effect (", paste0(sapply(vtr.batch, function(x) {return(paste0(x,": ",paste0(unique(sort(sobj@meta.data[[x]])),collapse=", ")))}), collapse="; "), ") was regressed with the other potential biases.")
  )
}else if(integration.method == "Liger"){
  sobj@misc$parameters$Materials_and_Methods$Integration_Norm_DimRed_Eval <- paste0("Datasets were integrated using the Liger method. Datasets were individually normalized by ",norm.method,", as described before. ",
    "Data were merged by the merge() function from Seurat (version ",sobj@misc$technical_info$Seurat,"), and the integration was performed by the Liger method included in SeuratWrappers package (version ",sobj@misc$technical_info$SeuratWrappers,"). The functions used were RunOptimizeALS() to make dimension reduction and RunQuantileNorm() to build a shared factor neighborhood graph to jointly cluster cells and corrects this clusters, with the batch effect (", paste0(sapply(vtr.batch, function(x) {return(paste0(x,": ",paste0(unique(sort(sobj@meta.data[[x]])),collapse=", ")))}), collapse="; "), ") for split.by parameter and default for other parameters."
  )
}else if(integration.method == "Harmony"){
  sobj@misc$parameters$Materials_and_Methods$Integration_Norm_DimRed_Eval <- paste0("Datasets were integrated using the Harmony method. Datasets were merged by the merge() function from Seurat(version ",sobj@misc$technical_info$Seurat,"), and",
    if(norm.method == 'SCTransform')  paste0("the SCTransform normalization method (Hafemeister C, Satija R. Normalization and variance stabilization of single-cell RNA-seq data using regularized negative binomial regression. Genome Biol. 2019;20 10.1186/s13059-019-1874-1.) was used to normalize, scale, select ",features.n," Highly Variable Genes", MM_tmp2,"."),
    if(norm.method == 'LogNormalize')  paste0(features.n," Highly Variable Genes (HVG) were identified using the FindVariableFeatures() method from Seurat applied on data transformed by its LogNormalize method."),
    if(dimred.method == 'pca'){
      if(norm.method == 'SCTransform') paste0(" Person residuals from this regression were used for dimension reduction by Principal Component Analysis (PCA).")
      if(norm.method == 'LogNormalize') paste0(" HVG were scaled and and centered, providing person residuals used for dimension reduction by Principal Component Analysis (PCA).")
    },
    if(dimred.method == 'scbfa') paste0(" As the scBFA dimension reduction method (version ",sobj@misc$technical_info$scBFA,") is meant to be applied on a subset of the count matrix, we followed the authors recommendation and applied it on the HVG. ", MM_tmp3),
    " The reduced ",MM_tmp," spaces are used as input for the HarmonyMatrix() function implemented in Harmony package (version ",sobj@misc$technical_info$Harmony,") where the batch effect (", paste0(sapply(vtr.batch, function(x) {return(paste0(x,": ",paste0(unique(sort(sobj@meta.data[[x]])),collapse=", ")))}), collapse="; "), ") was regressed. The batch-corrected shared space output by harmony is used for clustering."
  )
  MM_tmp <- "Harmony"
}
sobj@misc$parameters$Materials_and_Methods$Integration_Norm_DimRed_Eval <- paste0( sobj@misc$parameters$Materials_and_Methods$Integration_Norm_DimRed_Eval, "The number of ",MM_tmp," dimensions to keep for further analysis was evaluated by assessing a range of reduced ",MM_tmp," spaces using ",dims.min," to ",dims.max," dimensions, with a step of ",dims.steps,". For each generated ",MM_tmp," space, Louvain clustering of cells was performed using a range of values for the resolution parameter from ",res.min," to ",res.max," with a step of ",res.steps,". The optimal space was manually evaluated as the one combination of kept dimensions and clustering resolution resolving the best structure (clusters homogeneity and compacity) in a Uniform Manifold Approximation and Projection space (UMAP). Additionaly, we used the clustree method (version ",sobj@misc$technical_info$clustree,") to assess if the selected optimal space corresponded to a relatively stable position in the clustering results tested for these dimensions / resolution combinations.")
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)

### Saving reduced normalized object
cat("\nSaving object...\n")
sobj@misc$params$analysis_type <- paste0("Integrated analysis; Method: ", integration.method)
sobj@misc$params$sobj_creation$Rsession <- utils::capture.output(devtools::session_info())
sobj@misc$params$species <- species
sobj@misc$params$name.int <- name.int
sobj@misc$params$names.ge <- names.GE
Seurat::Project(sobj) <- name.int
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)
save(sobj, file = paste0(norm.dim.red.dir, '/', paste(c(name.int, norm_vtr, dimred_vtr), collapse = '_'), '.rda'), compress = "bzip2")

### Correlating reduction dimensions with biases and markers expression
cat("\nCorrelation of dimensions...\n")
dimensions.eval(sobj = sobj, reduction = red.name, eval.markers = eval.markers, slot = 'data', out.dir = norm.dim.red.dir, nthreads = floor(nthreads/2))
gc()

### Testing multiple clustering parameters (nb dims kept + Louvain resolution)
cat("\nEvaluation of multiple clustering parameters...\n")
clustering.eval.mt(sobj = sobj, reduction = red.name, dimsvec = seq.int(dims.min, dims.max, dims.steps), resvec = seq(res.min,res.max,res.steps), out.dir = norm.dim.red.dir, solo.pt.size = solo.pt.size, BPPARAM = cl)
