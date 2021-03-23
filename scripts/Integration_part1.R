## INTEGRATION PROTOCOL (Seurat or scbfa)
## To Do: test Harmony and LIGER integration
#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.dir.rda", help="List of .rda file with individual analysis"),
  make_option("--output.dir.int", help="Output path"),
  make_option("--sample.name.int", help="Name of the integration (advice: include integration method name)"),
  make_option("--eval.markers", help="Genes to evaluate to check normalization and dimension reduction"),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  # Load data
  make_option("--min.cells", help="Number minimum of cells to keep a dataset"),
  make_option("--remove.mt.genes", help="Remove mitochondrial genes of the analysis (TRUE/FALSE)"),
  make_option("--remove.rb.genes", help="Remove ribosomal genes of the analysis (TRUE/FALSE)"),
  make_option("--remove.st.genes", help="Remove mecanical stress genes of the analysis (TRUE/FALSE)"),
  # Integration
  make_option("--integration.method", help="Name of integration method (Seurat, scbfa, Harmony or Liger)"),
  make_option("--batch.vtr", help="List of batch effect to regress into Harmony, Liger, scbfa or bpca correction ('orig.ident')"),
  # Normalization and dimension reduction
  make_option("--features.n", help="Number of High Variable Genes to consider"),
  make_option("--norm.method", help="Name of normalization method (LogNormalize or SCTransform)"),
  make_option("--dimred.method", help="Name of dimension reduction method (scbfa or bpca or pca or ica or mds)"),
  make_option("--vtr", help="List of biases to regress into normalisation or dimension reduction (percent_mt, percent_rb, nFeature_RNA, percent_st, Cyclone.Phase, and all other column name in metadata)"),
  make_option("--vtr.scale", help="TRUE to center biaises to regress (for scbfa and bpca only)"),
  make_option("--dims.max", help="Number max of dimensions to compute (depends on sample complexity and number of cells)"),
  make_option("--dims.min", help="Number min of dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--dims.steps", help="Steps for dimensions to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.max", help="Number max of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.min", help="Number min of resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  make_option("--res.steps", help="Steps for resolution to compute for evaluation (depends on sample complexity and number of cells)"),
  
  # Yaml parameters file to remplace all parameters before (to use R script without snakemake)
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
rda.list <- unlist(stringr::str_split(args$options$input.dir.rda, ","))
output.dir.int <- args$options$output.dir.int
sample.name.INT <- args$options$sample.name.int
eval.markers <- unlist(stringr::str_split(args$options$eval.markers, ","))
author.name <- args$options$author.name
author.mail <- args$options$author.mail
### Computational Parameters
nthreads <- if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Load data
min.cells <- if (!is.null(args$options$min.cells)) as.numeric(args$options$min.cells)
remove.mt.genes <- args$options$remove.mt.genes
remove.rb.genes <- args$options$remove.rb.genes
remove.st.genes <- args$options$remove.st.genes
# Integration
integration.method <- args$options$integration.method
batch.vtr <- unlist(stringr::str_split(args$options$batch.vtr, ","))
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

#### Get Missing Paramaters ####
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
# Load
if (is.null(remove.mt.genes)) remove.mt.genes <- FALSE
if (is.null(remove.rb.genes)) remove.rb.genes <- FALSE
if (is.null(remove.st.genes)) remove.st.genes <- FALSE
# Normalization and dimension reduction
if (is.null(features.n)) features.n <- 3000
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
multi.pt.size <- solo.pt.size <- 2
gradient.cols <- c("gold", "blue")
raw.methods <- c('scbfa', 'bpca')
all.methods <- c(raw.methods, 'pca', 'mds', 'ica')

#### Check non-optional parameters ####
if (is.null(rda.list)) stop("input.dir.rda parameter can't be empty!")
if (is.null(output.dir.int)) stop("output.dir.int parameter can't be empty!")
if (is.null(integration.method)) stop("integration.method parameter can't be empty!")


## Cheking parameters
if (((integration.method == "Seurat" )||(integration.method == "Liger")) && !is.null(norm.method)){
  message(paste0("With ",integration.method," integration, normaliztion of each sample is kept. Set norm.method to NULL."))
  norm.method <- NULL
}
if ((integration.method == "Liger") && !is.null(dimred.method)){
  message(paste0("With ",integration.method," integration, dimensions reduction is included. Set dimred.method to NULL."))
  dimred.method <- NULL
}
if (is.null(norm.method) && ((integration.method == "Harmony") || (integration.method == "scbfa"))) norm.method <- 'SCTransform'
if (is.null(dimred.method) && ((integration.method == "Seurat") || (integration.method == "Harmony"))) dimred.method <- 'pca'
if(is.null(dimred.method) &&(integration.method %in% raw.methods)) dimred.method <- integration.method
if((integration.method != dimred.method) && (integration.method %in% raw.methods)){
  message(paste("For ", integration.method, " integration, the same method must be used for dimension reduction! Set dimred.method to ", integration.method, "."))
  dimred.method <- integration.method
}
if (is.null(batch.vtr) && (integration.method == ("Harmony" || "Liger" || "scbfa" || "bpca"))) stop(paste0("batch.vtr parameter can't be empty for ",integration.method,"! (example: orig.ident)"))
if(integration.method == ("scbfa" || "bpca")) vtr <- c(vtr, batch.vtr)

## Sourcing functions
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

## RUN
######

print("#########################")
print(paste0("integration.method: ",integration.method))
print(paste0("sample.name.INT: ",sample.name.INT))
if (is.null(norm.method)) print(paste0("norm.method: ",norm.method))
if (is.null(dimred.method)) print(paste0("dimred.method: ",dimred.method))
print("#########################")

data.path <- paste0(output.dir.int,'/GROUPED_ANALYSIS/INTEGRATED/',sample.name.INT ,'/')
dir.create(data.path, recursive = TRUE, showWarnings = FALSE)

## Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

## Building the list of Seurat objects.
sobj.list <- sapply(seq_along(rda.list), function(x) {
  message(paste0("Loading '", rda.list[x], "' ..."))
  load(rda.list[x])
  ## Cleaning reductions and graphs
  sobj@reductions <- list()
  sobj@graphs <- list()
  ## Filtering MT genes when requested
  if (remove.mt.genes) {
    message("     Removing mt genes from the matrix ...")
    mt.genes <- readRDS(file = mt.genes.file)
    sobj <- sobj[!rownames(sobj) %in% mt.genes,]
  }
  ## Filtering RB genes when requested
  if (remove.rb.genes) {
    message("     Removing rb genes from the matrix ...")
    rb.genes <- readRDS(file = crb.genes.file)
    sobj <- sobj[!rownames(sobj) %in% rb.genes,]
  }
  ## Filtering ST genes when requested
  if (remove.st.genes) {
    message("     Removing st genes from the matrix ...")
    st.genes <- readRDS(file = str.genes.file)
    sobj <- sobj[!rownames(sobj) %in% st.genes,]
  }
  return(sobj)
})
names(sobj.list) <- vapply(sobj.list, Seurat::Project, 'a')
message(paste0("There are ", length(sobj.list), " samples."))
species <- sapply(seq_along(sobj.list), function(x) { sobj.list[[x]]@misc$params$species })
if(length(unique(species)) != 1) stop(paste0("We can't mix several species: ", paste0(species, collapse = ", ")))
species = species[1]


## Get assay parameter
if(integration.method == "scbfa"){ # redo normalisation
  assay <- 'RNA'
}
if(integration.method == "Harmony") { # redo normalisation
  assay <- 'RNA'
}
if(integration.method == "Seurat") { #keep normalisation
  n.meth <- sapply(seq_along(sobj.list), function(x) { sobj.list[[x]]@misc$params$normalization$normalization.method })
  if(length(unique(n.meth)) != 1) stop(paste0("We can't mix several normalisation method with ",integration.method,": ",paste0(n.meth, collapse = ", ")))
  if(n.meth[[1]] == "SCTransform") stop(paste0("Normalisation must be SCTransform with ",integration.method," integration: ", n.meth[[1]]))
  assay <- 'SCT'
}
if(integration.method == "Liger") { #keep normalisation
  n.meth <- sapply(seq_along(sobj.list), function(x) { sobj.list[[x]]@misc$params$normalization$normalization.method })
  if(length(unique(n.meth)) != 1) stop(paste0("We can't mix several normalisation method with ",integration.method,": ",paste0(n.meth, collapse = ", ")))
  assay <- sobj[1]@misc$params$normalization$assay.out
}


## Filtering low cells datasets # pour seurat surtout!
sobj.cells <- vapply(sobj.list, ncol, 1L)
if(any(sobj.cells < min.cells)) warning(paste0('Some datasets had less than ', min.cells, ' cells, thus were removed !'))
sobj.list <- sobj.list[sobj.cells >= min.cells]


### Add prefix for colnames of sample clustering and clean TCR/BCR
for (i in names(sobj.list)){
  # add prefix for colnames of sample clustering 
  to_rename=grep("_res\\.",colnames(sobj.list[[i]]@meta.data), value = TRUE)
  for (j in to_rename){
    sobj.list[[i]]@meta.data[[paste0(i,'_',j)]]=sobj.list[[i]]@meta.data[[j]]
    sobj.list[[i]]@meta.data[[j]]=NULL
  }
  # cleaning integrated sobj for TCR and BCR part
  TCR_BCR_col=grep("^TCR|^BCR", colnames(sobj.list[[i]]@meta.data), value = TRUE)
  if(length(TCR_BCR_col) > 0) sobj.list[[i]]@meta.data[TCR_BCR_col] <- NULL
}

## Seurat Integration
if((integration.method == "Seurat") || (integration.method == 'Liger')){
  ## Get scale.data for each sample
  message("Get scale.data for each sample:")
  sobj.list <- sapply(seq_along(sobj.list), function(x) {
    ## Scaling if necessary
    if (sum(dim(sobj.list[[x]]@assays[[assay]]@scale.data)) < 3) {
      #Check vtr
      scale.vtr.all <- NULL
      if(!any(is.na(sobj.list[[x]]@assays[[assay]]@misc$scaling$vtr))) {
        scale.vtr.all <- c(sobj.list[[x]]@assays[[assay]]@misc$scaling$vtr)
        message(paste0("Found scaling coveriate(s) '", paste(scale.vtr.all, collapse = "', '"), "' to regress from normalization ..."))
      }
      #Scaling
      Seurat::DefaultAssay(sobj.list[[x]]) <- assay
      if(assay == 'SCT') {
        sobj <- Seurat::ScaleData(object = sobj.list[[x]], features = rownames(sobj.list[[x]]@assays[[assay]]@counts),
                                  vars.to.regress = scale.vtr.all,
                                  do.scale = FALSE, scale.max = Inf, block.size = 750)
      } else {
        sobj <- Seurat::ScaleData(object = sobj.list[[x]], features = rownames(sobj.list[[x]]@assays[[assay]]@counts),
                                  vars.to.regress = scale.vtr.all,
                                  do.scale = if(integration.method == 'Liger') FALSE else TRUE, scale.max = 10, block.size = 1000)
      }
    }
    return(sobj)
  })
  ## Integration
  if(integration.method == "Seurat"){
    message(paste0(integration.method," integration..."))
    if (tolower(norm.method) == 'sctransform') int.norm.method <- 'SCT' else  int.norm.method <- 'LogNormalize'
    sobj.features <- Seurat::SelectIntegrationFeatures(object.list = sobj.list, nfeatures = 3000) #Sélection des marqueurs biologiques partagés
    sobj.list <- Seurat::PrepSCTIntegration(object.list = sobj.list, anchor.features = sobj.features) #Verifie que les résidus de Pearson ont bien été calculés
    sobj.anchors <- Seurat::FindIntegrationAnchors(object.list = sobj.list, normalization.method = int.norm.method, anchor.features = sobj.features) #CCA + L2normalisation; puis KNN; puis MNNs : identification des paires de cellules; filtrage des anchors; calcul des scores
    sobj <- Seurat::IntegrateData(anchorset = sobj.anchors, normalization.method = int.norm.method) #Calcul des poids; application des poids sur la matrice d'expression: intégration
    # Params
    Seurat::Project(sobj) <- sample.name.INT
    sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed
    DefaultAssay(sobj) <- "integrated"
    assay <- "integrated"
    sobj@assays[["integrated"]]@misc$scaling$vtr = NA
    # Cleaning
    rm(sobj.list, sobj.features, sobj.anchors)
    sobj@assays[[assay]]@scale.data <- matrix(nrow = 0, ncol = 0)
    gc()
    ## Dimensions reduction
    red.name <- paste(c("integrated", dimred.method, integration.method), collapse = '_')
    sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = "integrated", max.dims = max.dims, vtr = reduction.vtr, vtr.scale = vtr.scale, red.name = red.name)
    
  } else if(integration.method == 'Liger'){
    ## Merge data
    sobj <- merge(x = sobj.list[[1]], y = sobj.list[-1], add.cell.ids = names(sobj.list), project = sample.name.INT, merge.data = TRUE)
    Seurat::Project(sobj) <- sample.name.INT
    sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed
    print("scaling vtr:")
    print(sobj@assays[[assay]]@misc$scaling$vtr)
    #sobj@assays[[assay]]@misc$scaling$vtr = NA
    ## Cleaning
    rm(sobj.list)
    gc()
    ## Integration LIGER
    message(paste0(integration.method," integration..."))
    red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
    sobj <- SeuratWrappers::RunOptimizeALS(sobj, k = 20, lambda = 5, split.by = batch.vtr, rand.seed = sobj@misc$params$seed) #calcul les matrices (fait la red dim)
    sobj <- SeuratWrappers::RunQuantileNorm(sobj, resolution = 1, split.by = batch.vtr, reduction.name = red.name, rand.seed = sobj@misc$params$seed) #fait le SFN graph + clustering + correction de la red dim
    sobj@reductions[[red.name]]@misc$vtr = NA
  }
}

## scbfa/bpca (or Harmony integration beging)
if((integration.method %in% raw.methods) || (integration.method == 'Harmony')){
  ## Merge data
  sobj <- merge(x = sobj.list[[1]], y = sobj.list[-1], add.cell.ids = names(sobj.list), project = sample.name.INT, merge.data = TRUE)
  Seurat::Project(sobj) <- sample.name.INT
  sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed
  ## Cleaning
  rm(sobj.list)
  gc()
  ### Saving non-normalized object
  save(sobj, file = paste0(data.path, sample.name.INT, '_NON-NORMALIZED.rda'), compress = "bzip2")
  #Normalization
  sobj <- sc.normalization(sobj = sobj, assay = assay, normalization.method = norm.method, features.n = features.n, vtr = normalization.vtr)
  if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
  ## Dimensions reduction
  if(integration.method %in% raw.methods) message(paste0(integration.method," integration..."))
  red.name <- if(integration.method %in% raw.methods) paste(c(assay, dimred.method, integration.method), collapse = '_') else paste(c(assay, dimred.method), collapse = '_')
  sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = max.dims, vtr = reduction.vtr, vtr.scale = vtr.scale, red.name = red.name)
  
}

### Building reduced normalized output dir
norm_vtr = paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr)) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL, collapse = '_')
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[red.name]]@misc$vtr)) paste(sobj@reductions[[red.name]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
norm.dim.red.dir = paste0(data.path, norm_vtr, '/', dimred_vtr)
dir.create(path = norm.dim.red.dir, recursive = TRUE, showWarnings = TRUE)

## Harmony integration ending
if (integration.method == 'Harmony'){
  message(paste0(integration.method," integration..."))
  ## Scaling if necessary #NB: harmony needs scale.data
  if (sum(dim(sobj@assays[[assay]]@scale.data)) < 3) {
    #Check vtr
    scale.vtr.all <- NULL
    if(!any(is.na(sobj@assays[[assay]]@misc$scaling$vtr))) {
      scale.vtr.all <- c(sobj@assays[[assay]]@misc$scaling$vtr)
      message(paste0("Found scaling coveriate(s) '", paste(scale.vtr.all, collapse = "', '"), "' to regress from normalization ..."))
    }
    #Scaling
    Seurat::DefaultAssay(sobj) <- assay
    if(assay == 'SCT') {
      sobj <- Seurat::ScaleData(object = sobj, features = rownames(sobj@assays[[assay]]@counts),
                                vars.to.regress = scale.vtr.all,
                                do.scale = FALSE, scale.max = Inf, block.size = 750)
    } else {
      sobj <- Seurat::ScaleData(object = sobj, features = rownames(sobj@assays[[assay]]@counts),
                                vars.to.regress = scale.vtr.all,
                                do.scale = TRUE, scale.max = 10, block.size = 1000)
    }
  }
  ## Integration
  png(paste0(norm.dim.red.dir, '/harmony_convergence_plot.png'), width = 1000, height = 1000)
  sobj <- harmony::RunHarmony(sobj, batch.vtr, reduction = red.name, assay.use = assay, plot_convergence = TRUE, reduction.save = paste(c(assay, dimred.method, integration.method), collapse = '_')) #, do_pca=FALSE ??
  dev.off()
  red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
  
}

### Saving reduced normalized object
sobj@misc$params$analysis_type <- paste0("Integrated analysis; Method: ", integration.method)
sobj@misc$params$sobj_creation$Rsession <- utils::capture.output(devtools::session_info())
sobj@misc$params$species <- species
sobj@misc$params$sample.name.INT <- sample.name.INT
save(sobj, file = paste0(norm.dim.red.dir, '/', paste(c(sample.name.INT, norm_vtr, dimred_vtr), collapse = '_'), '.rda'), compress = "bzip2")

### Correlating reduction dimensions with biases and markers expression
dimensions.eval(sobj = sobj, reduction = red.name, eval.markers = eval.markers, slot = 'data', out.dir = norm.dim.red.dir, nthreads = floor(nthreads/2))
gc()

### Testing multiple clustering parameters (nb dims kept + Louvain resolution)
clustering.eval.mt(sobj = sobj, reduction = red.name, dimsvec = seq.int(3, max.dims, 2), resvec = seq(.1,1.2,.1), out.dir = norm.dim.red.dir, solo.pt.size = solo.pt.size, BPPARAM = cl)
