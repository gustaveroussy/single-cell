## INTEGRATION PROTOCOL (Seurat or scbfa)
## To Do: add Harmony and LIGER integration

# projectdir <- "/WORKDIR/P31_LOPO/"
# resdir <- "/WORKDIR/ressources/"
# scriptsdir="/WORKDIR/scripts/"

# projectdir <- "~/Bureau/Test_int/"
# resdir <- "/home/m_aglave/Bureau/SCRNASEQ/RESOURCES/"
projectdir <- "/WORKDIR/P31_LOPO/"
resdir <- "/WORKDIR/ressources/"
scriptsdir="/WORKDIR/scripts/"

### Vars
species <- 'homo_sapiens'
min.cells <- 0
assay <- 'SCT'
norm.method <- 'SCTransform' ## 'LogNormalize' or 'SCTransform'
features.n <- 3000
normalization.vtr <- NULL
reduction.vtr <- NULL #c('Cyclone.SmG2M.Score')
dimred.method <- 'pca' ## 'scbfa' or 'bpca' or 'pca' or 'ica' or 'mds'
vtr.scale=FALSE
max.dims <- 50
resvec <- seq(.1,1.2,.1)
multi.pt.size <- solo.pt.size <- 2
nthreads <- 8
gradient.cols <- c("gold", "blue")
remove.mt.genes <- FALSE
remove.rb.genes <- FALSE
remove.st.genes <- FALSE
integration.method <- 'seurat' #seurat, scbfa, bpca (Harmony, LIGER)
# harmony.vtr <- 'orig.ident'

# VARS : FXDA
# rda.list <- c(
#   "/home/m_aglave/Bureau/SCRNASEQ/RUN/IGR/200525_A00461_0112_AHKTJKDRXX/P31_FXDA_1/855m_IPI_GE/F200_C1000_M0-0.1_R0-1/DOUBLETSFILTER_both/LogNormalize/scbfa_Cyclone.Phase/scbfa_17_0.9/855m_IPI_GE_LogNormalize_scbfaCyclone.Phase_17_0.9.rda",
#   "/home/m_aglave/Bureau/SCRNASEQ/RUN/IGR/200525_A00461_0112_AHKTJKDRXX/P31_FXDA_1/855m_ISOTYPE_GE/F200_C1000_M0-0.1_R0-1/DOUBLETSFILTER_both/LogNormalize/scbfa/scbfa_17_0.9/855m_ISOTYPE_GE_LogNormalize_scbfa_17_0.9.rda",
#   "/home/m_aglave/Bureau/SCRNASEQ/RUN/IGR/200525_A00461_0112_AHKTJKDRXX/P31_FXDA_1/855m_MEDIUM_GE/F200_C1000_M0-0.1_R0-1/DOUBLETSFILTER_both/LogNormalize/scbfa/scbfa_19_0.9/855m_MEDIUM_GE_LogNormalize_scbfa_19_0.9.rda"
# )
# sample.name <- c("855m_IPI","855m_ISOTYPE","855m_MEDIUM")
# sample.name.GE <- paste0(sample.name, "_GE")
# sample.name.INT <- 'FXDA_INT_SCT_pca_Harmony'
# eval.markers <- NULL
# ctrl.genes <- c('GAPDH')

# # VARS : P31_LOPO
run.name <- "IGR/200811_A00461_0118_AHTN5CDMXX"
project.name <- "P31_LOPO"
data.path <- paste0(projectdir,'GROUPED_ANALYSIS/INTEGRATED/')
markfile1 <- paste0(projectdir, 'communication/in/signature_AR_NEPC_from_Beltran_NatMed2016.xlsx')
mark1.xl <- openxlsx::read.xlsx(markfile1, sheet = 1, startRow = 1, fillMergedCells = TRUE, colNames = TRUE)
mark1.xl <- mark1.xl[order(mark1.xl[,2]),]
markers1 <- setNames(mark1.xl[,1], mark1.xl[,2])
markfile2 <- paste0(projectdir, 'communication/in/signature_AR_NEPC_from_porteur_project.xlsx')
mark2.xl <- openxlsx::read.xlsx(markfile2, sheet = 1, startRow = 1, fillMergedCells = TRUE, colNames = TRUE)
mark2.xl <- mark2.xl[order(mark2.xl[,2]),]
markers2 <- setNames(mark2.xl[,1], mark2.xl[,2])
markfile3 <- paste0(projectdir, 'communication/in/signature_AR_NEPC_from_Karthaus_Science 2020.xlsx')
mark3.xl <- openxlsx::read.xlsx(markfile3, sheet = 1, startRow = 1, fillMergedCells = TRUE, colNames = TRUE)
mark3.xl <- mark3.xl[order(mark3.xl[,2]),]
markers3 <- setNames(mark3.xl[,1], mark3.xl[,2])
markers <- c(markers1, markers2, markers3)
ctrl.genes <- c('GAPDH')
eval.markers <- NULL
rda.list <- c(
  paste0(projectdir,"MR009_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb/pca/pca_15_0.3/MR009_P4_graft_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_pca_15_0.3.rda"),
  paste0(projectdir,"MR041_P2_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb/pca/pca_15_0.4/MR041_P2_graft_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_pca_15_0.4.rda"),
  paste0(projectdir,"MR041R_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb/pca/pca_15_0.2/MR041R_P4_graft_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_pca_15_0.2.rda"),
  paste0(projectdir,"MR077_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt/pca/pca_15_0.5/MR077_P4_graft_GE_SCTransformnFeature_RNA_percent_mt_pca_15_0.5.rda"),
  paste0(projectdir,"MR084_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt/pca/pca_19_0.1/MR084_P4_graft_GE_SCTransformnFeature_RNA_percent_mt_pca_19_0.1.rda"),
  paste0(projectdir,"MR123_P3_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_st/pca/pca_19_0.3/MR123_P3_graft_GE_SCTransformnFeature_RNA_percent_mt_percent_st_pca_19_0.3.rda"),
  paste0(projectdir,"MR151R_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st/pca/pca_23_0.4/MR151R_P4_graft_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st_pca_23_0.4.rda"),
  paste0(projectdir,"MR170_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt/pca/pca_15_0.2/MR170_P4_graft_GE_SCTransformnFeature_RNA_percent_mt_pca_15_0.2.rda"),
  paste0(projectdir,"MR178R_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_st/pca/pca_13_0.2/MR178R_P4_graft_GE_SCTransformnFeature_RNA_percent_mt_percent_st_pca_13_0.2.rda"),
  paste0(projectdir,"MR182_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st/pca/pca_13_0.1/MR182_P4_graft_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st_pca_13_0.1.rda"),
  paste0(projectdir,"MR191_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransform/pca/pca_21_0.2/MR191_P4_graft_GE_SCTransform_pca_21_0.2.rda"),
  paste0(projectdir,"MR283_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_st/pca/pca_17_0.5/MR283_P4_graft_GE_SCTransformnFeature_RNA_percent_st_pca_17_0.5.rda")
)
sample.name <- c("MR009_P4_graft","MR041_P2_graft","MR041R_P4_graft","MR077_P4_graft","MR084_P4_graft","MR123_P3_graft","MR151R_P4_graft","MR170_P4_graft","MR178R_P4_graft","MR182_P4_graft","MR191_P4_graft","MR283_P4_graft")
sample.name.GE <- paste0(sample.name, "_GE")
sample.name.INT <- 'LOPO_01_INT_seurat'


## Setting cells annotation bases
## Gene lists loading
if (species == "homo_sapiens") {
  species.rdx <- 'hg'
  mt.genes.file <- paste0(resdir, "GENELISTS/homo_sapiens_mito_symbols_20191001.rds")
  crb.genes.file <- paste0(resdir, "GENELISTS/homo_sapiens_cribo_symbols_20191015.rds")
  cc.genes.file <- paste0(resdir, "GENELISTS/homo_sapiens_Seurat_cc.genes_20191031.rds")
  cc.pairs.file <- paste0(resdir, "GENELISTS/homo_sapiens_cyclone_pairs_symbols_20191001.rds")
  str.genes.file <- paste0(resdir, "GENELISTS/homo_sapiens_stress_symbols_20200224.rds")
  singler.setnames <- c("HumanPrimaryCellAtlasData", "NovershternHematopoieticData", "DatabaseImmuneCellExpressionData", "MonacoImmuneData")
  clustifyr.setnames <- c("pbmc_avg", "hema_microarray_matrix", "gtex_bulk_matrix")
}
if (species == "mus_musculus") {
  species.rdx <- 'mm'
  mt.genes.file <- paste0(resdir, "GENELISTS/mus_musculus_mito_symbols_20191015.rds")
  crb.genes.file <- paste0(resdir, "GENELISTS/mus_musculus_cribo_symbols_20191015.rds")
  cc.genes.file <- paste0(resdir, "GENELISTS/mus_musculus_Seurat_cc.genes_20191031.rds")
  cc.pairs.file <- paste0(resdir, "GENELISTS/mus_musculus_cyclone_pairs_symbols_20191015.rds")
  str.genes.file <- paste0(resdir, "GENELISTS/mus_musculus_stress_symbols_20200224.rds")
  singler.setnames <- c("MouseRNAseqData", "ImmGenData")
  clustifyr.setnames <- c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq")
}
if (species == "rattus_norvegicus") {
  species.rdx <- 'rn'
  mt.genes.file <- paste0(resdir, "GENELISTS/rattus_norvegicus_mito_symbols_20200315.rds")
  crb.genes.file <- paste0(resdir, "GENELISTS/rattus_norvegicus_cribo_symbols_20200315.rds")
  cc.genes.file <- paste0(resdir, "GENELISTS/rattus_norvegicus_Seurat_cc.genes_20200315.rds")
  cc.pairs.file <- paste0(resdir, "GENELISTS/rattus_norvegicus_cyclone_pairs_symbols_20200315.rds")
  str.genes.file <- paste0(resdir, "GENELISTS/rattus_norvegicus_stress_symbols_20200315.rds")
  singler.setnames <- c("MouseRNAseqData", "ImmGenData")
  clustifyr.setnames <- c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq")
}
gmt.file <- paste0(resdir, 'DATABASE/MSIGDB/7.1/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt')

## Sourcing functions
# source("~/Bureau/scrnaseq_10X_3p_proto/bustools2seurat_preproc_functions.R")
source("/WORKDIR/scripts/bustools2seurat_preproc_functions.R")

raw.methods <- c('scbfa', 'bpca')
all.methods <- c(raw.methods, 'pca', 'mds', 'ica')

## Cheking parameters
if(integration.method %in% raw.methods && integration.method != dimred.method) stop('For scbfa/bpca, the same method must be used for integration and dimension reduction!')
if(integration.method == 'Harmony' && dimred.method != 'pca') stop('For Harmony integration, the dimension reduction method must be "pca"!')

## RUN
######
data.path <- paste0(projectdir,'GROUPED_ANALYSIS/INTEGRATED/',sample.name.INT ,'/')
dir.create(data.path, recursive = TRUE)

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
if(integration.method == "seurat"){
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
                                  do.scale = TRUE, scale.max = 10, block.size = 1000)
      }
    }
    return(sobj)
  })
  ## Integration
  message("Seurat integration...")
  if (tolower(norm.method) == 'sctransform') int.norm.method <- 'SCT' else  int.norm.method <- 'LogNormalize'
  sobj.features <- SelectIntegrationFeatures(object.list = sobj.list, nfeatures = 3000) #Sélection des marqueurs biologiques partagés
  sobj.list <- PrepSCTIntegration(object.list = sobj.list, anchor.features = sobj.features) #Verifie que les résidus de Pearson ont bien été calculés
  sobj.anchors <- FindIntegrationAnchors(object.list = sobj.list, normalization.method = int.norm.method, anchor.features = sobj.features) #CCA + L2normalisation; puis KNN; puis MNNs : identification des paires de cellules; filtrage des anchors; calcul des scores
  sobj <- IntegrateData(anchorset = sobj.anchors, normalization.method = int.norm.method) #Calcul des poids; application des poids sur la matrice d'expression: intégration
  # Params
  Seurat::Project(sobj) <- sample.name.INT
  sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed
  DefaultAssay(sobj) <- "integrated"
  ori.assay <- assay
  assay <- "integrated"
  sobj@assays[[assay]]@misc$scaling$vtr = NA
  # Cleaning
  rm(sobj.list, sobj.features, sobj.anchors)
  sobj@assays[[ori.assay]]@scale.data <- matrix(nrow = 0, ncol = 0)
  gc()
  ## Dimensions reduction
  red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
  sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = max.dims, vtr = reduction.vtr, vtr.scale = vtr.scale, red.name = red.name)
}

## scbfa/bpca (or Harmony integration beging)
if((integration.method %in% raw.methods) || (integration.method == 'Harmony')){
  if(integration.method %in% raw.methods) message(paste0(integration.method," integration..."))
  ## Merge data
  sobj <- merge(x = sobj.list[[1]], y = sobj.list[-1], add.cell.ids = names(sobj.list), project = sample.name.INT, merge.data = TRUE)
  Seurat::Project(sobj) <- sample.name.INT
  sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed
  ## Cleaning
  rm(sobj.list)
  gc()
  ### Saving normalized object
  save(sobj, file = paste0(data.path, sample.name.INT, '_NON-NORMALIZED.rda'), compress = "bzip2")
  #Normalization
  sobj <- sc.normalization(sobj = sobj, assay = assay, normalization.method = norm.method, features.n = features.n, vtr = normalization.vtr)
  if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
  ## Dimensions reduction
  red.name <- if(integration.method %in% raw.methods) paste(c(assay, dimred.method, integration.method), collapse = '_') else paste(c(assay, dimred.method), collapse = '_')
  sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = max.dims, vtr = reduction.vtr, vtr.scale = vtr.scale, red.name = red.name)
}


### Building reduced normalized output dir
norm_vtr = paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr)) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL)
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[red.name]]@misc$vtr)) paste(sobj@reductions[[red.name]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
norm.dim.red.dir = paste0(data.path, norm_vtr, '/', dimred_vtr)
dir.create(path = norm.dim.red.dir, recursive = TRUE, showWarnings = TRUE)

### Saving reduced normalized object
save(sobj, file = paste0(norm.dim.red.dir, '/', paste(c(sample.name.INT, norm_vtr, dimred_vtr), collapse = '_'), '.rda'), compress = "bzip2")

### Correlating reduction dimensions with biases and markers expression
dimensions.eval(sobj = sobj, reduction = red.name, eval.markers = eval.markers, slot = 'data', out.dir = norm.dim.red.dir, nthreads = floor(nthreads/2))
gc()

## PAUSE : Observe Correlation Bias !

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
  sobj <- harmony::RunHarmony(sobj, harmony.vtr, reduction = red.name, assay.use = assay, plot_convergence = TRUE, reduction.save = paste(c(assay, dimred.method, integration.method), collapse = '_'))
  red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
}
  
### Testing multiple clustering parameters (nb dims kept + Louvain resolution)
clustering.eval.mt(sobj = sobj, reduction = red.name, dimsvec = seq.int(3, max.dims, 2), resvec = seq(.1,1.2,.1), out.dir = norm.dim.red.dir, solo.pt.size = solo.pt.size, BPPARAM = cl)

## PAUSE : Observe ALL UMAPS !
