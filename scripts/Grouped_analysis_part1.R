## GROUPED PROTOCOL (no integration)

projectdir <- "/WORKDIR/B21001_CAAL_01/"
resdir <- "/WORKDIR/ressources/"
scriptsdir="/WORKDIR/scripts/"

if(TRUE){
### B21001_CAAL_01
species <- "mus_musculus" # Only 'homo_sapiens', 'mus_musculus' supported yet.
min.cells <- 0
#assay <- 'SCT'
assay <- 'RNA'
norm.method <- 'SCTransform' ## 'LogNormalize' or 'SCTransform'
features.n <- 3000
normalization.vtr <- c('percent_mt', 'percent_rb', 'nFeature_RNA', 'percent_st') #c('percent_mt', 'percent_rb', 'nFeature_RNA', 'percent_st')
reduction.vtr <- NULL #c('percent_mt', 'percent_rb', 'nFeature_RNA', 'percent_st')
dimred.method <- 'pca' ## 'scbfa' or 'bpca' or 'pca' or 'ica' or 'mds'
vtr.scale <- FALSE
max.dims <- 50
resvec <- seq(.1,1.2,.1)
multi.pt.size <- solo.pt.size <- 2
nthreads <- 4
gradient.cols <- c("gold", "blue")
gmt.file <- paste0(resdir, '/DATABASE/MSIGDB/7.1/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt')
remove.mt.genes <- FALSE
remove.rb.genes <- FALSE
remove.st.genes <- FALSE
#Keep_Norm <- TRUE #bool (if keep individual normalization (necessary for seurat integration, but for scbfa we can choose))
Keep_Norm <- FALSE
project.name <- "B21001_CAAL_01"
markfile <- paste0(projectdir, 'com/input/20211222_gene_intestine_formatted.xlsx')
ctrl.genes <- c('GAPDH')
eval.markers <- NULL
rda.list <- c(
  paste0(projectdir,"data_output/1_GE/F200_C1500_M0-0.2_R0-1/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st/pca/pca_33_1.2/1_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st_pca_33_1.2.rda"),
  paste0(projectdir,"data_output/3_GE/F200_C1500_M0-0.2_R0-1/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st/pca/pca_33_1.2/3_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st_pca_33_1.2.rda"),
  paste0(projectdir,"data_output/4_GE/F200_C1500_M0-0.2_R0-1/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st/pca/pca_25_1.2/4_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st_pca_25_1.2.rda"),
  paste0(projectdir,"data_output/921_GE/F200_C1500_M0-0.2_R0-1/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st/pca/pca_27_1.1/921_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st_pca_27_1.1.rda"),
  paste0(projectdir,"data_output/922_GE/F200_C1500_M0-0.2_R0-1/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st/pca/pca_31_1/922_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st_pca_31_1.rda"),
  paste0(projectdir,"data_output/923_GE/F200_C1500_M0-0.2_R0-1/DOUBLETSFILTER_both/SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st/pca/pca_29_1.2/923_GE_SCTransformnFeature_RNA_percent_mt_percent_rb_percent_st_pca_29_1.2.rda")
)
sample.name <- c("1","3","4","921","922","923")
sample.name.GE <- paste0(sample.name, "_GE")
sample.name.INT <- 'CAAL_01_GlobalNorm'
#sample.name.INT <- 'CAAL_01_IndivNorm'
}

if(FALSE){
# VARS : P31_LOPO
### Vars
species <- 'homo_sapiens'
min.cells <- 0
#assay <- 'SCT'
assay <- 'RNA'
norm.method <- 'SCTransform' ## 'LogNormalize' or 'SCTransform'
features.n <- 3000
normalization.vtr <- c('Cyclone.Phase', 'percent_rb', 'percent_st') ## 'Cyclone.SmG2M.Score'
reduction.vtr <- NULL #c('percent_mt', 'percent_rb', 'nFeature_RNA', 'percent_st')
dimred.method <- 'pca' ## 'scbfa' or 'bpca' or 'pca' or 'ica' or 'mds'
vtr.scale <- FALSE
max.dims <- 50
resvec <- seq(.1,1.2,.1)
multi.pt.size <- solo.pt.size <- 2
nthreads <- 4
gradient.cols <- c("gold", "blue")
gmt.file <- paste0(resdir, '/DATABASE/MSIGDB/7.1/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt')
remove.mt.genes <- FALSE
remove.rb.genes <- FALSE
remove.st.genes <- FALSE
#Keep_Norm <- TRUE #bool (if keep individual normalization (necessary for seurat integration, but for scbfa we can choose))
Keep_Norm <- FALSE
# run.name <- "IGR/200811_A00461_0118_AHTN5CDMXX"
project.name <- "B20036_YOLO_01"
markfile <- c(paste0(projectdir, 'com/input/signature_AR_NEPC_from_Beltran_NatMed2016.xlsx'),paste0(projectdir, 'com/input/signature_AR_NEPC_from_porteur_project.xlsx'),paste0(projectdir, 'com/input/signature_AR_NEPC_from_Karthaus_Science2020.xlsx'),paste0(projectdir, 'com/input/20211221_Liste_genes_pour_scRNAseq.xlsx'))
ctrl.genes <- c('GAPDH')
eval.markers <- NULL
rda.list <- c(
  paste0(projectdir,"data_output/MR009_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_mt/pca/pca_39_1.2/MR009_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_mt_pca_39_1.2.rda"),
  paste0(projectdir,"data_output/MR009_RI_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_rb/pca/pca_27_1/MR009_RI_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_rb_pca_27_1.rda"),
  paste0(projectdir,"data_output/MR041_P2_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_mt_percent_rb/pca/pca_17_0.9/MR041_P2_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_mt_percent_rb_pca_17_0.9.rda"),
  paste0(projectdir,"data_output/MR041R_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_mt_percent_rb/pca/pca_17_0.9/MR041R_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_mt_percent_rb_pca_17_0.9.rda"),
  paste0(projectdir,"data_output/MR059_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA/pca/pca_19_1.2/MR059_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_pca_19_1.2.rda"),
  paste0(projectdir,"data_output/MR077_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_mt_percent_rb/pca/pca_21_0.8/MR077_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_mt_percent_rb_pca_21_0.8.rda"),
  paste0(projectdir,"data_output/MR084_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_rb/pca/pca_29_1.2/MR084_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_rb_pca_29_1.2.rda"),
  paste0(projectdir,"data_output/MR123_P3_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_mt/pca/pca_21_0.8/MR123_P3_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_mt_pca_21_0.8.rda"),
  paste0(projectdir,"data_output/MR150_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_rb/pca/pca_21_0.5/MR150_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_rb_pca_21_0.5.rda"),
  paste0(projectdir,"data_output/MR151R_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_rb/pca/pca_35_0.9/MR151R_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_rb_pca_35_0.9.rda"),
  paste0(projectdir,"data_output/MR170_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA/pca/pca_23_0.6/MR170_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_pca_23_0.6.rda"),
  paste0(projectdir,"data_output/MR178R_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_percent_mt/pca/pca_33_0.5/MR178R_P4_graft_GE_SCTransformCyclone.Phase_percent_mt_pca_33_0.5.rda"),
  paste0(projectdir,"data_output/MR182_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase_nFeature_RNA_percent_mt_percent_rb/pca/pca_31_0.9/MR182_P4_graft_GE_SCTransformCyclone.Phase_nFeature_RNA_percent_mt_percent_rb_pca_31_0.9.rda"),
  paste0(projectdir,"data_output/MR191_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase/pca/pca_23_0.7/MR191_P4_graft_GE_SCTransformCyclone.Phase_pca_23_0.7.rda"),
  paste0(projectdir,"data_output/MR283_P4_graft_GE/F200_C1000_M0-0.2_R0-1_retain1000/DOUBLETSFILTER_both/SCTransformCyclone.Phase/pca/pca_27_0.9/MR283_P4_graft_GE_SCTransformCyclone.Phase_pca_27_0.9.rda")
)
sample.name <- c("MR009_P4_graft","MR009_RI_P4_graft","MR041_P2_graft","MR041R_P4_graft","MR059_P4_graft","MR077_P4_graft","MR084_P4_graft","MR123_P3_graft","MR150_P4_graft","MR151R_P4_graft","MR170_P4_graft","MR178R_P4_graft","MR182_P4_graft","MR191_P4_graft","MR283_P4_graft")
sample.name.GE <- paste0(sample.name, "_GE")
sample.name.INT <- 'YOLO_01_GlobalNorm'
#sample.name.INT <- 'YOLO_01_IndivNorm'
}

## Setting cells annotation bases
## Gene lists loading
# rootdir <- "/home/job/WORKSPACE/SINGLECELL/"
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

## Sourcing functions
# source("~/Bureau/scrnaseq_10X_3p_proto/bustools2seurat_preproc_functions.R")
source("/WORKDIR/scripts/bustools2seurat_preproc_functions.R")

raw.methods <- c('scbfa', 'bpca')
all.methods <- c(raw.methods, 'pca', 'mds', 'ica')

## Cheking parameters
if(Keep_Norm && norm.method != 'SCTransform') stop('Keeping normalization is only possible for "SCTransform"!')
if(Keep_Norm && assay !='SCT') stop('Keeping normalization is only possible for assay "SCT"!')
if(Keep_Norm && !is.null(normalization.vtr)) message('To keep normalization, the normalization.vtr parameter will not be used...')

## RUN
######
print("#########################")
print("Grouped_analysis_part1")
print(paste0("Keep_Norm: ",Keep_Norm))
print(paste0("sample.name.INT: ",sample.name.INT))
print("#########################")

data.path <- paste0(projectdir,'data_output/GROUPED_ANALYSIS/NO_INTEGRATED/',sample.name.INT ,'/')
dir.create(data.path, recursive = TRUE)

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)
set.seed(1337)

### Building the list of Seurat objects.
sobj.list <- sapply(seq_along(rda.list), function(x) {
  message(paste0("Loading '", rda.list[x], "' ..."))
  load(rda.list[x])
  ## Cleaning other assays, reductions and graphs
  Seurat::DefaultAssay(sobj) <- assay
  if(!Keep_Norm) sobj@assays = sobj@assays[assay]
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
  if(sample.name.GE[x] %in% c("1_GE","3_GE","4_GE")){
    sobj@meta.data$conditions=rep("RET", length(sobj@meta.data$orig.ident))
  }else{
    sobj@meta.data$conditions=rep("WT", length(sobj@meta.data$orig.ident))    
  }
  
  return(sobj)
})
names(sobj.list) <- vapply(sobj.list, Seurat::Project, 'a')
message(paste0("There are ", length(sobj.list), " samples."))

### Filtering low cells datasets # pour seurat surtout!
sobj.cells <- vapply(sobj.list, ncol, 1L)
if(any(sobj.cells < min.cells)) warning(paste0('Some datasets had less than ', min.cells, ' cells, thus were removed !'))
sobj.list <- sobj.list[sobj.cells >= min.cells]


### Add prefix for colnames of sample clustering and clean TCR/BCR
for (i in names(sobj.list)){
  # add prefix for colnames of sample clustering 
  to_rename=grep("_res\\.",colnames(sobj.list[[i]]@meta.data), value = TRUE)
  sobj.list[[i]]@meta.data[[paste0(i,'_',to_rename)]]=sobj.list[[i]]@meta.data[[to_rename]]
  sobj.list[[i]]@meta.data[[to_rename]]=NULL
  
  # cleaning integrated sobj for TCR and BCR part
  TCR_BCR_col=grep("^TCR|^BCR", colnames(sobj.list[[i]]@meta.data), value = TRUE)
  if(length(TCR_BCR_col) > 0) sobj.list[[i]]@meta.data[TCR_BCR_col] <- NULL
}

### If Keeping normlization by SCT ## Very long! 1h per sample!
if(Keep_Norm && norm.method == 'SCTransform' && assay == 'SCT'){
  # Get scale.data for each sample
  sobj.list <- sapply(seq_along(sobj.list), function(x) {
    ## Scaling if necessary
    if (sum(dim(sobj.list[[x]]@assays[[assay]]@scale.data)) < 3) {
      message(paste0("Scaling ",Seurat::Project(sobj.list[[x]])))
      #Check vtr
      scale.vtr.all <- NULL
      if(!any(is.na(sobj.list[[x]]@assays[[assay]]@misc$scaling$vtr))) {
        scale.vtr.all <- c(sobj.list[[x]]@assays[[assay]]@misc$scaling$vtr)
        message(paste0("Found scaling coveriate(s) '", paste(scale.vtr.all, collapse = "', '"), "' to regress from normalization ..."))
      }
      #Scaling
      sobj <- Seurat::ScaleData(object = sobj.list[[x]], features = rownames(sobj.list[[x]]@assays[[assay]]@counts),
                                vars.to.regress = scale.vtr.all, do.scale = FALSE, scale.max = Inf, block.size = 750)
    }
    return(sobj)
  })
}

### Merge
sobj <- merge(x = sobj.list[[1]], y = sobj.list[-1], add.cell.ids = names(sobj.list), project = sample.name.INT, merge.data = TRUE)
Seurat::Project(sobj) <- sample.name.INT
sobj@misc$params$seed <- sobj.list[[1]]@misc$params$seed

### Clean
rm(sobj.list)
gc()

### Complete the normalization
if(Keep_Norm && norm.method == 'SCTransform' && assay =='SCT'){
  ## Get HVG (because merge delete HVG)
  Seurat::VariableFeatures(sobj[[assay]]) <- rownames(sobj[[assay]]@scale.data)
  sobj@assays[[assay]]@misc$params$normalization <- list(normalization.method = 'SCTransform', assay.ori = 'SCT', assay.out = 'SCT', features.used = NA)
  sobj@assays[[assay]]@misc$scaling$vtr <- NA
  ## Saving normalized object
  save(sobj, file = paste0(data.path, sample.name.INT, '_SCTKept.rda'), compress = "bzip2")
}else{
  ### Saving normalized object
  save(sobj, file = paste0(data.path, sample.name.INT, '_NON-NORMALIZED.rda'), compress = "bzip2")
  ## Normalisation
  sobj <- sc.normalization(sobj = sobj, assay = assay, normalization.method = norm.method, features.n = features.n, vtr = normalization.vtr)
  if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
}

### Reduction dimension
sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = max.dims, vtr = reduction.vtr, vtr.scale = vtr.scale)

### Building reduced normalized output dir
norm_vtr = paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr)) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL)
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr)) paste(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
norm.dim.red.dir = paste0(data.path, norm_vtr, '/', dimred_vtr)
dir.create(path = norm.dim.red.dir, recursive = TRUE, showWarnings = TRUE)

### Saving reduced normalized object
sobj@misc$params$analysis_type=paste0("Grouped analysis; Keep individual normalization: ", Keep_Norm)
sobj@assays$RNA@misc$params$Rsession <- utils::capture.output(devtools::session_info())
save(sobj, file = paste0(norm.dim.red.dir, '/', paste(c(sample.name.INT, norm_vtr, dimred_vtr), collapse = '_'), '.rda'), compress = "bzip2")

### Correlating reduction dimensions with biases and markers expression
dimensions.eval(sobj = sobj, reduction = paste0(assay, "_", dimred.method), eval.markers = eval.markers, slot = 'data', out.dir = norm.dim.red.dir, nthreads = floor(nthreads/2))
gc()

### Testing multiple clustering parameters (nb dims kept + Louvain resolution)
clustering.eval.mt(sobj = sobj, reduction = paste0(assay, "_", dimred.method), dimsvec = seq.int(3, max.dims, 2), resvec = seq(.1,1.2,.1), out.dir = norm.dim.red.dir, solo.pt.size = solo.pt.size, BPPARAM = cl)
