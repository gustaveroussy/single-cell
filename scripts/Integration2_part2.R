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





load("/WORKDIR/P31_LOPO/GROUPED_ANALYSIS/INTEGRATED/LOPO_01_INT_seurat/SCTransform/pca/LOPO_01_INT_seurat_SCTransform_pca.rda")
assay <- "integrated"
red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')

### Building reduced normalized output dir
norm_vtr = paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr)) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL)
#norm_vtr = norm.method
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[red.name]]@misc$vtr)) paste(sobj@reductions[[red.name]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
#dimred_vtr = dimred.method
norm.dim.red.dir = paste0(data.path, norm_vtr, '/', dimred_vtr)


### Manually Choose optimal dims + resolution from clustree / uMAPS
keep.dims <- 19
keep.res <- 0.4

### Building clustered output directory
clust.dir <- paste(norm.dim.red.dir, paste(c(dimred.method, keep.dims, keep.res), collapse = '_'), sep = '/')
dir.create(path = clust.dir, recursive = TRUE, showWarnings = TRUE)

### Replotting final clusters
sobj <- louvain.cluster(sobj = sobj, reduction = red.name, max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size, algorithm = 1)

## Setting ident name and RNA.reduction
ident.name <- paste0(paste0(red.name, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))
INT.reduction <- paste(c(red.name, keep.dims, 'umap'), collapse = '_')
sobj@reduction[[red.name]]@misc$from.assay <- assay

### uMAP plot with harmony-regressed variable(s)
blockpix = 600
png(filename = paste0(clust.dir, '/', paste(c(sample.name.INT, red.name, 'uMAP.png'), collapse = '_')), width = 1000, height = 1000)
print(Seurat::DimPlot(object = sobj, reduction = paste(c(red.name, keep.dims, 'umap'), collapse = '_'), order = sample(x = 1:ncol(sobj), size = ncol(sobj), replace = FALSE), group.by = 'orig.ident', pt.size = solo.pt.size) + ggplot2::ggtitle(paste0("uMAP for all samples ", paste(reduction.vtr, collapse = ', '))) + Seurat::DarkTheme())
dev.off()
grid.xy <- grid.scalers(length(unique(sobj@meta.data$orig.ident)))
png(filename = paste0(clust.dir, '/', paste(c(sample.name.INT, red.name, 'split', 'uMAP.png'), collapse = '_')), width = grid.xy[1]*blockpix, height = grid.xy[2]*blockpix)
print(Seurat::DimPlot(object = sobj, reduction = paste(c(red.name, keep.dims, 'umap'), collapse = '_'), group.by = ident.name, split.by = 'orig.ident', pt.size = solo.pt.size, ncol = grid.xy[1]) + ggplot2::ggtitle(paste0("uMAP with scBFA regression on ", paste(reduction.vtr, collapse = ', '))) + Seurat::DarkTheme())
dev.off()

### Technical plots
technical.plot(sobj = sobj, ident = ident.name, out.dir = clust.dir, multi.pt.size = multi.pt.size)

### Finding markers
sobj <- find.markers.quick(sobj = sobj, ident = ident.name, test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, adjp.p.max = 5E-02, topn = 10, heatmap.cols = gradient.cols, out.dir = clust.dir)

### Automatic cell type annotation
sobj <- cells.annot(sobj = sobj, ident = ident.name, singler.setnames = singler.setnames, clustifyr.setnames = clustifyr.setnames, sr.minscore = .25, cfr.minscore = .35, out.dir = clust.dir, solo.pt.size = solo.pt.size)

### Assessing clusters : Plotting control genes
sobj <- ctrl.umap.plot(sobj = sobj, ctrl.genes = ctrl.genes, ident = ident.name, out.dir = clust.dir, dimplot.cols = gradient.cols, multi.pt.size = 2)

### Assessing clusters : Plotting provided marker genes
if(!is.null(markers)) sobj <- markers.umap.plot(sobj = sobj, markers = markers, ident = ident.name, out.dir = clust.dir, dimplot.cols = gradient.cols, multi.pt.size = 2)

### Saving final object
GE_file=paste0(clust.dir, '/', paste(c(sample.name.INT, norm_vtr, dimred_vtr, keep.dims, keep.res), collapse = "_"))
save(sobj, file = paste0(GE_file, '.rda'), compress = "bzip2")

### Building cerebro binary (with default 'orig.ident' as sample.colname)
seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, file = GE_file, nthreads = nthreads, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')
seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, file = GE_file, nthreads = nthreads, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')


