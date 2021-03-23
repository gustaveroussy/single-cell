library(optparse)
option_list <- list(
  make_option(c("-s", "--samples"), default=NULL, help="Name of sample"),
  make_option(c("-r", "--rda"), default=NULL, help="path/file_name.rda"),
  make_option(c("-d", "--dim"), default=NULL, help="dimensions to keep"),
  make_option(c("-g", "--res"), default=NULL, help="resolution parameter to used") #granularity
)

parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)

# SAMPLE/DIR
sample.name.INT <- args$options$samples
rda_file <- args$options$rda
keep.dims <- as.numeric(args$options$dim)
keep.res <- as.numeric(args$options$res)

print("#####################################")
print(paste0("Sample: ", sample.name.INT))
print(paste0("RDA file: ", rda_file))
print(paste0("Dimension: ", keep.dims))
print(paste0("Resolution: ", keep.res))
print("#####################################")

if(FALSE){
### B21001_CAAL_01
projectdir <- "/WORKDIR/B21001_CAAL_01/"
resdir <- "/WORKDIR/ressources/"
scriptsdir="/WORKDIR/scripts/"
species <- "mus_musculus" # Only 'homo_sapiens', 'mus_musculus' supported yet.
multi.pt.size <- solo.pt.size <- 2
nthreads <- 4
gradient.cols <- c("gold", "blue")
gmt.file <- paste0(resdir, '/DATABASE/MSIGDB/7.1/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt')
remove.mt.genes <- FALSE
remove.rb.genes <- FALSE
remove.st.genes <- FALSE
project.name <- "B21001_CAAL_01"
markfile <- paste0(projectdir, 'com/input/20211222_gene_intestine_formatted.xlsx')
ctrl.genes <- c('GAPDH')
eval.markers <- NULL
}

if(TRUE){
# VARS : P31_LOPO
projectdir <- "/WORKDIR/B20036_YOLO_01/"
resdir <- "/WORKDIR/ressources/"
scriptsdir="/WORKDIR/scripts/"
species <- 'homo_sapiens'
multi.pt.size <- solo.pt.size <- 2
nthreads <- 8
gradient.cols <- c("gold", "blue")
gmt.file <- paste0(resdir, '/DATABASE/MSIGDB/7.1/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt')
remove.mt.genes <- FALSE
remove.rb.genes <- FALSE
remove.st.genes <- FALSE
project.name <- "B20036_YOLO_01"
markfile <- c(paste0(projectdir, 'com/input/signature_AR_NEPC_from_Beltran_NatMed2016.xlsx'),paste0(projectdir, 'com/input/signature_AR_NEPC_from_porteur_project.xlsx'),paste0(projectdir, 'com/input/signature_AR_NEPC_from_Karthaus_Science2020.xlsx'),paste0(projectdir, 'com/input/20211221_Liste_genes_pour_scRNAseq.xlsx'))
ctrl.genes <- c('GAPDH')
eval.markers <- NULL
}

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
#get genes markers
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
## Sourcing functions
source("/WORKDIR/scripts/bustools2seurat_preproc_functions.R")

raw.methods <- c('scbfa', 'bpca')
all.methods <- c(raw.methods, 'pca', 'mds', 'ica')


## RUN
######
data.path <- paste0(projectdir,'data_output/GROUPED_ANALYSIS/NO_INTEGRATED/',sample.name.INT ,'/')

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### load data
load(rda_file)

### Basic normalization and dimension reduction
if(length(grep("SCTransform", rda_file))==1) norm.method <- 'SCTransform' else norm.method <- "LogNormalize"
if(tolower(norm.method) == 'sctransform') assay <- 'SCT' else assay <- 'RNA'
if(length(grep("pca", rda_file))==1) dimred.method <- 'pca' else dimred.method <- 'scbfa'

### path
norm_vtr = paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr)) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL)
dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr)) paste(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
norm.dim.red.dir = paste0(data.path, norm_vtr, '/', dimred_vtr)



### Building clustered output directory
clust.dir <- paste(norm.dim.red.dir, paste(c(dimred.method, keep.dims, keep.res), collapse = '_'), sep = '/')


if(TRUE){
dir.create(path = clust.dir, recursive = TRUE, showWarnings = TRUE)

### Replotting final clusters
sobj <- louvain.cluster(sobj = sobj, reduction = paste0(assay, "_", dimred.method), max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size, algorithm = 1)

## Setting ident name and RNA.reduction
ident.name <- paste0(paste0(assay, "_", dimred.method, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))

### uMAP plot by sample
blockpix = 600
png(filename = paste0(clust.dir, '/', paste(c(sample.name.INT, assay, dimred.method, 'uMAP.png'), collapse = '_')), width = 1000, height = 1000)
print(Seurat::DimPlot(object = sobj, reduction = paste(c(assay, dimred.method, keep.dims, 'umap'), collapse = '_'), order = sample(x = 1:ncol(sobj), size = ncol(sobj), replace = FALSE), group.by = 'orig.ident', pt.size = solo.pt.size) + ggplot2::ggtitle("uMAP for all samples ") + Seurat::DarkTheme())
dev.off()
grid.xy <- grid.scalers(length(unique(sobj@meta.data$orig.ident)))
png(filename = paste0(clust.dir, '/', paste(c(sample.name.INT, assay, dimred.method, 'split', 'uMAP.png'), collapse = '_')), width = grid.xy[1]*blockpix, height = grid.xy[2]*blockpix)
print(Seurat::DimPlot(object = sobj, reduction = paste(c(assay, dimred.method, keep.dims, 'umap'), collapse = '_'), group.by = ident.name, split.by = 'orig.ident', pt.size = solo.pt.size, ncol = grid.xy[1]) + ggplot2::ggtitle(paste0("uMAP split on samples")) + Seurat::DarkTheme())
dev.off()

### Technical plots
technical.plot(sobj = sobj, ident = ident.name, out.dir = clust.dir, multi.pt.size = multi.pt.size)

### Finding markers
sobj <- find.markers.quick(sobj = sobj, ident = ident.name, test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, adjp.p.max = 5E-02, topn = 10, heatmap.cols = gradient.cols, out.dir = clust.dir)

### Automatic cell type annotation
sobj <- cells.annot(sobj = sobj, ident = ident.name, singler.setnames = singler.setnames, clustifyr.setnames = clustifyr.setnames, sr.minscore = .25, cfr.minscore = .35, out.dir = clust.dir, solo.pt.size = solo.pt.size)

### Assessing clusters : Plotting control genes
if(!is.null(ctrl.genes)) sobj <- ctrl.umap.plot(sobj = sobj, ctrl.genes = ctrl.genes, ident = ident.name, out.dir = clust.dir, dimplot.cols = gradient.cols, multi.pt.size = 2)

### Assessing clusters : Plotting provided marker genes
if(!is.null(markers)) sobj <- markers.umap.plot(sobj = sobj, markers = markers, ident = ident.name, out.dir = clust.dir, dimplot.cols = gradient.cols, multi.pt.size = 2)

### Saving final object
#Keep_Norm=TRUE
#sobj@misc$params$analysis_type=paste0("Grouped analysis; Keep individual normalization: ", Keep_Norm)
GE_file=paste0(clust.dir, '/', paste(c(sample.name.INT, norm_vtr, dimred_vtr, keep.dims, keep.res), collapse = "_"))
save(sobj, file = paste0(GE_file, '.rda'), compress = "bzip2")

### Building cerebro binary (with default 'orig.ident' as sample.colname)
seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, file = GE_file, nthreads = nthreads, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')
seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, file = GE_file, nthreads = nthreads, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')

}else{
  ident.name <- paste0(paste0(red.name, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))
  GE_file=paste0(clust.dir, '/', paste(c(sample.name.INT, norm_vtr, dimred_vtr, keep.dims, keep.res), collapse = "_"))
  load(paste0(GE_file, '.rda'))
  print("Building cerebro binary")
  seurat2cerebro_1.3(sobj = sobj, ident = ident.name, groups='conditions', remove.other.reductions = FALSE, remove.other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, file = GE_file, nthreads = nthreads, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 5E-02, test = 'wilcox')
  seurat2cerebro_1.3(sobj = sobj, ident = ident.name, groups='conditions', remove.other.reductions = FALSE, remove.other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, file = GE_file, nthreads = nthreads, min_pct = .75, thresh_logFC = .5, thresh_p_val = 5E-02, test = 'wilcox')

}

