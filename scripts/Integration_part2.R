## INTEGRATION PROTOCOL (Seurat or scbfa)
## To Do: test Harmony and LIGER integration

library(optparse)
option_list <- list(
  make_option(c("-s", "--samples"), default=NULL, help="Name of sample"),
  make_option(c("-r", "--rda"), default=NULL, help="path/file_name.rda"),
  make_option(c("-d", "--dim"), default=NULL, help="dimensions to keep"),
  make_option(c("-g", "--res"), default=NULL, help="resolution parameter to used"), #granularity
  make_option(c("-t", "--tcr"), default=NULL, help="CellRanger TCR results")
)

parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)

# SAMPLE/DIR
sample.name.INT <- args$options$samples
rda_file <- args$options$rda
keep.dims <- as.numeric(args$options$dim)
keep.res <- as.numeric(args$options$res)
vdj.input.file <- unlist(stringr::str_split(args$options$tcr, ","))

print("#####################################")
print(paste0("Sample: ", sample.name.INT))
print(paste0("RDA file: ", rda_file))
print(paste0("Dimension: ", keep.dims))
print(paste0("Resolution: ", keep.res))
if(!is.null(vdj.input.file)) print(paste0("vdj.input.file: ", vdj.input.file))
print("#####################################")



if(FALSE){
# B21001_CAAL_01
projectdir <- "/WORKDIR/B21001_CAAL_01/"
resdir <- "/WORKDIR/ressources/"
scriptsdir="/WORKDIR/scripts/"
### Vars
species <- 'mus_musculus'
min.cells <- 0 #nb de cellules min dans un dataset
assay <- 'SCT'
norm.method <- 'SCTransform' ## 'LogNormalize' or 'SCTransform'
features.n <- 3000
normalization.vtr <- NULL
reduction.vtr <- NULL #c('Cyclone.SmG2M.Score')
dimred.method <- 'pca' ## 'scbfa' or 'bpca' or 'pca' or 'ica' or 'mds'
vtr.scale <- FALSE
max.dims <- 50
multi.pt.size <- solo.pt.size <- 2
nthreads <- 4
gradient.cols <- c("gold", "blue")
remove.mt.genes <- FALSE
remove.rb.genes <- FALSE
remove.st.genes <- FALSE
integration.method <- 'Seurat' #Seurat, scbfa, bpca (Harmony, Liger)
harmony.liger.vtr <- NULL # on peut en mettre plusieurs pour Harmony, mais un seul pour Liger?
project.name <- "B21001_CAAL_01"
eval.markers <- NULL # list of gene to evaluate the dimensions reduction
ctrl.genes <- c('GAPDH', eval.markers) # list of controle genes
markfile <- paste0(projectdir, 'com/input/20211222_gene_intestine_formatted.xlsx')
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
}

## Setting cells annotation bases
## Gene lists loading
if (species == "homo_sapiens") {
  species.rdx <- 'hg'
  singler.setnames <- c("HumanPrimaryCellAtlasData", "NovershternHematopoieticData", "DatabaseImmuneCellExpressionData", "MonacoImmuneData")
  clustifyr.setnames <- c("pbmc_avg", "hema_microarray_matrix", "gtex_bulk_matrix")
}
if (species == "mus_musculus") {
  species.rdx <- 'mm'
  singler.setnames <- c("MouseRNAseqData", "ImmGenData")
  clustifyr.setnames <- c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq")
}
if (species == "rattus_norvegicus") {
  species.rdx <- 'rn'
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
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

#fixed parameters
raw.methods <- c('scbfa', 'bpca')
all.methods <- c(raw.methods, 'pca', 'mds', 'ica')

#path
data.path <- paste0(projectdir,'data_output/GROUPED_ANALYSIS/INTEGRATED/',sample.name.INT ,'/')

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### load data
load(rda_file)

if(is.null(vdj.input.file)){
  ### Basic normalization and dimension reduction
  if(length(grep("SCTransform", rda_file))==1) norm.method <- 'SCTransform'
  if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
  if(length(grep("pca", rda_file))==1) dimred.method <- 'pca'
  if(integration.method == "Seurat") assay <- "integrated"
  
  ### path
  if(integration.method == "Seurat") red.name <- paste(c("integrated", dimred.method, integration.method), collapse = '_')
  if(integration.method %in% c(raw.methods,'Liger')) red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
  if(integration.method == "Harmony") red.name <- paste(c(assay, dimred.method), collapse = '_')
  norm_vtr = paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr)) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL, collapse = '_')
  dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[red.name]]@misc$vtr)) paste(sobj@reductions[[red.name]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
  norm.dim.red.dir = paste0(data.path, norm_vtr, '/', dimred_vtr)
  if(integration.method == "Harmony")   red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')

  ### Building clustered output directory
  clust.dir <- paste(norm.dim.red.dir, paste(c(dimred.method, keep.dims, keep.res), collapse = '_'), sep = '/')
  dir.create(path = clust.dir, recursive = TRUE, showWarnings = TRUE)
  
  ### Replotting final clusters
  sobj <- louvain.cluster(sobj = sobj, reduction = red.name, max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size, algorithm = 1)
  
  ## Setting ident name and RNA.reduction
  ident.name <- paste0(paste0(red.name, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))
  INT.reduction <- paste(c(red.name, keep.dims, 'umap'), collapse = '_')
  sobj@reductions[[red.name]]@misc$from.assay <- assay #necessary for seurat integration
  
  ### uMAP plot by sample
  blockpix = 600
  png(filename = paste0(clust.dir, '/', paste(c(sample.name.INT, red.name, 'uMAP.png'), collapse = '_')), width = 1000, height = 1000)
  print(Seurat::DimPlot(object = sobj, reduction = paste(c(red.name, keep.dims, 'umap'), collapse = '_'), order = sample(x = 1:ncol(sobj), size = ncol(sobj), replace = FALSE), group.by = 'orig.ident', pt.size = solo.pt.size) + ggplot2::ggtitle("uMAP for all samples ") + Seurat::DarkTheme())
  dev.off()
  grid.xy <- grid.scalers(length(unique(sobj@meta.data$orig.ident)))
  png(filename = paste0(clust.dir, '/', paste(c(sample.name.INT, red.name, 'split', 'uMAP.png'), collapse = '_')), width = grid.xy[1]*blockpix, height = grid.xy[2]*blockpix)
  print(Seurat::DimPlot(object = sobj, reduction = paste(c(red.name, keep.dims, 'umap'), collapse = '_'), group.by = ident.name, split.by = 'orig.ident', pt.size = solo.pt.size, ncol = grid.xy[1]) + ggplot2::ggtitle(paste0("uMAP split on samples")) + Seurat::DarkTheme())
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
  sobj@misc$params$analysis_type=paste0("Integrated analysis; Method: ", integration.method)
  GE_file=paste0(clust.dir, '/', paste(c(sample.name.INT, norm_vtr, dimred_vtr, keep.dims, keep.res), collapse = "_"))
  save(sobj, file = paste0(GE_file, '.rda'), compress = "bzip2")
  
  ### Building cerebro binary (with default 'orig.ident' as sample.colname
  # print("Building cerebro binary")
  #v1.2.2
  # seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, file = GE_file, nthreads = nthreads, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')
  # seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, file = GE_file, nthreads = nthreads, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')
  #v1.3
  # seurat2cerebro_1.3(sobj = sobj, ident = ident.name, groups='conditions', remove.other.reductions = FALSE, remove.other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, file = GE_file, nthreads = nthreads, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 5E-02, test = 'wilcox')
  # seurat2cerebro_1.3(sobj = sobj, ident = ident.name, groups='conditions', remove.other.reductions = FALSE, remove.other.idents = FALSE, species = species.rdx, gmt.file = gmt.file, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, file = GE_file, nthreads = nthreads, min_pct = .75, thresh_logFC = .5, thresh_p_val = 5E-02, test = 'wilcox')
}

if(!is.null(vdj.input.file)){
  
  ### Basic normalization and dimension reduction
  if(length(grep("SCTransform", rda_file))==1) norm.method <- 'SCTransform'
  if(tolower(norm.method) == 'sctransform') assay <- 'SCT'
  if(length(grep("pca", rda_file))==1) dimred.method <- 'pca'
  
  ### path
  if(integration.method == "Seurat") red.name <- paste(c("integrated", dimred.method, integration.method), collapse = '_')
  if(integration.method %in% c(raw.methods,'Liger')) red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
  if(integration.method == "Harmony") red.name <- paste(c(assay, dimred.method), collapse = '_')
  norm_vtr = paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr)) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL, collapse = '_')
  dimred_vtr = paste0(c(dimred.method, if(!is.na(sobj@reductions[[red.name]]@misc$vtr)) paste(sobj@reductions[[red.name]]@misc$vtr, collapse = '_') else NULL), collapse = '_')
  norm.dim.red.dir = paste0(data.path, norm_vtr, '/', dimred_vtr)
  if(integration.method == "Harmony")   red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
  
  ### Building clustered output directory
  clust.dir <- paste(norm.dim.red.dir, paste(c(dimred.method, keep.dims, keep.res), collapse = '_'), sep = '/')
  
  ## Setting ident name and RNA.reduction
  ident.name <- paste0(paste0(red.name, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))
  INT.reduction <- paste(c(red.name, keep.dims, 'umap'), collapse = '_')
  
  ### Saving final object
  GE_file=paste0(clust.dir, '/', paste(c(sample.name.INT, norm_vtr, dimred_vtr, keep.dims, keep.res), collapse = "_"))
  
  load(paste0(GE_file, '.rda'))

  ### OTHERS : TCR
  require(patchwork)
  suppressMessages(require(Seurat))
  library(dplyr)
  set.seed(sobj@misc$params$seed)
  
  ## output directory
  output_path_TCR <- paste0(clust.dir, "/TCR_results/")
  
  
  if(FALSE){
    ## searching clusters results and representation
    table(sobj[[ident.name]])
    Seurat::Idents(sobj) <- sobj[[ident.name]]
    
    ## global variables
    list_type_clT <- c("gene+nt", "gene", "nt", "aa")
    list_type_contig <- c("nt","aa")
    caption='"gene"? - use the genes comprising the TCR
"nt" - use the nucleotide sequence of the CDR3 region
"aa" - use the amino acid sequence of the CDR3 region
"gene+nt" - use the genes comprising the TCR + the nucleotide sequence of the CDR3 region for T cells. This is the proper definition of clonotype.'
    
    
    
    
    ## GLOBAL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print("Global_analysis")
    global_output <- paste0(output_path_TCR, "Global_analysis")
    dir.create(path = global_output, recursive = TRUE, showWarnings = TRUE)
    
    ## Loading input data and Combining contigs
    cr_res <- lapply(seq_along(vdj.input.file), load.sc.tcr.bcr, sobj=sobj, vdj.input.file=vdj.input.file, sample.name=sample.name.GE)
    tcr.combined <- scRepertoire::combineTCR(df = cr_res, samples = sample.name, ID = rep("TCR", length(sample.name)), cells = "T-AB")
    
    ## Quantification of unique contig analysis
    Quantif.unique.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption=caption, sample.name=sample.name.INT)
    
    ## Abundance analysis
    ### Plots
    for(x in list_type_clT) assign(paste0("plot_abundanceContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::abundanceContig(tcr.combined, cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
    ### Save
    png(paste0(global_output,'/abundanceContig.png'), width = 800, height = 300)
    (plot_abundanceContig_gene_nt | plot_abundanceContig_gene | plot_abundanceContig_nt | plot_abundanceContig_aa ) +
      plot_annotation(title = sample.name.INT, subtitle = paste0("(",dim(sobj)[2]," cells)"), caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
    dev.off()
    
    ## Contigs Length analysis
    for(x in list_type_contig){
      ### This should give multimodal plot
      assign(paste0("plot_lengthContig_",x,"_comb"), scRepertoire::lengthContig(tcr.combined, cloneCall=x, chains = "combined") + Seurat::NoLegend())
      ### Plots the A and B chains distribution separately
      assign(paste0("plot_lengthContig_",x,"_sin"), scRepertoire::lengthContig(tcr.combined, cloneCall=x, chains = "single") + Seurat::NoLegend())
    }
    ### Save
    png(paste0(global_output,'/lengthContig.png'), width = 800, height = 800)
    (plot_lengthContig_nt_comb | plot_lengthContig_nt_sin) / (plot_lengthContig_aa_comb | plot_lengthContig_aa_sin ) +
      plot_annotation(title = sample.name.INT, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
    dev.off()
    
    ## Clonal Homeostasis analysis
    Homeo.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name.INT)

    ## Clonal Proportions analysis
    Prop.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name.INT)
    
    ## Diversity analysis
    Div.g(combined = tcr.combined, list_type_clT = list_type_clT, out.dir = global_output, caption = caption, sample.name = sample.name.INT)
    

    ## Combine TCR data with seurat object
    ### Corresponding barcode
    for(i in 1:length(tcr.combined)) tcr.combined[[i]]$barcode <- gsub(pattern = "TCR", replacement = "GE", tcr.combined[[i]]$barcode)
    ### Combination par echantillon
    sobj <- scRepertoire::combineExpression(df = tcr.combined, sc = sobj, cloneCall="aa")
    sobj@meta.data$Frequency_indiv=sobj@meta.data$Frequency
    ### Combination tout echantillon confondu
    tcr.combined_unlist <- do.call("rbind", tcr.combined)
    sobj <- scRepertoire::combineExpression(df = tcr.combined_unlist, sc = sobj, cloneCall="aa")
    sobj@meta.data$Frequency_all=sobj@meta.data$Frequency
    sobj@meta.data$Frequency=NULL
    rm(tcr.combined,tcr.combined_unlist)
    ### Spliting CTstrict (into separate columns for TRA-V/J/C, TRB-V/J/C and corresponding sequences, with 2 possible clonotypes) and save as metadata
    ### and Adding length of TR sequence to meta.data
    sobj <- split.CTstrict.tcr(sobj)
    ###for all samples
    #### Plots
    for(x in c("TRAV_1","TRAJ_1","TRAC_1","TRBV_1","TRBJ_1","TRBC_1","TRAV_2","TRAJ_2","TRAC_2","TRBV_2","TRBJ_2","TRBC_2")) assign(paste0("dimplot_",x), tryCatch( {  print(Seurat::DimPlot(sobj, group.by = x, reduction = INT.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x))  },  error=function(err) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x)  } ))
    for(x in c("TRA_nt_1_len","TRA_nt_2_len","TRB_nt_1_len","TRB_nt_2_len")) assign(paste0("featureplot_",x), tryCatch( {  print(Seurat::FeaturePlot(sobj, features = x, reduction = INT.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x)))  },  error=function(err) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))  } ))
    #### Save plots
    png(paste0(global_output,'/cloneType_',sample.name.INT,'.png'), width =1400, height = 3000)
    ( wrap_elements( (dimplot_TRAV_1 / dimplot_TRAJ_1 / dimplot_TRAC_1 / dimplot_TRBV_1 / dimplot_TRBJ_1 / dimplot_TRBC_1 / featureplot_TRA_nt_1_len / featureplot_TRB_nt_1_len) +
                       plot_annotation(title = 'TR clonotype 1', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) |
        wrap_elements( (dimplot_TRAV_2 / dimplot_TRAJ_2 /dimplot_TRAC_2 / dimplot_TRBV_2 / dimplot_TRBJ_2 / dimplot_TRBC_2 / featureplot_TRA_nt_2_len / featureplot_TRB_nt_2_len) +
                         plot_annotation(title = 'TR clonotype 2', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) +
        plot_annotation(title = sample.name.INT, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) )
    dev.off()
    ### by samples
    for (i in seq(sample.name)){
      #### selection des data du sample
      cells_sample=rownames(sobj@meta.data[sobj@meta.data$orig.ident == sample.name.GE[i],])
      sub_sobj=subset(sobj, cells = cells_sample)
      #### Plots
      for(x in c("TRAV_1","TRAJ_1","TRAC_1","TRBV_1","TRBJ_1","TRBC_1","TRAV_2","TRAJ_2","TRAC_2","TRBV_2","TRBJ_2","TRBC_2")) assign(paste0("dimplot_",x), tryCatch( {  print(Seurat::DimPlot(sub_sobj, group.by = x, reduction = INT.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(x)) },  error=function(err) { patchwork::plot_spacer() + Seurat::DarkTheme() + ggplot2::ggtitle(x)  } ))
      for(x in c("TRA_nt_1_len","TRA_nt_2_len","TRB_nt_1_len","TRB_nt_2_len")) assign(paste0("featureplot_",x), tryCatch( {  print(Seurat::FeaturePlot(sub_sobj, features = x, reduction = INT.reduction) + Seurat::DarkTheme() + ggplot2::ggtitle(paste0("Nucleotidic length ", x))) },  error=function(err) { patchwork::plot_spacer() + Seurat::DarkTheme() +ggplot2::ggtitle(paste0("Nucleotidic length ", x))  } ))
      #### Save plots
      png(paste0(global_output,'/cloneType_',sample.name[i],'.png'), width =1400, height = 3000)
      print( wrap_elements( (dimplot_TRAV_1 / dimplot_TRAJ_1 / dimplot_TRAC_1 / dimplot_TRBV_1 / dimplot_TRBJ_1 / dimplot_TRBC_1 / featureplot_TRA_nt_1_len / featureplot_TRB_nt_1_len) +
                              plot_annotation(title = 'TR clonotype 1', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) |
               wrap_elements( (dimplot_TRAV_2 / dimplot_TRAJ_2 /dimplot_TRAC_2 / dimplot_TRBV_2 / dimplot_TRBJ_2 / dimplot_TRBC_2 / featureplot_TRA_nt_2_len / featureplot_TRB_nt_2_len) +
                                plot_annotation(title = 'TR clonotype 2', theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) ) +
               plot_annotation(title = sample.name[i], theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold"))) )
      dev.off()
    }
    
    ## Frequency analysis
    require(dplyr)
    sobj <- Freq.g(sobj=sobj, out.dir = global_output, sample.name=sample.name.INT, reduction=INT.reduction, freq_col="Frequency_all")
    for (i in seq(sample.name)){
      #### selection des data du sample
      cells_sample=rownames(sobj@meta.data[sobj@meta.data$orig.ident == sample.name.GE[i],])
      sub_sobj=subset(sobj, cells = cells_sample)
      #### Analysis
      #Top 10 frequencies, and top 10 to top 20 frequencies
      top20_freq = sub_sobj@meta.data %>% select(Frequency_indiv,CTaa,highlight_aa_all) %>% distinct() %>% arrange(desc(Frequency_indiv)) %>% na.omit() %>% top_n(n = 20, wt = Frequency_indiv)
      top20_freq = top20_freq[1:20,]
      rownames(top20_freq)=top20_freq$highlight_aa_all
      sub_sobj$highlight_aa_top10_freq <- ifelse(sub_sobj$highlight_aa_all %in% top20_freq$highlight[1:10], sub_sobj$highlight_aa_all, NA)
      sub_sobj$highlight_aa_top11to20_freq <- ifelse(sub_sobj$highlight_aa_all %in% top20_freq$highlight[11:length(top20_freq$highlight)], sub_sobj$highlight_aa_all, NA)
      #UMAP of top 10 frequencies
      png(paste0(global_output,'/Frequency_top_10_umap',sample.name[i],'.png'), width = 800, height = (400+350))
      print(patchwork::wrap_elements( (Seurat::DimPlot(sub_sobj, reduction = INT.reduction, group.by = "highlight_aa_top10_freq")  + Seurat::DarkTheme()) / gridExtra::tableGrob(top20_freq[1:10,c("Frequency_indiv","CTaa")], theme = gridExtra::ttheme_default(base_size = 10)) +
                                        plot_annotation(title = paste0(sample.name[i],": Top 10 Clonotypes (by frequencies)"),  subtitle = paste0("(",dim(sobj)[2]," cells)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0.5, face="bold"))) +
                                        plot_layout(heights = c(2, 1))))
      dev.off()
      #UMAP of top 11 to 20 frequencies
      png(paste0(global_output,'/Frequency_top11to20_umap',sample.name[i],'.png'), width = 800, height = (400+350))
      print(patchwork::wrap_elements( (Seurat::DimPlot(sub_sobj, reduction = INT.reduction, group.by = "highlight_aa_top11to20_freq")  + Seurat::DarkTheme()) / gridExtra::tableGrob(top20_freq[11:length(top20_freq$CTaa),c("Frequency_indiv","CTaa")], theme = gridExtra::ttheme_default(base_size = 10)) +
                                        plot_annotation(title = paste0(sample.name[i], ": Top 11 to 20 Clonotypes (by frequencies)"),  subtitle = paste0("(",dim(sobj)[2]," cells)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0.5, face="bold"))) +
                                        plot_layout(heights = c(2, 1))))
      dev.off()
    }
    
    ## Physicochemical properties of the CDR3
    Physicochemical_properties.g(sobj=sobj, list_type_clT = list_type_clT, out.dir = global_output, sample.name=sample.name.INT, type='TCR')

    
    
    ## CLUSTERS LEVEL ANALYSIS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    print("Clusters_analysis")
    clusters_output <- paste0(output_path_TCR, "Clusters_analysis")
    dir.create(path = clusters_output, recursive = TRUE, showWarnings = TRUE)
    
    for (i in seq(sample.name)){
      #create directory
      sample_output=paste0(clusters_output, "/", sample.name[i])
      dir.create(path = sample_output, recursive = TRUE, showWarnings = TRUE)
      
      #### selection des data du sample
      cells_sample=rownames(sobj@meta.data[sobj@meta.data$orig.ident == sample.name.GE[i],])
      sub_sobj=subset(sobj, cells = cells_sample)
      
      ## Filter cells that have no value for the x concerned + conversion to a list by clusters
      ## Need to filter, otherwise the functions count the 'NA' as a sequence.
      for(x in list_type_clT){
        if(x=="gene+nt") y="strict" else y=x
        filtred_sobj = sub_sobj[,!is.na(sub_sobj@meta.data[paste0("CT", y)])]
        assign(paste0("filtred_metadata_", sub("\\+","_",x)), scRepertoire::expression2List(sc = filtred_sobj, group = "seurat_clusters"))
      }
      rm(filtred_sobj)
      
      ## Quantification of unique contig analysis
      Quantif.unique.c(sobj = sub_sobj, ident.name=ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption=caption, sample.name=sample.name[i],filtred_metadata_aa=filtred_metadata_aa,filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)
      
      ## Abundance analysis
      ### Plots
      for(x in list_type_clT) assign(paste0("plot_cluster_abundanceContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::abundanceContig(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, scale = F) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) ))
      ### Save
      png(paste0(sample_output,'/abundanceContig',sample.name[i],'.png'), width = 2000, height = 600)
      print(plot_cluster_abundanceContig_gene_nt | plot_cluster_abundanceContig_gene | plot_cluster_abundanceContig_nt | plot_cluster_abundanceContig_aa ) +
        plot_annotation(title = sample.name[i], caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold")))
      dev.off()
      
      ## Clonal Homeostasis analysis
      sub_sobj=Homeo.c(sobj = sub_sobj, ident.name=ident.name, list_type_clT = list_type_clT, out.dir = sample_output, caption=caption, sample.name=sample.name[i], filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)
      
      ## Clonal Proportions analysis
      sub_sobj=Prop.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption=caption, sample.name=sample.name[i], filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)
      
      ## Diversity analysis
      sub_sobj=Div.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption=caption, sample.name=sample.name[i], filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)
      
      ## Frequency analysis
      Freq.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption=caption, sample.name=sample.name[i], ident.name=ident.name, reduction=INT.reduction, freq_col="Frequency_indiv", filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)
      
      ## Clonal Overlap analysis
      if(length(levels(Seurat::Idents(sobj)))!=1) Overlap.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption=caption, sample.name=sample.name[i], filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt)
      
      ### Physico-chemical properties of the CDR3
      Physicochemical_properties.c(sobj = sub_sobj, list_type_clT = list_type_clT, out.dir = sample_output, caption=caption, sample.name=sample.name[i], ident.name=ident.name, filtred_metadata_aa=filtred_metadata_aa, filtred_metadata_nt=filtred_metadata_nt, filtred_metadata_gene=filtred_metadata_gene, filtred_metadata_gene_nt=filtred_metadata_gene_nt, type='TCR')
    }
    
    #renamme TCR columns with 'TCR_' prefix
    toMatch <- c("^CTgene","^CTnt","^CTaa","^CTstrict","^Frequency","^cloneType","^TRAV_1","^TRAJ_1","^TRAC_1","^TRAV_2","^TRAJ_2","^TRAC_2","^TRA_nt_1","^TRA_nt_2","^TRBV_1","^TRBJ_1","^None_1","^TRBC_1","^TRBV_2","^TRBJ_2","^None_2","^TRBC_2","^TRB_nt_1","^TRB_nt_2","^TRA_nt_1_len","^TRA_nt_2_len","^TRB_nt_1_len","^TRB_nt_2_len","^highlight_aa_all","^highlight_aa_top")
    matches <- grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data))
    colnames(sobj@meta.data)[matches] = paste0("TCR_", grep(paste(toMatch,collapse="|"), colnames(sobj@meta.data), value=TRUE))
    
    ### Saving GE_ADT_TCR object
    GE_TCR_file = paste0(output_path_TCR, basename(GE_file), '_TCR')
    save(sobj, file = paste0(GE_TCR_file, '.rda'), compress = "bzip2")
    print("TCR done!")
    
  }else{
    GE_TCR_file = paste0(output_path_TCR, basename(GE_file), '_TCR')
    load(paste0(GE_TCR_file, '.rda'))
    
    ### Building cerebro binary
    seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = FALSE, clusters.colnames = ident.name, species = species.rdx, gmt.file = gmt.file, file = GE_TCR_file, nthreads = nthreads, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')
    seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = FALSE, clusters.colnames = ident.name, species = species.rdx, gmt.file = gmt.file, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, file = GE_TCR_file, nthreads = nthreads,  min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')
    #TCR view
    #Fusion of *highlight_aa_top10_freq TCR_highlight_aa_top11to20_freq *highlight_aa_top11to20_freq (only TCR)
    col.names.top <- grep("highlight_aa_top", names(sobj@meta.data), value=TRUE)
    col.name.all <- grep("highlight_aa_all", names(sobj@meta.data), value=TRUE)
    for (i in 1:length(sobj@meta.data[[col.name.all]])){
      if (is.na(sobj@meta.data[[col.names.top[1]]][i])){
        if (is.na(sobj@meta.data[[col.names.top[2]]][i])){
          sobj@meta.data$TCR_highlight_aa_top20_freq[i] <- NA
        }else{
          sobj@meta.data$TCR_highlight_aa_top20_freq[i] <- sobj@meta.data[[col.names.top[2]]][i]
        }
      }else if (is.na(sobj@meta.data[[col.names.top[2]]][i])){
        sobj@meta.data$TCR_highlight_aa_top20_freq[i] <- sobj@meta.data[[col.names.top[1]]][i]
      }else{
          stop(paste0("Error :", col.names.top[1]," AND ", col.names.top[2]," are full in ",i,"! "))
        }
    }

    seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = TRUE, clusters.colnames = "TCR_highlight_aa_top20_freq", species = species.rdx, gmt.file = gmt.file, file = GE_TCR_file, nthreads = nthreads, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')
    seurat2cerebro(sobj = sobj, ident = ident.name, other.reductions = FALSE, other.idents = TRUE, clusters.colnames = "TCR_highlight_aa_top20_freq", species = species.rdx, gmt.file = gmt.file, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, file = GE_TCR_file, nthreads = nthreads,  min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox')
  }
 
 }