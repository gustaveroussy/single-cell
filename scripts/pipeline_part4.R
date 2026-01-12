#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.ge", help="Input filtred, normalized and reducted seurat object (in .rda format)."),
  make_option("--output.dir.ge", help="Output path"),
  make_option("--author.name", help="Name of author of the analysis"),
  make_option("--author.mail", help="Email of author of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder"),
  ### Analysis Parameters
  # Metadata
  make_option("--metadata.file", help="csv file with the metadata to add in the seurat object"),
  # Clustering
  make_option("--keep.dims", help="Number of dimension to keep for clustering (from 0 to keep.dims)"),
  make_option("--keep.res", help="Resolution value for clustering"),
  # Annotation
  make_option("--custom.sce.ref", help="List of .RData files containing SingleCellExpriment objects with your reference"),
  make_option("--custom.markers.ref", help="List of .xlsx files containing your reference"),
  make_option("--cfr.minscore", help="Minimum correlation score for clustifyr to consider"),
  make_option("--sr.minscore", help="Minimum correlation score for SingleR to consider"),
  #skip steps
  make_option("--skip.technical_plots", help="Allow to skip the plotting of thechnical biases on umap"),
  make_option("--skip.annotation", help="Allow to skip the automatic annotation step"),
  make_option("--skip.markers_identification", help="Allow to skip the identification of marker genes for each cluster"),
  # Markerfile: umap + violin
  make_option("--markfile", help="Genes to plot on umap (format: 2 columns named Genes and Signatures)"),
  make_option("--markers.pt.size", help="Adjust point size to plot genes from the markfile"),
  make_option("--markers.order", help="Boolean determining whether to plot cells in order of expression of markfile genes (can be useful if cells expressing given gene are getting buried).")
)
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)


### Sourcing functions ####
if(is.null(args$options$pipeline.path)) stop("--pipeline.path parameter must be set!")
source(paste0(args$options$pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))


#### Get Paramaters ####
### Formatting and extract from args
args <- convert_NULL_BOOL(args)
### Project
list.author.name <- splitByComma.ifnotNULL(author.name)
list.author.mail <- splitByComma.ifnotNULL(author.mail)
### Computational Parameters
nthreads <- asNumeric.ifnotNULLelse(nthreads, 1)
### Analysis Parameters
# Clustering
keep.dims <- asNumeric.ifnotNULLelse(keep.dims, NULL)
keep.res <- asNumeric.ifnotNULLelse(keep.res, NULL)
# Annotation
custom.sce.ref <- splitByComma.ifnotNULL(custom.sce.ref)
custom.markers.ref <- splitByComma.ifnotNULL(custom.markers.ref)
cfr.minscore <- asNumeric.ifnotNULLelse(cfr.minscore, 0.35)
sr.minscore <- asNumeric.ifnotNULLelse(sr.minscore, 0.25)
# Skip steps
if (is.null(skip.technical_plots)) skip.technical_plots <- FALSE
if (is.null(skip.annotation)) skip.annotation <- FALSE
if (is.null(skip.markers_identification)) skip.markers_identification <- FALSE
# Markerfile: umap + violin
markfile <- splitByComma.ifnotNULL(markfile)
markers.pt.size <- asNumeric.ifnotNULLelse(markers.pt.size, 2)
if (is.null(markers.order)) markers.order <- FALSE
### Metadata
metadata.file <- splitByComma.ifnotNULL(metadata.file)
#### Fixed parameters
solo.pt.size <- 3
multi.pt.size <- 2
gradient.cols <- c("gold", "blue")
### Clean
rm(option_list,parser,args)

#### Check non-optional parameters ####
if (is.null(input.rda.ge)) stop("input.rda.ge parameter can't be empty!")
if (is.null(output.dir.ge)) stop("output.dir.ge parameter can't be empty!")
if (is.null(keep.dims)) stop("keep.dims parameter can't be empty!")
if (is.null(keep.res)) stop("keep.res parameter can't be empty!")

### Load data
load(input.rda.ge)

#### Get Missing Paramaters ####
sample.name.ge <- sobj@misc$params$sample.name.ge
species <- sobj@misc$params$species
assay <- Seurat::DefaultAssay(sobj)
dimred.method <- sobj@assays[[assay]]@misc$params$reductions$method
if (species == "homo_sapiens") {
  singler.setnames <- c("HumanPrimaryCellAtlasData", "BlueprintEncodeData", "NovershternHematopoieticData", "DatabaseImmuneCellExpressionData", "MonacoImmuneData")
  clustifyr.setnames <- c("pbmc_avg", "ref_hema_microarray", "ref_cortex_dev","ref_pan_indrop") # ref_hema_microarray same as hema_microarray_matrix
  scrnaseq.setnames <- c("BaronPancreasData(human)","MuraroPancreasData","SegerstolpePancreasData")
}
if (species == "mus_musculus") {
  singler.setnames <- c("MouseRNAseqData", "ImmGenData")
  clustifyr.setnames <- c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq")
  scrnaseq.setnames <- c("BaronPancreasData(mouse)","ZeiselBrainData") #,"ShekharRetinaData"
}
  
################################################################################
## MAIN: CLUSTERING, ANNOTATION, MARKER GENES
################################################################################

### printing parameters
print("###########################################")
print(paste0("Sample : ",sample.name.ge))
print(paste0("RDA file : ",input.rda.ge))
print(paste0("Output directory : ", output.dir.ge))
print(paste0("Dimension: ", keep.dims))
print(paste0("Resolution: ", keep.res))
print("###########################################")

### Creating parallel instance and Seurat multithreading
RcppParallel::setThreadOptions(numThreads = nthreads) #PCA, scaling
future::plan("multicore", workers = nthreads)  #FindMarkers, SCTransform, integration
#options(future.globals.maxSize = 8 * 1024^3)  #for a Seurat object of 8Go (avoid error like: "Error: The total size of the global variables exported to the workers; for future expression (‘globals’) exceeds the maximum allowed size.")
#options(future.rng.onMisuse = "ignore")       #remove useless warnings
options(SeuratCommand.umap.threads = nthreads) #for RunUMAP() (replace the n.thread option of the function)

### Add authors and metadata
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)
if(!is.null(metadata.file)) sobj <- add_metadata_sobj(sobj=sobj, metadata.file = metadata.file)

### Building clustered output directory
clust.dir <- paste(output.dir.ge, paste0("dims", keep.dims, "_res", keep.res), sep = '/')
dir.create(path = clust.dir, recursive = TRUE, showWarnings = FALSE)

### Replotting final clusters
cat("\nClustering...\n")
sobj <- louvain.cluster(sobj = sobj, reduction = paste0(assay, "_", dimred.method), max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size)

### Setting ident name
ident.name <- paste0(paste0(assay, "_", dimred.method, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))

### Technical plots
if (!skip.technical_plots){
    cat("\nSaving technical plots...\n")
    technical.plot(sobj = sobj, ident = ident.name, out.dir = clust.dir, multi.pt.size = multi.pt.size)
}

### Normalize RNA assay for DEG
sobj <- Seurat::NormalizeData(sobj, normalization.method = 'LogNormalize', assay = "RNA")

### Finding markers
if (!skip.markers_identification){
    cat("\nFinding markers...\n")
    sobj <- find.markers.quick(sobj = sobj, ident = ident.name, test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, adjp.p.max = 5E-02, topn = 10, heatmap.cols = gradient.cols, out.dir = clust.dir)
}

### Automatic cell type annotation
if (!skip.annotation){
    cat("\nAutomatic cell type annotation...\n")
    sobj <- cells.annot(sobj = sobj, ident = ident.name, singler.setnames = singler.setnames, clustifyr.setnames = clustifyr.setnames, scrnaseq.setnames = scrnaseq.setnames, custom.sce.ref = custom.sce.ref, custom.markers.ref = custom.markers.ref, sr.minscore = sr.minscore, cfr.minscore = cfr.minscore, out.dir = clust.dir, solo.pt.size = solo.pt.size, nthreads = nthreads)
}

### Assessing clusters : Plotting provided marker genes
cat("\nPlotting provided marker genes...\n")
if (!is.null(markfile)){
  markers <- get.markers.from.markersfiles(markfiles = markfile)
  sobj <- markers.umap.plot(sobj = sobj, markers = markers, ident = ident.name, out.dir = clust.dir, dimplot.cols = gradient.cols, markers.pt.size = markers.pt.size, markers.order = markers.order)
}

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$part4_Clust_Markers_Annot <- paste0(
  if (!skip.annotation) "An automatic annotation of cell types was perfom by SingleR (version ",sobj@misc$technical_info$SingleR,") (with fine-tuning step) and ClustifyR (version ",sobj@misc$technical_info$clustifyr,"), using packages built-in references. It labels clusters (or cells) from a dataset based on similarity (Spearman correlation score) to a reference dataset with known labels. The labels with a correlation score greater than ",sr.minscore," for SingleR or greater than ",cfr.minscore," for ClustifyR were kept. The annotation was also made CelliD (version ",sobj@misc$technical_info$CelliD,") with genes signatures from pangloa database. ",
  if (!skip.markers_identification) "Marker genes for Louvain clusters were identified through a 'one versus others' differential analysis using the Wilcoxon test through the FindAllMarkers() function from Seurat, considering only genes with a minimum log fold-change of 0.5 in at least 75% of cells from one of the groups compared, and FDR-adjusted p-values <0.05 (Benjaminin-Hochberg method)."
)
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj=sobj, output.dir=clust.dir)

### Saving final object
cat("\nSaving object...\n")
GE_file <- paste0(clust.dir, '/', paste(c(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input.rda.ge)), keep.dims, keep.res), collapse = '_'))
save(sobj, file = paste0(GE_file, '.rda'), compress = "bzip2")
