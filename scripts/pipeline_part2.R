#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.ge", help="Input unfiltred seurat object (in .rda format)."),
  make_option("--output.dir.ge", help="Output path"),
  make_option("--author.name", help="Name of author of the analysis"),
  make_option("--author.mail", help="Email of author of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder"),
  ### Analysis Parameters
  # QC cell
  make_option("--pcmito.min", help="Threshold min for percentage of mitochondrial RNA (below this threshold the cells are eliminated)"),
  make_option("--pcmito.max", help="Threshold max for percentage of mitochondrial RNA (above this threshold the cells are eliminated)"),
  make_option("--pcribo.min", help="Threshold min for percentage of ribosomal RNA (below this threshold the cells are eliminated)"),
  make_option("--pcribo.max", help="Threshold max for percentage of ribosomal RNA (above this threshold the cells are eliminated)"),
  make_option("--min.features", help="Threshold min for number of genes (below this threshold the cells are eliminated)"),
  make_option("--min.counts", help="Threshold min for number of UMI (below this threshold the cells are eliminated)"),
  # QC gene
  make_option("--min.cells", help="Include genes expressed in at least this many cells (minimum cells covering)"),
  # Doublets
  make_option("--doublets.filter.method", help="Method used to filter doublets (scDblFinder, scds, all). To not filter set the parameter to none."),
  ### Databases
  ##Metadata
  make_option("--metadata.file", help="csv file with the metadata to add in the seurat object"),
  # QC
  make_option("--cc.seurat.file", help="RDS file with list of cell cycle genes for seurat"),
  make_option("--cc.cyclone.file", help="RDS file with list of cell cycle genes for cyclone")
)
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)


### Sourcing functions
if(is.null(args$options$pipeline.path)) stop("--pipeline.path parameter must be set!")
source(paste0(args$options$pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))


#### Get Paramaters ####
### Formatting: convert "NULL"/"FALSE"/"TRUE" (in character) into NULL/FALSE/TRUE
args <- convert_NULL_BOOL(args)
### Project
list.author.name <- splitByComma.ifnotNULL(author.name)
list.author.mail <- splitByComma.ifnotNULL(author.mail)
### Computational Parameters
nthreads <- asNumeric.ifnotNULLelse(nthreads, 1)
### Analysis Parameters
# QC cell
pcmito.min <- asNumeric.ifnotNULLelse(pcmito.min, 0)
pcmito.max <- asNumeric.ifnotNULLelse(pcmito.max, 20)
pcribo.min <- asNumeric.ifnotNULLelse(pcribo.min, 0)
pcribo.max <- asNumeric.ifnotNULLelse(pcribo.max, 100)
min.features <- asNumeric.ifnotNULLelse(min.features, 200)
min.counts <- asNumeric.ifnotNULLelse(min.counts, 1000)
# QC gene
min.cells <- asNumeric.ifnotNULLelse(min.cells, 5)
# Doublets
if (is.null(doublets.filter.method)) doublets.filter.method <- "all"
### Databases
metadata.file <- splitByComma.ifnotNULL(metadata.file)
### Fixed parameters
assay <- 'RNA'
### Cleaning
rm(option_list,parser,args)

#### Check non-optional parameters ####
if (is.null(input.rda.ge)) stop("input.rda.ge parameter can't be empty!")
if (is.null(output.dir.ge)) stop("output.dir.ge parameter can't be empty!")

### Load data
load(input.rda.ge)

#### Get Missing Paramaters ####
### Project
sample.name.ge <- sobj@misc$params$sample.name.ge
species <- sobj@misc$params$species
### Databases
if (species == "homo_sapiens") {
  if (is.null(cc.seurat.file)) cc.seurat.file <-  paste0(pipeline.path,"/resources/GENELISTS/homo_sapiens_Seurat_cc.genes_20191031.rds")
  if (is.null(cc.cyclone.file)) cc.cyclone.file <- paste0(pipeline.path,"/resources/GENELISTS/homo_sapiens_cyclone_pairs_symbols_20191001.rds")
}
if (species == "mus_musculus") {
  if (is.null(cc.seurat.file)) cc.seurat.file <- paste0(pipeline.path,"/resources/GENELISTS/mus_musculus_Seurat_cc.genes_20191031.rds")
  if (is.null(cc.cyclone.file)) cc.cyclone.file <- paste0(pipeline.path,"/resources/GENELISTS/mus_musculus_cyclone_pairs_symbols_20191015.rds")
}


################################################################################
## MAIN: HALF-FILTERED ANALYSIS + DOUBLETS ANALYSIS
################################################################################

### printing parameters:
print("###########################################")
print(paste0("sample : ",sample.name.ge))
print(paste0("input.rda.ge : ",input.rda.ge))
print(paste0("output.dir.ge : ",output.dir.ge))
print("###########################################")

### Seurat multithreading
RcppParallel::setThreadOptions(numThreads = nthreads)

### Add authors and metadata
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)
if(!is.null(metadata.file)) sobj <- add_metadata_sobj(sobj=sobj, metadata.file = metadata.file)

### Filtering cells
cat("\nFiltering cells...\n")
sobj <- cells.QC.filter(sobj = sobj, min.features = min.features, min.counts = min.counts, pcmito.range = c(pcmito.min, pcmito.max), pcribo.range = c(pcribo.min, pcribo.max))
sobj@misc$excel$After_QC_cells_filtering$estim_cells <- ncol(sobj)

### Cell cycle prediction
cat("\nCell cycle prediction...\n")
sobj <- cell.cycle.predict(sobj = sobj, assay = assay, cc.cyclone.file = cc.cyclone.file, cc.seurat.file = cc.seurat.file, nbin = 10, nthreads = nthreads)

### Filtering features (based on minimum cells covering)
cat("\nFiltering features...\n")
sobj <- features.filter(sobj = sobj, min.cells = min.cells)

### Identification of doublets
cat("\nIdentification of doublets...\n")
sobj <- find.doublets(sobj = sobj, assay = assay)

### Building output directory
filtered.dir <- paste0(output.dir.ge, paste0('F', min.features, '_C', min.counts, '_M', paste(c(pcmito.min, pcmito.max), collapse = '-'), '_R', paste(c(pcribo.min, pcribo.max), collapse = '-'),'_G', min.cells))


################################################################################
## MAIN:  DOUBLETS KEPT
################################################################################
cat("\nDoublets kept analysis\n")

### Building doublet kept output directory
doublet.kept.dir <- paste0(filtered.dir, '/DOUBLETSKEPT/')
dir.create(path = doublet.kept.dir, recursive = TRUE, showWarnings = FALSE)

### Computing basic metrics : percentage of counts in the top features + mito + ribo + stress + nb features + nb counts
cat("\nComputing QC metrics with kept doublets...\n")
sobj <- QC.metrics(sobj = sobj, assay = assay, pcmito.range = c(pcmito.min, pcmito.max), pcribo.range = c(pcribo.min, pcribo.max), min.features = min.features, min.counts = min.counts)

### QC histograms
QC.hist(sobj = sobj, assay = assay, out.dir = doublet.kept.dir)

### Saving stat for sheet
if(doublets.filter.method == 'none'){
  message("Saving stat...\n")
  sobj@misc$excel$Final_Cells_Quality$nb_genes <- dim(sobj)[1]
  sobj@misc$excel$Final_Cells_Quality$nb_cells <- dim(sobj)[2]
  sobj@misc$excel$QC_genomic_core_facility$nb_reads_per_cells_CellRanger_like <- sobj@misc$excel$Kallisto_Bustools_alignment$Total_reads / sobj@misc$excel$Droplet_Quality$estimated_cells
  sobj@misc$excel$QC_genomic_core_facility$nb_reads_per_cells_final <- sobj@misc$excel$Kallisto_Bustools_alignment$Total_reads / sobj@misc$excel$Final_Cells_Quality$nb_cells
  sobj@misc$params$doublets$method_filtering ="none"
  save_stat(sobj = sobj, sample.name = sample.name.ge, title = sample.name.ge, out.dir = doublet.kept.dir)
  ### Materials and Methods
   sobj@misc$parameters$Materials_and_Methods$part2_Filtering_droplets <- paste0("The count matrix was filtered to exclude genes detected in less than ",min.cells," cells, cells with less than ",min.counts," UMIs or less than ",min.features," detected genes, as well as cells with mitochondrial transcript proportions higher than ",pcmito.max,"%",
                                                              if(!pcmito.min==0) paste0("and less than ",pcmito.min, "%"),".",
                                                              if(!pcribo.max==1 && !pcribo.min==0) paste0("as well as cells with ribosomal transcripts proportion higher than ",pcribo.max,"% and less than ",pcribo.min, "%. The proportion of mechanical stress-response gene counts (Thesis of Leo Machado entitled 'From skeletal muscle stem cells to tissue atlases: new tools to investigate and circumvent dissociation-induced stress', 2019) were also estimated but not used to filter cells.") else "The proportion of ribosomal gene counts and the proportion of mechanical stress-response gene counts (Thesis of Leo Machado entitled 'From skeletal muscle stem cells to tissue atlases: new tools to investigate and circumvent dissociation-induced stress', 2019) were also estimated but not used to filter cells.",
"Cell cycle scoring of each cell was performed using two methods : the CellcycleScoring() function from the Seurat package (version ",sobj@misc$technical_info$Seurat,"), and the cyclone() function from Scran (version ",sobj@misc$technical_info$scran,").
Barcodes corresponding to doublet cells were identified using two methods: scDblFinder (version ",sobj@misc$technical_info$scDblFinder,") using default parameters, and scds (version ",sobj@misc$technical_info$scds,") with its ", sobj@misc$params$doublets$scds_method," method using default parameters. However boublets were not discarded.")
  sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
}

### Saving doublets kept non-normalized object
cat("\nSaving object...\n")
save(sobj, file = paste0(doublet.kept.dir, '/', sample.name.ge, '_DOUBLETSKEPT_NON-NORMALIZED.rda'), compress = "bzip2")

### Basic normalization, dimension reduction, clustering, technical plots and save
cat("\nComputing complete analysis with doublets kept...\n")
norm.red.plot.quick(sobj = sobj, sample.name.ge = sample.name.ge, pre.out.dir = doublet.kept.dir, file.name = 'DOUBLETSKEPT')


################################################################################
## MAIN:  DOUBLETS REMOVED
################################################################################

if(doublets.filter.method != 'none'){
  cat("\nDoublets removed analysis\n")

  ### Building doublet filter output directory
  doublet.filtered.dir <-  paste0(filtered.dir, '/DOUBLETSFILTER_', doublets.filter.method, '/')
  dir.create(path = doublet.filtered.dir, recursive = TRUE, showWarnings = FALSE)
  
  ### Filter doublets
  cat("\nFiltering doublets...\n")
  sobj <- filter.doublets(sobj = sobj, method = doublets.filter.method)
  
  ### Computing basic metrics : percentage of counts in the top features + mito + ribo + stress + nb features + nb counts
  cat("\nComputing QC metrics with removed doublets...\n")
  sobj <- QC.metrics(sobj = sobj, assay = assay, pcmito.range = c(pcmito.min, pcmito.max), pcribo.range = c(pcribo.min, pcribo.max), min.features = min.features, min.counts = min.counts)

  ### QC histograms
  QC.hist(sobj = sobj, assay = assay, out.dir = doublet.filtered.dir)

  ### Saving stat for sheet
  cat("\nSaving stat...\n")
  sobj@misc$excel$Final_Cells_Quality$nb_genes <- dim(sobj)[1]
  sobj@misc$excel$Final_Cells_Quality$nb_cells <- dim(sobj)[2]
  sobj@misc$excel$QC_genomic_core_facility$nb_reads_per_cells_CellRanger_like <- sobj@misc$excel$Kallisto_Bustools_alignment$Total_reads / sobj@misc$excel$Droplet_Quality$estimated_cells
  sobj@misc$excel$QC_genomic_core_facility$nb_reads_per_cells_final <- sobj@misc$excel$Kallisto_Bustools_alignment$Total_reads / sobj@misc$excel$Final_Cells_Quality$nb_cells
  save_stat(sobj = sobj, sample.names = sample.name.ge, title = sample.name.ge, out.dir = doublet.filtered.dir)

  ### Materials and Methods
  sobj@misc$parameters$Materials_and_Methods$part2_Filtering <- paste0("The count matrix was filtered to exclude genes detected in less than ",min.cells," cells, cells with less than ",min.counts," UMIs or less than ",min.features," detected genes, as well as cells with mitochondrial transcript proportions higher than ",pcmito.max,"%",
                                                              if(!pcmito.min==0) paste0(" and less than ",pcmito.min, "%"),
                                                              if(!pcribo.max==100 || !pcribo.min==0) paste0(" as well as cells with ribosomal transcripts proportion higher than ",pcribo.max,"% and less than ",pcribo.min,"%. The proportion of mechanical stress-response gene counts (Thesis of Leo Machado entitled 'From skeletal muscle stem cells to tissue atlases: new tools to investigate and circumvent dissociation-induced stress', 2019) were also estimated but not used to filter cells.") else ". The proportion of ribosomal gene counts and the proportion of mechanical stress-response gene counts (Thesis of Leo Machado entitled 'From skeletal muscle stem cells to tissue atlases: new tools to investigate and circumvent dissociation-induced stress', 2019) were also estimated but not used to filter cells. ",
                                                              "Cell cycle scoring of each cell was performed using two methods: the CellcycleScoring() function from the Seurat package (version ",sobj@misc$technical_info$Seurat,"), and the cyclone() function from Scran (version ",sobj@misc$technical_info$scran,"). ",
                                                              if(doublets.filter.method == 'all') paste0("Barcodes corresponding to doublet cells were identified and discarded using the union of two methods: scDblFinder (version ",sobj@misc$technical_info$scDblFinder,") using default parameters, and scds (version ",sobj@misc$technical_info$scds,") with its ", sobj@misc$params$doublets$scds_method," method using default parameters. We manually verified that the cells identified as doublets did not systematically correspond to cells in G2M phase.") else if(doublets.filter.method == 'scDblFinder') paste0("Barcodes corresponding to doublet cells were identified and discarded using scDblFinder (version ",sobj@misc$technical_info$scDblFinder,") using default parameters. We manually verified that the cells identified as doublets did not systematically correspond to cells in G2M phase.") else if(doublets.filter.method == 'scds') paste0("Barcodes corresponding to doublet cells were identified and discarded using scds (version ",sobj@misc$technical_info$scds,") with its ", sobj@misc$params$doublets$scds_method," method using default parameters. We manually verified that the cells identified as doublets did not systematically correspond to cells in G2M phase."))
  sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
  
  ### Saving non-normalized object
  cat("\nSaving object...\n")
  save(sobj, file = paste0(doublet.filtered.dir, sample.name.ge, '_FILTERED_NON-NORMALIZED.rda'), compress = "bzip2")
}
