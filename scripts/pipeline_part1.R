#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.dir.ge", help="Input path of the alignment"),
  make_option("--output.dir.ge", help="Output path"),
  make_option("--sample.name.ge", help="Name of RNA-seq sample"),
  make_option("--species", help="Species of sample ( Only 'homo_sapiens', 'mus_musculus' supported yet)"),
  make_option("--author.name", help="Name of author of the analysis"),
  make_option("--author.mail", help="Email of author of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder"),
  ### Analysis Parameters
  # Emptydrops
  make_option("--emptydrops.fdr", help="FDR threshold for emptydrops tool"),
  make_option("--droplets.limit", help="Number min of droplets to run emptydrops"),
  make_option("--emptydrops.retain", help="All droplets above this value is considered as a cell"),
  # QC cell
  make_option("--pcmito.min", help="Threshold min for percentage of mitochondrial RNA (below this threshold the cells are eliminated)"),
  make_option("--pcmito.max", help="Threshold max for percentage of mitochondrial RNA (above this threshold the cells are eliminated)"),
  make_option("--pcribo.min", help="Threshold min for percentage of ribosomal RNA (below this threshold the cells are eliminated)"),
  make_option("--pcribo.max", help="Threshold max for percentage of ribosomal RNA (above this threshold the cells are eliminated)"),
  make_option("--min.features", help="Threshold min for number of genes (below this threshold the cells are eliminated)"),
  make_option("--min.counts", help="Threshold min for number of UMI (below this threshold the cells are eliminated)"),
  # QC gene
  make_option("--min.cells", help="Include genes expressed in at least this many cells (minimum cells covering)"),
  ### Databases
  # Metadata
  make_option("--metadata.file", help="csv file with the metadata to add in the seurat object"),
  # QC
  make_option("--mt.genes.file", help="RDS file with list of mitochondrial genes"),
  make_option("--crb.genes.file", help="RDS file with list of ribosomal genes"),
  make_option("--str.genes.file", help="RDS file with list of stress genes")
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
if (is.null(species)) species <- "homo_sapiens"
list.author.name <- splitByComma.ifnotNULL(author.name)
list.author.mail <- splitByComma.ifnotNULL(author.mail)
### Computational Parameters
nthreads <- asNumeric.ifnotNULLelse(nthreads, 1)
### Analysis Parameters
# Emptydrops
emptydrops.fdr <-  asNumeric.ifnotNULLelse(emptydrops.fdr, 1E-03)
droplets.limit <- asNumeric.ifnotNULLelse(droplets.limit, 1E+05)
emptydrops.retain <- asNumeric.ifnotNULLelse(emptydrops.retain)
# QC cell
pcmito.min <- asNumeric.ifnotNULLelse(pcmito.min, 0)
pcmito.max <- asNumeric.ifnotNULLelse(pcmito.max, 20)
pcribo.min <- asNumeric.ifnotNULLelse(pcribo.min, 0)
pcribo.max <- asNumeric.ifnotNULLelse(pcribo.max, 100)
min.features <- asNumeric.ifnotNULLelse(min.features, 200)
min.counts <- asNumeric.ifnotNULLelse(min.counts, 1000)
# QC gene
min.cells <- asNumeric.ifnotNULLelse(min.cells, 5)
### Databases
metadata.file <- splitByComma.ifnotNULL(metadata.file)
if (species == "homo_sapiens") {
  if (is.null(mt.genes.file)) mt.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/homo_sapiens_mito_symbols_20191001.rds")
  if (is.null(crb.genes.file)) crb.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/homo_sapiens_cribo_symbols_20191015.rds")
  if (is.null(str.genes.file)) str.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/homo_sapiens_stress_symbols_20200224.rds")
}
if (species == "mus_musculus") {
  if (is.null(mt.genes.file)) mt.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/mus_musculus_mito_symbols_20191015.rds")
  if (is.null(crb.genes.file)) crb.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/mus_musculus_cribo_symbols_20191015.rds")
  if (is.null(str.genes.file)) str.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/mus_musculus_stress_symbols_20200224.rds")
}
### Fixed parameters
my.seed <- 1337L
assay <- 'RNA'
### Cleaning
rm(option_list,parser,args)


#### Check non-optional parameters ####
if (is.null(input.dir.ge)) stop("input.dir.ge parameter can't be empty!")
if (is.null(output.dir.ge)) stop("output.dir.ge parameter can't be empty!")
if (is.null(sample.name.ge)) stop("sample parameter can't be empty!")


################################################################################
## MAIN: QC_droplets ANALYSIS
################################################################################

### Printing parameters
print("###########################################")
print(paste0("sample.name.ge : ", sample.name.ge))
print(paste0("input.dir.ge : ", input.dir.ge))
print(paste0("output.dir.ge : ", output.dir.ge))
print("###########################################")

### Creating parallel instance and Seurat multithreading
cl <- create.parallel.instance(nthreads = nthreads)

### Building output directory
unfiltred.dir <- paste0(output.dir.ge, 'QC_droplets', if(!is.null(emptydrops.retain)) '_retain', emptydrops.retain, '/')
dir.create(path = unfiltred.dir, recursive = TRUE, showWarnings = FALSE)

### Loading raw count matrix + Filtering duplicated cell barcodes + Removing empty droplets
sobj <- load.sc.data(data.path = input.dir.ge, sample.name = sample.name.ge, assay = assay, droplets.limit = droplets.limit, emptydrops.fdr = emptydrops.fdr, emptydrops.retain = emptydrops.retain, BPPARAM = cl, my.seed = my.seed, out.dir = unfiltred.dir, min.count = min.count)

### Add metadata
if(!is.null(metadata.file)) sobj <- add_metadata_sobj(sobj=sobj, metadata.file = metadata.file)

### Save project parameters
sobj@misc$params$sample.name.ge <- sample.name.ge
sobj@misc$params$species <- species
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = author.name, list.author.mail = author.mail)
sobj@misc$params$analysis_type <- "Individual analysis"

### Computing basic metrics : percentage of counts in the top features + mito + ribo + stress + nb features + nb counts
sobj <- QC.metrics(sobj = sobj, assay = assay, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, pcmito.range = c(pcmito.min, pcmito.max), pcribo.range = c(pcribo.min, pcribo.max), min.features = min.features, min.counts = min.counts)

### QC histograms
QC.hist(sobj = sobj, assay = assay, out.dir = unfiltred.dir)

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$part1_Droplets_QC <- paste0("Cell barcode by symbol count table were loaded in R (version ", getRversion(), ") using the BUSpaRse package (version ", sobj@misc$technical_info$BUSpaRse,"). To call real cells from empty droplets, we used the emptyDrops() function from the dropletUtils package (version ", sobj@misc$technical_info$DropletUtils,"), which assesses whether the RNA content associated with a cell barcode is significantly distinct from the ambient background RNA present within each sample. Barcodes with p-value < ", emptydrops.fdr," (Benjamini-Hochberg-corrected) were considered as legitimate cells for further analysis.")
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)

### Saving non-normalized object
save(sobj, file = paste0(unfiltred.dir, sample.name.ge, '_QC_NON-NORMALIZED.rda'), compress = "bzip2")
