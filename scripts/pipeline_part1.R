#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.dir.ge", help="Input path of the alignment"),
  make_option("--output.dir.ge", help="Output path"),
  make_option("--sample.name.ge", help="Name of RNA-seq sample"),
  make_option("--species", help="Species of sample ( Only 'homo_sapiens', 'mus_musculus' supported yet)"),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  # Emptydrops
  make_option("--emptydrops.fdr", help="FDR threshold for emptydrops tool"),
  make_option("--droplets.limit", help="Number min of droplets to run emptydrops"),
  make_option("--emptydrops.retain", help="All droplets above this value is considered as a cell"),
  # Translate ENSG into Gene Symbol
  make_option("--translation", help="TRUE for translate ENSG into Gene Symbol"),
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
  # QC
  make_option("--mt.genes.file", help="RDS file with list of mitochondrial genes"),
  make_option("--crb.genes.file", help="RDS file with list of ribosomal genes"),
  make_option("--str.genes.file", help="RDS file with list of stress genes"),
  # Translation into gene Symbols
  make_option("--translation.file", help="Text file with correspondance ENSG / genes symbol"),
  # Yaml parameters file to remplace all parameters before (to use R script without snakemake)
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
input.dir.ge <- args$options$input.dir.ge
output.dir.ge <- args$options$output.dir.ge
sample.name.GE <- args$options$sample.name.ge
species <- args$options$species
author.name <- args$options$author.name
author.mail <- args$options$author.mail
### Computational Parameters
nthreads <- if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Emptydrops
emptydrops.fdr <-  if (!is.null(args$options$emptydrops.fdr) && (args$options$emptydrops.fdr != "NULL")) as.numeric(args$options$emptydrops.fdr)
droplets.limit <- if (!is.null(args$options$droplets.limit) && (args$options$droplets.limit != "NULL")) as.numeric(args$options$droplets.limit)
emptydrops.retain <- if (!is.null(args$options$emptydrops.retain) && (args$options$emptydrops.retain != "NULL")) as.numeric(args$options$emptydrops.retain)
# Translate ENSG into Gene Symbol
translation <- args$options$translation
# QC cell
pcmito.min <- if (!is.null(args$options$pcmito.min) && (args$options$pcmito.min != "NULL")) as.numeric(args$options$pcmito.min)
pcmito.max <- if (!is.null(args$options$pcmito.max) && (args$options$pcmito.max != "NULL")) as.numeric(args$options$pcmito.max)
pcribo.min <- if (!is.null(args$options$pcribo.min) && (args$options$pcribo.min != "NULL")) as.numeric(args$options$pcribo.min)
pcribo.max <- if (!is.null(args$options$pcribo.max) && (args$options$pcribo.max != "NULL")) as.numeric(args$options$pcribo.max)
min.features <- if (!is.null(args$options$min.features) && (args$options$min.features != "NULL")) as.numeric(args$options$min.features)
min.counts <- if (!is.null(args$options$min.counts) && (args$options$min.counts != "NULL")) as.numeric(args$options$min.counts)
# QC gene
min.cells <- if (!is.null(args$options$min.cells) && (args$options$min.cells != "NULL")) as.numeric(args$options$min.cells)
### Databases
## Gene lists loading
mt.genes.file <- args$options$mt.genes.file
crb.genes.file <- args$options$crb.genes.file
str.genes.file <- args$options$str.genes.file
# Translation into gene Symbols
translation.file <- args$options$translation.file

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
    if(i %in% c("emptydrops.fdr", "droplets.limit", "min.features", "min.counts", "features.n")) assign(i, as.numeric(yaml_options[[i]])) else assign(i, yaml_options[[i]])
  }
  rm(yaml_options, i)
}
### Clean
rm(option_list,parser,args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Get Missing Paramaters ####
### Project
if (is.null(species)) species <- "homo_sapiens"
### Computational Parameters
if (is.null(nthreads)) nthreads <- 1
### Analysis Parameters
# Emptydrops
if (is.null(emptydrops.fdr)) emptydrops.fdr <- 1E-03
if (is.null(droplets.limit)) droplets.limit <- 1E+05
# Translate ENSG into Gene Symbol
if (is.null(translation)) translation <- TRUE
# QC cell
if (is.null(pcmito.min)) pcmito.min <- 0
if (is.null(pcmito.max)) pcmito.max <- 0.2
pcmito.range <- c(pcmito.min, pcmito.max)
if (is.null(pcribo.min)) pcribo.min <- 0
if (is.null(pcribo.max)) pcribo.max <- 1
pcribo.range <- c(pcribo.min, pcribo.max)
if (is.null(min.features)) min.features <- 200
if (is.null(min.counts)) min.counts <- 1000
# QC gene
if (is.null(min.cells)) min.cells <- 5
### Databases
if (species == "homo_sapiens") {
  # QC
  if (is.null(mt.genes.file)) mt.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/homo_sapiens_mito_symbols_20191001.rds")
  if (is.null(crb.genes.file)) crb.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/homo_sapiens_cribo_symbols_20191015.rds")
  if (is.null(str.genes.file)) str.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/homo_sapiens_stress_symbols_20200224.rds")
  # Translation into gene Symbols
  if (is.null(translation.file)) translation.file <- paste0(pipeline.path, "/resources/GENE_CONVERT/EnsemblToGeneSymbol_Homo_sapiens.GRCh38.txt")
}
if (species == "mus_musculus") {
  # QC
  if (is.null(mt.genes.file)) mt.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/mus_musculus_mito_symbols_20191015.rds")
  if (is.null(crb.genes.file)) crb.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/mus_musculus_cribo_symbols_20191015.rds")
  if (is.null(str.genes.file)) str.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/mus_musculus_stress_symbols_20200224.rds")
  # Translation into gene Symbols
  if (is.null(translation.file)) translation.file <- paste0(pipeline.path, "/resources/GENE_CONVERT/EnsemblToGeneSymbol_Mus_musculus.GRCm38.txt")
}

if (species == "rattus_norvegicus") {
  # QC
  if (is.null(mt.genes.file)) mt.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/rattus_norvegicus_mito_symbols_20200315.rds")
  if (is.null(crb.genes.file)) crb.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/rattus_norvegicus_cribo_symbols_20200315.rds")
  if (is.null(str.genes.file)) str.genes.file <- paste0(pipeline.path, "/resources/GENELISTS/rattus_norvegicus_stress_symbols_20200315.rds")
  # Translation into gene Symbols
  if (is.null(translation.file)) translation.file <- paste0(pipeline.path, "/resources/GENE_CONVERT/EnsemblToGeneSymbol_Rattus_norvegicus.Rnor_6.0.txt")
}

#### Fixed parameters ####
### Computational Parameters
my.seed <- 1337L
### Analysis Parameters
assay <- 'RNA'

#### Check non-optional parameters ####
if (is.null(input.dir.ge)) stop("input.dir.ge parameter can't be empty!")
if (is.null(output.dir.ge)) stop("output.dir.ge parameter can't be empty!")
if (is.null(sample.name.GE)) stop("sample parameter can't be empty!")

#### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

#########
## MAIN
#########

## A. QC_droplets ANALYSIS
##-------------------

### printing parameters:
print("###########################################")
print(paste0("sample.name.GE : ", sample.name.GE))
print(paste0("input.dir.ge : ", input.dir.ge))
print(paste0("output.dir.ge : ", output.dir.ge))
print("###########################################")

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Building output directory
unfiltred.dir <- paste0(output.dir.ge, 'QC_droplets', if(!is.null(emptydrops.retain)) '_retain', emptydrops.retain, '/')
dir.create(path = unfiltred.dir, recursive = TRUE, showWarnings = FALSE)

### Loading raw count matrix + Filtering duplicated cell barcodes + Removing empty droplets
sobj <- load.sc.data(data.path = input.dir.ge, sample.name = sample.name.GE, assay = assay, droplets.limit = droplets.limit, emptydrops.fdr = emptydrops.fdr, emptydrops.retain = emptydrops.retain, translation = translation, translation.file = translation.file, BPPARAM = cl, my.seed = my.seed, out.dir = unfiltred.dir)

### Save project parameters
sobj@misc$params$sample.name.GE <- sample.name.GE
sobj@misc$params$species <- species
sobj@misc$params$author.name <- author.name
sobj@misc$params$author.mail <- author.mail
sobj@misc$params$analysis_type <- "Individual analysis"

### Computing basic metrics : percentage of counts in the top features + mito + ribo + stress + nb features + nb counts
sobj <- QC.metrics(sobj = sobj, assay = assay, mt.genes.file = mt.genes.file, crb.genes.file = crb.genes.file, str.genes.file = str.genes.file, pcmito.range = pcmito.range, pcribo.range = pcribo.range, min.features = min.features, min.counts = min.counts, BPPARAM = cl)

### QC histograms
QC.hist(sobj = sobj, assay = assay, out.dir = unfiltred.dir)

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$part1_Droplets_QC <- paste0("Cell barcode by symbol count table were loaded in R (version ", getRversion(), ") using the BUSpaRse package (version ", sobj@misc$technical_info$BUSpaRse,").
To call real cells from empty droplets, we used the emptyDrops() function from the dropletUtils (version ", sobj@misc$technical_info$DropletUtils,") package, which assesses whether the RNA content associated with a cell barcode is significantly distinct from the ambient background RNA present within each sample. Barcodes with p-value < ", emptydrops.fdr," (Benjamini-Hochberg-corrected) were considered as legitimate cells for further analysis.")
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)

### Saving non-normalized object
save(sobj, file = paste0(unfiltred.dir, sample.name.GE, '_QC_NON-NORMALIZED.rda'), compress = "bzip2")
