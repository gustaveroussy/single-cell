#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.ge", help="Input unfiltred seurat object (in .rda format)."),
  make_option("--output.dir.ge", help="Output path"),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
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
  make_option("--doublets.filter.method", help="Method used to filter doublets (scDblFinder, scds, both). To not filter set the parameter to none."),
  ### Databases
  # QC
  make_option("--cc.seurat.file", help="RDS file with list of cell cycle genes"),
  make_option("--cc.cyclone.file", help="RDS file with list of cell cycle genes"),
  ### Yaml parameters file to remplace all parameters before (usefull to use R script without snakemake)
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
input.rda.ge <- args$options$input.rda.ge
output.dir.ge <- args$options$output.dir.ge
author.name <- args$options$author.name
author.mail <- args$options$author.mail
### Computational Parameters
nthreads <- if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# QC cell
pcmito.min <- if (!is.null(args$options$pcmito.min)) as.numeric(args$options$pcmito.min)
pcmito.max <- if (!is.null(args$options$pcmito.max)) as.numeric(args$options$pcmito.max)
pcribo.min <- if (!is.null(args$options$pcribo.min)) as.numeric(args$options$pcribo.min)
pcribo.max <- if (!is.null(args$options$pcribo.max)) as.numeric(args$options$pcribo.max)
min.features <- if (!is.null(args$options$min.features)) as.numeric(args$options$min.features)
min.counts <- if (!is.null(args$options$min.counts)) as.numeric(args$options$min.counts)
# QC gene
min.cells <- if (!is.null(args$options$min.cells)) as.numeric(args$options$min.cells)
# Doublets
doublets.filter.method <- args$options$doublets.filter.method
### Databases
## Gene lists loading
cc.seurat.file <- args$options$cc.seurat.file
cc.cyclone.file <- args$options$cc.cyclone.file
### Yaml parameters file to remplace all parameters before (usefull to use R script without snakemake)
if (!is.null(args$options$yaml)){
  yaml_options <- yaml::yaml.load_file(args$options$yaml)
  for(i in names(yaml_options)) if(i %in% c("min.features", "min.counts", "features.n")) assign(i, as.numeric(yaml_options[[i]]))else assign(i, yaml_options[[i]])
  rm(yaml_options)
}
### Clean
rm(option_list,parser,args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(input.rda.ge)) stop("input.rda.ge parameter can't be empty!")
if (is.null(output.dir.ge)) stop("output.dir.ge parameter can't be empty!")

### Load data
load(input.rda.ge)

### Save project parameters
if (!is.null(author.name) && !tolower(author.name) %in% tolower(sobj@misc$params$author.name)) sobj@misc$params$author.name <- c(sobj@misc$params$author.name, author.name)
if (!is.null(author.mail) && !tolower(author.mail) %in% tolower(sobj@misc$params$author.mail)) sobj@misc$params$author.mail <- c(sobj@misc$params$author.mail, author.mail)

#### Get Missing Paramaters ####
### Project
sample.name.GE <- sobj@misc$params$sample.name.GE
species <- sobj@misc$params$species
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
# QC cell
if (is.null(pcmito.min)) pcmito.min <- 0
if (is.null(pcmito.min)) pcmito.min <- 0.2
pcmito.range <- c(pcmito.min, pcmito.max)
if (is.null(pcribo.min)) pcribo.min <- 0
if (is.null(pcribo.max)) pcribo.max <- 1
pcribo.range <- c(pcribo.min, pcribo.max)
if (is.null(min.features)) min.features <- 200
if (is.null(min.counts)) min.counts <- 1000
# QC gene
if (is.null(min.cells)) min.cells <- 5
# Doublets
if (is.null(doublets.filter.method)) doublets.filter.method <- 'all'
### Databases
if (species == "homo_sapiens") {
  # QC
  if (is.null(cc.seurat.file)) cc.seurat.file <-  paste0(pipeline.path,"/resources/GENELISTS/homo_sapiens_Seurat_cc.genes_20191031.rds")
  if (is.null(cc.cyclone.file)) cc.cyclone.file <- paste0(pipeline.path,"/resources/GENELISTS/homo_sapiens_cyclone_pairs_symbols_20191001.rds")
}
if (species == "mus_musculus") {
  # QC
  if (is.null(cc.seurat.file)) cc.seurat.file <- paste0(pipeline.path,"/resources/GENELISTS/mus_musculus_Seurat_cc.genes_20191031.rds")
  if (is.null(cc.cyclone.file)) cc.cyclone.file <- paste0(pipeline.path,"/resources/GENELISTS/mus_musculus_cyclone_pairs_symbols_20191015.rds")
}
if (species == "rattus_norvegicus") {
  # QC
  if (is.null(cc.seurat.file)) cc.seurat.file <- paste0(pipeline.path,"/resources/GENELISTS/rattus_norvegicus_Seurat_cc.genes_20200315.rds")
  if (is.null(cc.cyclone.file)) cc.cyclone.file <- paste0(pipeline.path,"/resources/GENELISTS/rattus_norvegicus_cyclone_pairs_symbols_20200315.rds")
}

#### Fixed parameters ####
### Analysis Parameters
assay <- 'RNA'

#### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

#########
## MAIN
#########

## B. HALF-FILTERED ANALYSIS + DOUBLETS ANALYSIS
##-------------------

### printing parameters:
print("###########################################")
print(paste0("sample : ",sample.name.GE))
print(paste0("input.rda.ge : ",input.rda.ge))
print(paste0("output.dir.ge : ",output.dir.ge))
print("###########################################")

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Filtering cells
cat("\nFiltering cells...\n")
sobj <- cells.QC.filter(sobj = sobj, min.features = min.features, min.counts = min.counts, pcmito.range = pcmito.range, pcribo.range = pcribo.range)

### Cell cycle prediction
cat("\nCell cycle prediction...\n")
sobj <- cell.cycle.predict(sobj = sobj, assay = assay, cc.cyclone.file = cc.cyclone.file, cc.seurat.file = cc.seurat.file, nbin = 10, BPPARAM = cl)

### Filtering features (based on minimum cells covering)
cat("\nFiltering features...\n")
sobj <- features.filter(sobj = sobj, min.cells = min.cells)

### Identification of doublets
cat("\nIdentification of doublets...\n")
sobj <- find.doublets(sobj = sobj, assay = assay)

### Building output directory
filtered.dir <- paste0(output.dir.ge, paste0('F', min.features, '_C', min.counts, '_M', paste(pcmito.range, collapse = '-'), '_R', paste(pcribo.range, collapse = '-'),'_G', min.cells))


### C. DOUBLETS KEPT:
##-------------------

### Building doublet kept output directory
doublet.kept.dir <- paste0(filtered.dir, '/DOUBLETSKEPT/')
dir.create(path = doublet.kept.dir, recursive = TRUE, showWarnings = FALSE)

### Computing basic metrics : percentage of counts in the top features + mito + ribo + stress + nb features + nb counts
cat("\nComputing QC metrics with kept doublets...\n")
sobj <- QC.metrics(sobj = sobj, assay = assay, pcmito.range = pcmito.range, pcribo.range = pcribo.range, min.features = min.features, min.counts = min.counts, BPPARAM = cl)

### QC histograms
QC.hist(sobj = sobj, assay = assay, out.dir = doublet.kept.dir)

### Saving stat for sheet
if(doublets.filter.method == 'none'){
  message("Saving stat...\n")
  sobj@misc$excel$Final_Cells_Quality$nb_genes <- dim(sobj)[1]
  sobj@misc$excel$Final_Cells_Quality$nb_cells <- dim(sobj)[2]
  save_stat(sobj = sobj, sample.name = sample.name.GE, title = sample.name.GE, out.dir = doublet.kept.dir)
  ### Materials and Methods
  sobj@misc$Materials_and_Methods$part2_Filtering_droplets <- paste0("The count matrix was filtered to exclude genes detected in less than ",min.cells," cells, cells with less than ",min.counts," UMIs or less than ",min.features," detected genes, as well as cells with mitochondrial transcripts proportion higher than ",pcmito.range[2]*100,"%",
                                                              if(!pcmito.range[1]==0) paste0("and less than ",pcmito.range[1]*100, "%"),".",
                                                              if(!pcribo.range[2]==1 && !pcribo.range[1]==0) paste0("as well as cells with ribosomal transcripts proportion higher than ",pcribo.range[2]*100,"% and less than ",pcribo.range[1]*100, "%. The proportion of mechanical stress-response gene counts (Thesis of Léo Machado entitled « From skeletal muscle stem cells to tissue atlases: new tools to investigate and circumvent dissociation-induced stress », 2019) were also estimated but not used to filter cells.") else "The proportion of ribosomal gene counts and the proportion of mechanical stress-response gene counts (Thesis of Léo Machado entitled « From skeletal muscle stem cells to tissue atlases: new tools to investigate and circumvent dissociation-induced stress », 2019) were also estimated but not used to filter cells.",
"Cell cycle scoring of each cell was performed using two methods : the CellcycleScoring() function from the Seurat package (version ",sobj@misc$technical_info$Seurat,"), and the cyclone() function from Scran (version ",sobj@misc$technical_info$scran,").

Barcodes corresponding to doublet cells were identified using the union of two methods: scDblFinder (version ",sobj@misc$technical_info$scDblFinder,") using default parameters, and scds (version ",sobj@misc$technical_info$scds,") with its hybrid method using default parameters. However boublets were not discarded.")
  sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
}

### Saving doublets kept non-normalized object
cat("\nSaving object...\n")
save(sobj, file = paste0(doublet.kept.dir, '/', sample.name.GE, '_DOUBLETSKEPT_NON-NORMALIZED.rda'), compress = "bzip2")

### Basic normalization, dimension reduction, clustering, technical plots and save
cat("\nComputing complete analysis with doublets kept...\n")
norm.red.plot.quick(sobj = sobj, sample.name.GE = sample.name.GE, pre.out.dir = doublet.kept.dir, file.name = 'DOUBLETSKEPT')


### D. DOUBLETS REMOVED:
##-------------------

if(doublets.filter.method != 'none'){

  ### Building doublet filter output directory
  doublet.filtered.dir <-  paste0(filtered.dir, '/DOUBLETSFILTER_', doublets.filter.method, '/')
  dir.create(path = doublet.filtered.dir, recursive = TRUE, showWarnings = FALSE)
  
  ### Filter doublets
  cat("\nFiltering doublets...\n")
  sobj <- filter.doublets(sobj = sobj, method = doublets.filter.method)
  
  ### Computing basic metrics : percentage of counts in the top features + mito + ribo + stress + nb features + nb counts
  cat("\nComputing QC metrics with removed doublets...\n")
  sobj <- QC.metrics(sobj = sobj, assay = assay, pcmito.range = pcmito.range, pcribo.range = pcribo.range, min.features = min.features, min.counts = min.counts, BPPARAM = cl)

  ### QC histograms
  QC.hist(sobj = sobj, assay = assay, out.dir = doublet.filtered.dir)

  ### Saving stat for sheet
  cat("\nSaving stat...\n")
  sobj@misc$excel$Final_Cells_Quality$nb_genes <- dim(sobj)[1]
  sobj@misc$excel$Final_Cells_Quality$nb_cells <- dim(sobj)[2]
  save_stat(sobj = sobj, sample.name = sample.name.GE, title = sample.name.GE, out.dir = doublet.filtered.dir)

  ### Materials and Methods
  sobj@misc$parameters$Materials_and_Methods$part2_Filtering <- paste0("The count matrix was filtered to exclude genes detected in less than ",min.cells," cells, cells with less than ",min.counts," UMIs or less than ",min.features," detected genes, as well as cells with mitochondrial transcripts proportion higher than ",pcmito.range[2]*100,"%",
                                                              if(!pcmito.range[1]==0) paste0("and less than ",pcmito.range[1]*100, "%"),".",
                                                              if(!pcribo.range[2]==1 && !pcribo.range[1]==0) paste0("as well as cells with ribosomal transcripts proportion higher than ",pcribo.range[2]*100,"% and less than ",pcribo.range[1]*100, "%. The proportion of mechanical stress-response gene counts (Thesis of Léo Machado entitled « From skeletal muscle stem cells to tissue atlases: new tools to investigate and circumvent dissociation-induced stress», 2019) were also estimated but not used to filter cells.") else "The proportion of ribosomal gene counts and the proportion of mechanical stress-response gene counts (Thesis of Léo Machado entitled « From skeletal muscle stem cells to tissue atlases: new tools to investigate and circumvent dissociation-induced stress», 2019) were also estimated but not used to filter cells.",
                                                              "Cell cycle scoring of each cell was performed using two methods: the CellcycleScoring() function from the Seurat package (version ",sobj@misc$technical_info$Seurat,"), and the cyclone() function from Scran (version ",sobj@misc$technical_info$scran,").

Barcodes corresponding to doublet cells were identified and discarded using the union of two methods: scDblFinder (version ",sobj@misc$technical_info$scDblFinder,") using default parameters, and scds (version ",sobj@misc$technical_info$scds,") with its hybrid method using default parameters. We manualy verified that the cells identified as doublets did not systematically correspond to cells in G2M phase.")
  sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
  
  ### Saving non-normalized object
  cat("\nSaving object...\n")
  save(sobj, file = paste0(doublet.filtered.dir, sample.name.GE, '_FILTERED_NON-NORMALIZED.rda'), compress = "bzip2")
}
