#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.ge", help="Input seurat object (in .rda format) to convert in cerebro object."),
  make_option("--author.name", help="Name of author of the analysis"),
  make_option("--author.mail", help="Email of author of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  # Cerebro
  make_option("--version", help="Version of cerebro to use (v1.2 or v1.3 (default))."),
  make_option("--groups", help="Column name (in meta.data) to define clusters/comparisons for Cerebro object (usefull for TCR/BCR part). The last RNA clustering and samples information are already included."),
  make_option("--remove.other.reductions", help="Remove all other reductions present in seurat object (keep only final umap)"),
  make_option("--remove.other.idents", help="Remove all other clustering present in seurat object (keep only the last clustering)"),
  make_option("--remove.mt.genes", help="Remove mitochondrial genes (to see better the other genes)"),
  make_option("--remove.crb.genes", help="Remove ribosomal genes (to see better the other genes)"),
  make_option("--remove.str.genes", help="Remove stress genes (to see better the other genes)"),
  make_option("--only.pos.DE", help="Keep only positive DE genes from customized differential expression analysis (for genes markers identification is always only positive)."),
  make_option("--remove.custom.DE", help="Remove results from customized differential expression analysis."),
  ### Databases
  # Metadata
  make_option("--metadata.file", help="csv file with the metadata to add in the seurat object"),
  # Cerebro
  make_option("--gmt.file", help="GMT file for cerebro"),
  ### Yaml parameters file to remplace all parameters before (usefull to use R script without snakemake)
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
input.rda.ge <- args$options$input.rda.ge
list.author.name <- if (!is.null(args$options$author.name)) unlist(stringr::str_split(args$options$author.name, ","))
list.author.mail <- if (!is.null(args$options$author.mail)) unlist(stringr::str_split(args$options$author.mail, ","))
### Computational Parameters
nthreads <- args$options$nthreads
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Cerebro
version <- args$options$version
groups <- args$options$groups
remove.other.reductions <- args$options$remove.other.reductions
remove.other.idents <- args$options$remove.other.idents
remove.mt.genes <- args$options$remove.mt.genes
remove.crb.genes <- args$options$remove.crb.genes
remove.str.genes <- args$options$remove.str.genes
only.pos.DE <- args$options$only.pos.DE
remove.custom.DE <- args$options$remove.custom.DE
### Databases
# Metadata
metadata.file <-  if (!is.null(args$options$metadata.file)) unlist(stringr::str_split(args$options$metadata.file, ","))
# Cerebro
gmt.file <- args$options$gmt.file
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
    if(i %in% c("nthreads")) assign(i, as.numeric(yaml_options[[i]])) else assign(i, yaml_options[[i]])
  }
  rm(yaml_options, i)
}
### Clean
rm(option_list,parser,args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(input.rda.ge)) stop("input.rda.ge parameter can't be empty!")

### Load data
load(input.rda.ge)

### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

### Save project parameters
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)

#### Get Missing Paramaters ####
### Project
species <- sobj@misc$params$species
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
# Normalization and dimension reduction
dimred.method <- sobj@misc$params$reductions$method
## Clustering
GE_file <- sub("\\.rda$", "", input.rda.ge)
keep.dims <- sobj@misc$params$clustering$max.dim
keep.res <- sobj@misc$params$clustering$resolution
ident.name <- sobj@misc$params$clustering$ident
RNA.reduction <- sobj@misc$params$clustering$umap
# Cerebro
if (is.null(version)) version <- "v1.3"
if (is.null(remove.other.reductions)) remove.other.reductions <- FALSE
if (is.null(remove.other.idents)) remove.other.idents <- FALSE
if (is.null(remove.mt.genes)) remove.mt.genes <- FALSE
if (is.null(remove.crb.genes)) remove.crb.genes <- FALSE
if (is.null(remove.str.genes)) remove.str.genes <- FALSE
if (is.null(only.pos.DE)) only.pos.DE <- FALSE
if (is.null(remove.custom.DE)) remove.custom.DE <- FALSE
if (is.null(groups)){
  groups_v1.2 <- NULL
}else{
  groups_v1.2 <- unlist(stringr::str_split(groups, ","))
}
### Databases
# Cerebro
if (is.null(gmt.file)) gmt.file <- paste0(pipeline.path, "/resources/DATABASE/MSIGDB/7.1/msigdb_v7.1_GMTs/msigdb.v7.1.symbols.gmt")

#### Fixed parameters ####
# Cerebro
if (species == "homo_sapiens") species.rdx <- 'hg'
if (species == "mus_musculus") species.rdx <- 'rn'
if (species == "rattus_norvegicus") species.rdx <- 'rn'

#########
## MAIN
#########

### printing parameters:
print("###########################################")
print(paste0("input.rda.ge : ", input.rda.ge))
print(paste0("clustering : ", ident.name))
print("###########################################")

### Add metadata
if(!is.null(metadata.file)) sobj <- add_metadata_sobj(sobj=sobj, metadata.file = metadata.file)

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$Cerebro <- paste0("Cerebro (version ",utils::packageVersion('cerebroApp'),"), cell report browser, is an AppShiny which allows users to interactively visualize various parts of single cell transcriptomics analysis without requiring bioinformatics expertise. This package is also used to identify most expressed genes, and to compute pathway enrichment on marker genes (based on the Enrichr API) and Gene Set Enrichment Analysis (uses the Gene Set Variation Analysis method in combination with additional statistics as published by Diaz-Mejia et. al.(Diaz-Mejia JJ, Meng EC, Pico AR et al. Evaluation of methods to assign cell type labels to cell clusters from single-cell RNA-sequencing data [version 3; peer review: 2 approved, 1 approved with reservations]. F1000Research 2019, 8(ISCB Comm J):296 (https://doi.org/10.12688/f1000research.18490.3))) from the MSigDB (H collection: hallmark gene sets).")
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj=sobj, output.dir=dirname(input.rda.ge))

### Building cerebro binary
cat("\nBuilding cerebro object...\n")
if (version == "v1.2"){
    seurat2cerebro(sobj = sobj, ident = ident.name, clusters.colnames = NULL, remove.other.reductions = remove.other.reductions, remove.other.idents = remove.other.idents, species = species.rdx, gmt.file = gmt.file, remove.mt.genes = remove.mt.genes, remove.crb.genes = remove.crb.genes, remove.str.genes = remove.str.genes, file = GE_file, nthreads = nthreads, only_pos = TRUE, only_pos_DE = only.pos.DE, remove.custom.DE = remove.custom.DE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 5E-02, test = 'wilcox')
  if(!is.null(groups_v1.2)){
    for (clusters.colnames in groups_v1.2){
      seurat2cerebro(sobj = sobj, ident = ident.name, clusters.colnames = clusters.colnames, remove.other.reductions = remove.other.reductions, remove.other.idents = remove.other.idents, species = species.rdx, gmt.file = gmt.file, remove.mt.genes = remove.mt.genes, remove.crb.genes = remove.crb.genes, remove.str.genes = remove.str.genes, file = GE_file, nthreads = nthreads, only_pos = TRUE, only_pos_DE = only.pos.DE, remove.custom.DE = remove.custom.DE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 5E-02, test = 'wilcox')
    }
  }
}else if (version == "v1.3"){
  seurat2cerebro_1.3(sobj = sobj, ident = ident.name, groups = groups, species = species.rdx, remove.other.reductions = remove.other.reductions, remove.other.idents = remove.other.idents, gmt.file = gmt.file, remove.mt.genes = remove.mt.genes, remove.crb.genes = remove.crb.genes, remove.str.genes = remove.str.genes, file = GE_file, nthreads = nthreads, min_pct = .75, thresh_logFC = .5, thresh_p_val = 5E-02, test = 'wilcox', only_pos = TRUE, remove.custom.DE = remove.custom.DE, only_pos_DE = only.pos.DE)
}