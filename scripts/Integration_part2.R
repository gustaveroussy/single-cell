## INTEGRATION PROTOCOL (Seurat or scbfa)
## To Do: test Harmony and LIGER integration

#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.int", help="Input filtred, normalized and reducted seurat object (in .rda format)."),
  make_option("--output.dir.int", help="Output path"),
  make_option("--markfile", help="Genes to plot on umap (# )format: 2 columns named Genes and Signature"),
  make_option("--author.name", help="Name of auhtor of the analysis"),
  make_option("--author.mail", help="Email of auhtor of the analysis"),
  ### Computational Parameters
  make_option("--nthreads", help="Number of threads to use"),
  make_option("--pipeline.path", help="Path to pipeline folder; it allows to change path if this script is used by snakemake and singularity, or singularity only or in local way. Example for singularity only: /WORKDIR/scRNAseq_10X_R4"),
  ### Analysis Parameters
  # Metadata
  make_option("--metadata.file", help="csv file with the metadata to add in the seurat objects"),
  # Clustering
  make_option("--keep.dims", help="Number of dimension to keep for clustering (from 0 to keep.dims)"),
  make_option("--keep.res", help="Resolution value for clustering"),
  # Annotation
  make_option("--custom.sce.ref", help="List of .RData files containing SingleCellExpriment objects with your reference"),
  make_option("--custom.markers.ref", help="List of .xlsx files containing your reference"),
  make_option("--cfr.minscore", help="Minimum correlation score for clustifyr to consider"),
  make_option("--sr.minscore", help="Minimum correlation score for SingleR to consider"),
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
input.rda.int <- args$options$input.rda.int
output.dir.int <- args$options$output.dir.int
markfile <- if (!is.null(args$options$markfile)) unlist(stringr::str_split(args$options$markfile, ","))
list.author.name <- if (!is.null(args$options$author.name)) unlist(stringr::str_split(args$options$author.name, ","))
list.author.mail <- if (!is.null(args$options$author.mail)) unlist(stringr::str_split(args$options$author.mail, ","))
### Computational Parameters
nthreads <-  if (!is.null(args$options$nthreads)) as.numeric(args$options$nthreads)
pipeline.path <- args$options$pipeline.path
### Analysis Parameters
# Metadata
metadata.file <- unlist(stringr::str_split(args$options$metadata.file, ","))
# Clustering
keep.dims <- if (!is.null(args$options$keep.dims)) as.numeric(args$options$keep.dims)
keep.res <- if (!is.null(args$options$keep.res)) as.numeric(args$options$keep.res)
# Annotation
custom.sce.ref <- if (!is.null(args$options$custom.sce.ref)) unlist(stringr::str_split(args$options$custom.sce.ref, ","))
custom.markers.ref <- if (!is.null(args$options$custom.markers.ref)) unlist(stringr::str_split(args$options$custom.markers.ref, ","))
cfr.minscore <- if (!is.null(args$options$cfr.minscore)) as.numeric(args$options$cfr.minscore)
sr.minscore <- if (!is.null(args$options$sr.minscore)) as.numeric(args$options$sr.minscore)
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
    assign(i, yaml_options[[i]])
    if(i %in% c("nthreads","keep.dims","keep.res","cfr.minscore","sr.minscore")) assign(i, as.numeric(yaml_options[[i]]))else assign(i, yaml_options[[i]])
    
  }
  rm(yaml_options, i)
}
### Clean
rm(option_list,parser,args)

#### Get path if snakemake/singularity/local ####
if(is.null(pipeline.path)) stop("--pipeline.path parameter must be set!")

#### Check non-optional parameters ####
if (is.null(input.rda.int)) stop("input.rda.int parameter can't be empty!")
if (is.null(output.dir.int)) stop("output.dir.int parameter can't be empty!")
if (is.null(keep.dims)) stop("keep.dims parameter can't be empty!")
if (is.null(keep.res)) stop("keep.res parameter can't be empty!")

### Load data
load(input.rda.int)

### Sourcing functions ####
source(paste0(pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))

### Save project parameters
sobj <- Add_name_mail_author(sobj = sobj, list.author.name = list.author.name, list.author.mail = list.author.mail)

#### Get Missing Paramaters ####
### Project
name.int <- sobj@misc$params$name.int
species <- sobj@misc$params$species
### Computational Parameters
if (is.null(nthreads)) nthreads <- 4
### Analysis Parameters
# Integration
integration.method <- sobj@misc$params$integration$method
# Normalization and dimension reduction
assay <- sobj@misc$params$integration$out.assay
norm.method <- sobj@assays[[assay]]@misc$params$normalization$normalization.method
dimred.method <- sobj@assays[[assay]]@misc$params$reductions$method
if(integration.method == "Seurat") red.name <- paste(c("integrated", dimred.method, integration.method), collapse = '_')
if(integration.method %in% c('scbfa', 'bpca','Liger')) red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
if(integration.method == "Harmony") red.name <- paste(c(assay, dimred.method), collapse = '_')
if(integration.method == "Harmony") red.name <- paste(c(assay, dimred.method, integration.method), collapse = '_')
# Annotation
if (is.null(cfr.minscore)) cfr.minscore <- 0.35
if (is.null(sr.minscore)) sr.minscore <- 0.25

#### Fixed parameters ####
# Annotation
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
if (species == "rattus_norvegicus") {
  singler.setnames <- c("MouseRNAseqData", "ImmGenData")
  clustifyr.setnames <- c("ref_MCA", "ref_tabula_muris_drop", "ref_tabula_muris_facs", "ref_moca_main", "ref_immgen", "ref_mouse.rnaseq")
  scrnaseq.setnames <- c("BaronPancreasData(mouse)","ZeiselBrainData") #,"ShekharRetinaData"
}

#### Fixed parameters ####
# Plots
solo.pt.size <- 3
multi.pt.size <- 2
gradient.cols <- c("gold", "blue")

#### Get genes markers ####
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

## RUN
######

print("#####################################")
print(paste0("Sample: ", name.int))
print(paste0("RDA file: ", input.rda.int))
print(paste0("Dimension: ", keep.dims))
print(paste0("Resolution: ", keep.res))
print("#####################################")

### Creating parallel instance
cl <- create.parallel.instance(nthreads = nthreads)

### Add metadata
if(!is.null(metadata.file)) sobj <- add_metadata_sobj(sobj=sobj, metadata.file = metadata.file)

### Building clustered output directory
clust.dir <- paste(output.dir.int, paste0("dims", keep.dims, "_res", keep.res), sep = '/')
dir.create(path = clust.dir, recursive = TRUE, showWarnings = FALSE)

### Replotting final clusters
cat("\nClustering...\n")
sobj <- louvain.cluster(sobj = sobj, reduction = red.name, max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size, algorithm = 1)

## Setting ident name and RNA.reduction
ident.name <- paste0(paste0(red.name, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))
INT.reduction <- paste(c(red.name, keep.dims, 'umap'), collapse = '_')
#sobj@reductions[[red.name]]@misc$from.assay <- assay

### uMAP plot by sample
cat("\nuMAP plot by sample...\n")
blockpix = 600
png(filename = paste0(clust.dir, '/', paste(c(name.int, red.name, 'uMAP.png'), collapse = '_')), width = 1000, height = 1000)
print(Seurat::DimPlot(object = sobj, reduction = paste(c(red.name, keep.dims, 'umap'), collapse = '_'), order = sample(x = 1:ncol(sobj), size = ncol(sobj), replace = FALSE), group.by = 'orig.ident', pt.size = solo.pt.size) + ggplot2::ggtitle("uMAP for all samples ") + Seurat::DarkTheme())
dev.off()
grid.xy <- grid.scalers(length(unique(sobj@meta.data$orig.ident)))
png(filename = paste0(clust.dir, '/', paste(c(name.int, red.name, 'split', 'uMAP.png'), collapse = '_')), width = grid.xy[1]*blockpix, height = grid.xy[2]*blockpix)
print(Seurat::DimPlot(object = sobj, reduction = paste(c(red.name, keep.dims, 'umap'), collapse = '_'), group.by = ident.name, split.by = 'orig.ident', pt.size = solo.pt.size, ncol = grid.xy[1]) + ggplot2::ggtitle(paste0("uMAP split on samples")) + Seurat::DarkTheme())
dev.off()

### Technical plots
cat("\nSaving technical plots...\n")
technical.plot(sobj = sobj, ident = ident.name, out.dir = clust.dir, multi.pt.size = multi.pt.size)

### Finding markers
cat("\nFinding markers...\n")
sobj <- find.markers.quick(sobj = sobj, ident = ident.name, test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, adjp.p.max = 5E-02, topn = 10, heatmap.cols = gradient.cols, out.dir = clust.dir)

### Automatic cell type annotation
cat("\nAutomatic cell type annotation...\n")
sobj <- cells.annot(sobj = sobj, ident = ident.name, singler.setnames = singler.setnames, clustifyr.setnames = clustifyr.setnames, scrnaseq.setnames = scrnaseq.setnames, custom.sce.ref = custom.sce.ref, custom.markers.ref = custom.markers.ref, sr.minscore = sr.minscore, cfr.minscore = cfr.minscore, out.dir = clust.dir, solo.pt.size = solo.pt.size, nthreads = nthreads)

### Assessing clusters : Plotting provided marker genes
cat("\nPlotting provided marker genes...\n")
if(!is.null(markers)) sobj <- markers.umap.plot(sobj = sobj, markers = markers, ident = ident.name, out.dir = clust.dir, dimplot.cols = gradient.cols, multi.pt.size = 2)

### Materials and Methods
sobj@misc$parameters$Materials_and_Methods$Integration_Clust_Markers_Annot <- paste0(
  "An automatic annotation of cell types was perfom by SingleR (version ",sobj@misc$technical_info$SingleR,") (with fine-tuning step) and ClustifyR (version ",sobj@misc$technical_info$clustifyr,"), using packages built-in references. It labels clusters (or cells) from a dataset based on similarity (Spearman correlation score) to a reference dataset with known labels. The labels with a correlation score greater than ",sr.minscore," for SingleR or greater than ",cfr.minscore," for ClustifyR were kept. The annotation was also made CelliD (version ",sobj@misc$technical_info$CelliD,") with genes signatures from pangloa database. ",
  "Marker genes for Louvain clusters and samples, were identified through a «one versus others» differential analysis using the Wilcoxon test through the FindAllMarkers() function from Seurat, considering only genes with a minimum log fold-change of 0.5 in at least 75% of cells from one of the groups compared, and FDR-adjusted p-values <0.05 (Benjaminin-Hochberg method)."
)
sobj@misc$parameters$Materials_and_Methods$References_packages <- find_ref(MandM = sobj@misc$parameters$Materials_and_Methods, pipeline.path = pipeline.path)
write_MandM(sobj=sobj, output.dir=clust.dir)

### Saving final object
cat("\nSaving object...\n")
GE_file=paste0(clust.dir, '/', paste(c(sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(input.rda.int)), keep.dims, keep.res), collapse = "_"))
save(sobj, file = paste0(GE_file, '.rda'), compress = "bzip2")
