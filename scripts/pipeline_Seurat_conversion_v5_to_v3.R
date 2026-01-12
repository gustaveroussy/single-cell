#### Read parameters ####
library(optparse)
option_list <- list(
  ### Project
  make_option("--input.rda.ge", help="Input seurat object (in .rda format) to convert in cerebro object."),
  ### Computational Parameters
  make_option("--pipeline.path", help="Path to pipeline folder")
)
parser <- OptionParser(usage="Rscript %prog [options]", description = " ", option_list = option_list)
args <- parse_args(parser, positional_arguments = 0)


### Sourcing functions ####
if(is.null(args$options$pipeline.path)) stop("--pipeline.path parameter must be set!")
source(paste0(args$options$pipeline.path, "/scripts/bustools2seurat_preproc_functions.R"))


#### Get Paramaters ####
### Formatting and extract from args
args <- convert_NULL_BOOL(args)
### Clean
rm(args)

#### Check non-optional parameters ####
if (is.null(input.rda.ge)) stop("input.rda.ge parameter can't be empty!")

#########
## MAIN
#########

### printing parameters:
print("###########################################")
print(paste0("input.rda.ge : ", input.rda.ge))
print("###########################################")

### Load data
load(input.rda.ge)
defaultAssayRNA_bool <- FALSE
if(Seurat::DefaultAssay(sobj) == "RNA") defaultAssayRNA_bool <- TRUE

### Change the name of the RNA V5 assay (RNA to RNA_V5)
sobj[["RNA_V5"]] <- sobj[["RNA"]]
Seurat::DefaultAssay(sobj) <-"RNA_V5"
sobj[["RNA"]] <- NULL

### Recreate a new RNA assay in version 3, from the version 5
sobj[["RNA"]] <- Seurat::CreateAssayObject(counts = Seurat::GetAssayData(sobj, assay="RNA_V5", slot="counts"),
                                   min.cells = 0,
                                   min.features = 0,
                                   key = 'rna_')
sobj[["RNA"]]@data <- Seurat::GetAssayData(sobj, assay="RNA_V5", layer="data")

### Remove RNA V5 assay
if(defaultAssayRNA_bool){
  Seurat::DefaultAssay(sobj) <-"RNA"
}else{
  Seurat::DefaultAssay(sobj) <-"SCT"
}
sobj[["RNA_V5"]] <- NULL

### Saving seurat v3 object
GE_file <- sub("\\.rda$", "", input.rda.ge)
save(sobj, file = paste0(GE_file, '_SeuratV3.rda'), compress = "bzip2")
