## THIS IS THE MASTER SCRIPT FOR 10X SCRNASEQ PREPROCESSING FROM COUNT MATRIX
## TO A NORMALIZED (EVEN COVARIATES-REGRESSED) DATASET SAVED AS A RICH SEURAT
## OBJECT.
## ---
## This recommended workflow comes into 3 to 4 steps :
## 1) A completely unfiltered first step, just to display the first space of ALL
##    identified cells, despite putatively many bias source present in the data.
## 2) A filtered step, except for identified doublets : just to visualize them.
## 3) A filtered + doublets removed step, but with "basic" normalization (ie :
##    no covariate regressed)
## 4) Depending on the necessity, same as 3) but with one or multiple covariates
##    regressed.
## ---
## The ACTUAL, POOR script structure contains :
## 1) A series of blocks of variables, which are run/project -specific (a
##    detailed example is provided in the first blocl, "P30-MAMA1")
## 2) A block of common variables (variables that are very rarely to get
##    modified, paths to resources, etc...)
## 3) The functions, where all the magic happens.
## 4) A series of code blocks, corresponding to examples of each of the 4
##    described steps.
## ---
##
## ---
## WARNINGS :
## 1) PLEASE ONLY CONSIDER FUNCTIONS AS REAL CODE. CODE BLOCKS ARE LOSELY
##    MAINTAINED AND SHOULD BE BASIC HELPERS FOR THE ORDER OF FUNCTION CALLS.
## 2) The "clustering.eval.mt" function is a multithreaded version of the
##    "clustering.eval" function. Everything was done to limit its memory
##    hunger (a minimal Seurat object is temporarily created, only necessary
##    objects are exported), but running it on 4 threads can take up to 3x the
##    memory of the monothreaded version !
## 3) The use of scBFA/BPCA implies a variation of the pipeline and its steps.
##    Please read the corresponding function header.

## INSTALL :
## R3.6.x & rstudio installed via conda -c r
## packages installed from BioConductor via BiocManager::install :
## Seurat, scran, DropletUtils, scds, scDblFinder, scBFA, (celda)

## GUILTY AUTHOR : Bastien JOB (bastien.job@gustaveroussy.fr)




## FUNCTIONS
############

## multithreading cluster
create.parallel.instance <- function(nthreads = 1) {
  doParallel::registerDoParallel(nthreads)
  cl <- BiocParallel::DoparParam()
  return(cl)
}

## Loading data into a Seurat object
# 1) Loading data
# 2) Filtering duplicated cell barcodes
# 3) Rename ensembl genes id by genes symbols
# 4) Remove empty droplets
# 5) Plotting saturation and Kneeplot
# 6) Creation of the Seurat object
load.sc.data <- function(data.path = NULL, sample.name = NULL, assay = 'RNA', droplets.limit = 1E+05, emptydrops.fdr = 1E+03, emptydrops.retain = NULL, return.matrix = FALSE, translation = FALSE, translation.file = NULL, BPPARAM = BiocParallel::SerialParam(), my.seed = 1337, out.dir = NULL, draw_plots = TRUE) {
  if (file.exists(data.path) && !is.null(sample.name)) {
    message("Loading data ...")

    ## Loading data
    source.format <- ""
    if(file.exists(paste0(data.path, "/matrix.mtx"))) { ### Cell Ranger
      source.format <- "CellRanger"
      scmat <- Seurat::Read10X(data.path)
      if ('Gene Expression' %in% names (scmat)) {
        message("Keeping only gene expression")
        scmat <- scmat[['Gene Expression']]
      }
    } else if(file.exists(paste0(data.path, "/", sample.name, ".mtx"))) { ### BUStools
      source.format <- "BUStools"
      scmat <- BUSpaRse::read_count_output(dir = data.path, name = sample.name, tcc = FALSE)
    } else if (file.exists(paste0(data.path, "/quants_mat.gz"))) { ### Alevin
      source.format <- "Alevin"
      scmat <- Seurat::ReadAlevin(data.path)
    } else if (file.exists(paste0(data.path, "/", sample.name, "_counts.tsv.gz"))) { ### UMI-tools
      source.format <- "UMIt-ools"
      scmat <- read.table(file = paste0(data.path, "/", sample.name, "_counts.tsv.gz"), header = TRUE, sep = "\t", quote = "", check.names = FALSE, row.names = 1)
    } else {
      stop(paste0("No data found in [", data.path, "] (wrong path ?)"))
    }
    message(paste0("Found ", source.format, " data"))

    scmat <- scmat[,order(colnames(scmat))]

    message('Droplets matrix dimensions :')
    droplets.nb <- ncol(scmat)
    print(dim(scmat))

    message('Total UMIs :')
    umi.total.nb <- sum(scmat)
    print(umi.total.nb)

    ## Rename ensembl genes id by genes symbols
    if (translation){
      message('Translation to genes symbols...')
      data = read.table(file = translation.file, header = FALSE, sep = " ")
      gene_name=vector()
      for(i in 1:nrow(scmat)) {
        index = grep(gsub("\\.[0-9]*$", "",rownames(scmat)[i]), as.vector(data[,1]))
        if(!is.na(index)) gene_name[i] = as.vector(data[index,2]) else gene_name[i] = rownames(scmat)[i]
      }
      ##deduplicate lines
      #identify duplicate genes names and position
      dup.genes <- unique(gene_name[duplicated(gene_name)])
      if(length(dup.genes) > 0) {
        dup.pos=grep(paste0("^",paste(dup.genes,collapse="$|^"),"$"), gene_name)
        message(paste0("Found ", length(dup.genes), ' (', sprintf("%.2f", length(dup.genes) / nrow(scmat) * 100), "%) replicated genes causes by translation! Summing ..."))
        #data not duplicated
        scmat_without_dup = scmat[-dup.pos,]
        rownames(scmat_without_dup)=gene_name[-dup.pos]
        #data duplicated
        dup_scmat = as.data.frame(as.matrix(scmat[dup.pos,]))
        dup_gene_names=gene_name[dup.pos]
        rm(scmat)
        #transform in non duplicated data
        dedup_scmat = Matrix::Matrix(as.matrix(rowsum(dup_scmat,group=dup_gene_names)), sparse = TRUE)
        #merge data
        scmat = rbind(scmat_without_dup,dedup_scmat)
        rm(scmat_without_dup,dedup_scmat)
      }else{
        rownames(scmat) = gene_name
        message('No replicated gene found.')
      }
    }

    ## df for plot saturation and Kneeplot beging
    if(draw_plots){
      library(dplyr)
      nb_umi_genes_by_barcode <- data.frame(nb_genes=Matrix::colSums(scmat>0), nb_umi=Matrix::colSums(scmat), barcodes=colnames(scmat))
      nb_umi_genes_by_barcode <- nb_umi_genes_by_barcode %>% arrange(desc(nb_umi,nb_genes)) %>% dplyr::mutate(num_barcode=seq.int(ncol(scmat)))
    }
    if (!is.null(droplets.limit) && ncol(scmat) > droplets.limit && !is.null(emptydrops.fdr)) {
      ## Removing empty droplets
      message("Removing empty droplets with emptyDrops")
      bc_rank <- DropletUtils::barcodeRanks(scmat)
      set.seed(my.seed)
      bc_rank2 <- DropletUtils::emptyDrops(scmat, BPPARAM = BPPARAM, retain=emptydrops.retain)
      scmat_filtered <- scmat[, which(bc_rank2$FDR < emptydrops.fdr)]
      if(is.null(dim(scmat_filtered))){
        message("emptyDrops find all droplets as empty! emptyDrops don't used!")
        emptydrops.fdr <- NULL #for graphs
        umi.kept.nb <- umi.total.nb
        plot_emptydrops <- FALSE
      }else{
        scmat <- scmat_filtered
        message('Droplets matrix dimensions (filtered) :')
        print(dim(scmat))
        message('Total UMIs (filtered) :')
        umi.kept.nb <- sum(scmat)
        print(umi.kept.nb)
        message('Fraction of UMIs in cells :')
        print(umi.kept.nb / umi.total.nb)
        plot_emptydrops <- TRUE
      }
      ## Cleaning
      rm(bc_rank2)
      rm(scmat_filtered)
    } else {
      umi.kept.nb <- umi.total.nb
      plot_emptydrops <- FALSE
    }

    ## plot saturation and Kneeplot
    if(draw_plots){
      nb_umi_genes_by_barcode$droplets_state = "Empty Droplets"
      nb_umi_genes_by_barcode[nb_umi_genes_by_barcode$barcodes %in% colnames(scmat), "droplets_state"] ="Full Droplets"
  
      kneeplot <- ggplot2::ggplot(nb_umi_genes_by_barcode, ggplot2::aes(y=nb_umi, x=num_barcode, color=droplets_state)) +
        ggplot2::geom_point() + ggplot2::ggtitle(sample.name) + ggplot2::theme(legend.title = ggplot2::element_blank()) +
        ggplot2::scale_y_log10(name = "Number of umi by droplet (log scale)") + ggplot2::scale_x_log10(name = "Droplet rank (log scale)") +
        ggplot2::expand_limits(x = 0, y = 0)
      if(plot_emptydrops){
        kneeplot <- kneeplot +
          ggplot2::geom_hline(ggplot2::aes(yintercept = bc_rank@metadata$knee, colour = "Knee Line"), linetype = "solid", show.legend = FALSE) +
          ggplot2::geom_hline(ggplot2::aes(yintercept = bc_rank@metadata$inflection, colour = "Inflection Line"), linetype = "solid", show.legend = FALSE) +
          ggplot2::scale_colour_manual(values = c("cyan3","royalblue4","slateblue3","deeppink3"), guide='legend') +
          ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(linetype = c(0, 0, 1, 1), shape = c(16, 16, 16, 16))))
      }else{
        kneeplot <- kneeplot +
          ggplot2::scale_colour_manual(values = c("cyan3","royalblue4"))
      }
      ggplot2::ggsave(paste0(out.dir, sample.name, "_kneeplot.png"), plot = kneeplot, width = 7, height = 5)
      ## Cleaning
      rm(bc_rank)
  
      saturation_plot <- ggplot2::ggplot(nb_umi_genes_by_barcode, ggplot2::aes(y = nb_genes ,x = nb_umi, color = droplets_state)) +
        ggplot2::geom_point() + ggplot2::ggtitle(sample.name) + ggplot2::theme(legend.title = ggplot2::element_blank()) +
        ggplot2::geom_smooth(colour = "red") +
        ggplot2::scale_y_log10(name = "Number of genes by droplet (log scale)") + ggplot2::scale_x_log10(name = "Number of umi by droplet (log scale)") +
        ggplot2::expand_limits(x = 0, y = 0) +
        ggplot2::scale_colour_manual(values = c("cyan3","royalblue3"))
      ggplot2::ggsave(paste0(out.dir, sample.name, "_saturation_plot.png"), plot = saturation_plot, width = 7, height = 5)
      rm(nb_umi_genes_by_barcode)
    }
    
    #return matrix
    if (return.matrix) return(scmat)

    ## Creation of the Seurat object and save parameters
    sobj <- Seurat::CreateSeuratObject(counts = scmat, project = sample.name, assay = assay)
    rm(scmat)
    sobj[[paste0('log_nCount_', assay)]] <- log(sobj[[paste0('nCount_', assay)]])
    
    ## Read run_info.json from alignement by Kellisto BUStools
    if(file.exists(paste0(data.path, "/run_info.json"))) {
      json_data <- rjson::fromJSON(file=paste0(data.path, "/run_info.json"))
 
    ## Save excel measures
    sobj@misc$excel$Kallisto_Bustools_alignment$Total_reads <- json_data$n_processed
    sobj@misc$excel$Kallisto_Bustools_alignment$Pseudo_aligned_reads <- json_data$n_pseudoaligned
    sobj@misc$excel$Kallisto_Bustools_alignment$Pseudo_aligned_reads_percent <- json_data$p_pseudoaligned
    sobj@misc$excel$Kallisto_Bustools_alignment$Pseudo_aligned_reads_to_unique_target_sequence <- json_data$n_unique
    sobj@misc$excel$Kallisto_Bustools_alignment$Pseudo_aligned_reads_to_unique_target_sequence_percent <- json_data$p_unique
    }
    sobj@misc$excel$Droplet_Quality$captured_droplet <- droplets.nb
    sobj@misc$excel$Droplet_Quality$total_number_UMI <- umi.total.nb
    sobj@misc$excel$Droplet_Quality$estimated_cells <- ncol(sobj)
    sobj@misc$excel$Droplet_Quality$estimated_UMI <- umi.kept.nb
    sobj@misc$excel$Droplet_Quality$fraction_read_in_cells <- umi.kept.nb / umi.total.nb

    ## Save parameters
    sobj@misc$params$sobj_creation$emptydrops.fdr <- emptydrops.fdr
    sobj@misc$params$sobj_creation$droplets.limit <- droplets.limit
    sobj@misc$params$sobj_creation$emptydrops.fdr <- emptydrops.fdr
    sobj@misc$params$sobj_creation$emptydrops.retain <- emptydrops.retain
    sobj@misc$params$sobj_creation$translation <- translation
    sobj@misc$params$sobj_creation$translation.file <- translation.file
    sobj@misc$params$sobj_creation$Rsession <- utils::capture.output(devtools::session_info())
    sobj@misc$params$seed <- my.seed
    
    ## Save command
    sobj@misc$pipeline_commands <- paste0("load.sc.data(data.path = ", data.path, ", sample.name = ", sample.name, ", assay = ", assay, ", droplets.limit = ", droplets.limit, ", emptydrops.fdr = ", emptydrops.fdr, ", emptydrops.retain = ", emptydrops.retain, ", return.matrix = ", return.matrix, ", translation = ", translation, ",  translation.file = ", translation.file, ", BPPARAM = BiocParallel::SerialParam(), my.seed = 1337, out.dir = ", out.dir, ")")
    
    ## Save packages versions
    if(file.exists(paste0(data.path, "/run_info.json"))) sobj@misc$technical_info$kallisto <- json_data$kallisto_version
    sobj@misc$technical_info$BUSpaRse <- utils::packageVersion('BUSpaRse')
    sobj@misc$technical_info$DropletUtils <- utils::packageVersion('DropletUtils')
    sobj@misc$technical_info$Seurat <- utils::packageVersion('Seurat')
    
    ## Save Materials&Methods
    if(file.exists(paste0(data.path, "/Materials_and_Methods.txt"))){
      tmp <- readr::read_tsv(paste0(data.path, "/Materials_and_Methods.txt"), col_names = FALSE)$X1
      tmp2 <- ""
      for (i in 1:length(tmp)) tmp2=paste(tmp2,tmp[i], sep="")
      sobj@misc$parameters$Materials_and_Methods$part0_Alignment <- tmp2
    }
    
    return(sobj)
  } else stop('Data source does not exist, or no sample name specified !')
}

## Basic QC metrics
QC.metrics <- function(sobj = NULL, assay ='RNA', mt.genes.file = NULL, crb.genes.file = NULL, str.genes.file = NULL, pcmito.range = c(0, .1), pcribo.range = c(0, 1), min.features = 200, min.counts = 1000, nbin = 10, BPPARAM = BiocParallel::SerialParam()) {
  if(!is.null(sobj)) {
    if(!is.null(mt.genes.file)) if (!file.exists(mt.genes.file)) stop('mt.genes.file not found !')
    if(!is.null(crb.genes.file)) if (!file.exists(crb.genes.file)) stop('crb.genes.file not found !')
    if(!is.null(str.genes.file)) if (!file.exists(str.genes.file)) stop('str.genes.file not found !')

    ## Save command
    sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("QC.metrics(sobj = sobj, mt.genes.file = ", mt.genes.file, ", crb.genes.file = ", crb.genes.file, ", str.genes.file = ", str.genes.file, ", pcmito.range = c(", pcmito.range[1], ",", pcmito.range[2], ") , pcribo.range = c(", pcribo.range[1], ",", pcribo.range[2], ") , min.features = ", min.features, ", min.counts = ", min.counts, ", nbin = ", nbin, ", BPPARAM = BiocParallel::SerialParam())"))

    ## Restoring seed
    my.seed <- sobj@misc$params$seed

    ### percentage of counts in the top features
    pcQC <- scater::perCellQCMetrics(Seurat::as.SingleCellExperiment(sobj), BPPARAM = BPPARAM)
    sobj@meta.data <- cbind(sobj@meta.data, as.data.frame(pcQC[,grep("percent", colnames(pcQC))]))

    ### MITO
    mito.symbols <- if (!is.null(mt.genes.file)) readRDS(mt.genes.file) else if (!is.null(sobj@misc$params$QC$mito.symbols)) sobj@misc$params$QC$mito.symbols else NULL
    if (!is.null(mito.symbols)) {
      ## Manual percent
      sobj@misc$params$QC$mito.symbols = mito.symbols
      sobj@misc$params$QC$pcmito.range = pcmito.range
      inmito <- rownames(sobj@assays$RNA@counts) %in% mito.symbols
      sobj$percent_mt <- as.vector(Matrix::colSums(sobj@assays$RNA@counts[inmito,]) / sobj$nCount_RNA)
      message('% Mitochondrial expression :')
      print(summary(sobj$percent_mt))
      pcmito_leftbound <- sobj$percent_mt >= pcmito.range[1]
      pcmito_rightbound <- sobj$percent_mt <= pcmito.range[2]
      sobj$pcmito_inrange <- pcmito_leftbound & pcmito_rightbound
      message(paste0(pcmito.range[1], ' <= % mito <= ', pcmito.range[2], ' :'))
      pcmito_factor <- factor(as.numeric(pcmito_leftbound) + abs(as.numeric(pcmito_rightbound)-1))
      levels(pcmito_factor)[levels(pcmito_factor) == 0] <- "out.left"
      levels(pcmito_factor)[levels(pcmito_factor) == 1] <- "in"
      levels(pcmito_factor)[levels(pcmito_factor) == 2] <- "out.right"
      print(table(pcmito_factor))
      ## Seurat AddModuleScore
      # sobj@meta.data$MTscore <- Seurat::AddModuleScore(object = sobj, features = list(mito.symbols), nbin = nbin, seed = my.seed, assay = 'RNA', name = 'SCORE')@meta.data$SCORE1
      if(is.null(sobj@misc$excel$Cells_Quality$mito_summary)) sobj@misc$excel$Cells_Quality$mito_summary <- summary(sobj$percent_mt) else sobj@misc$excel$Final_Cells_Quality$mito_summary <- summary(sobj$percent_mt)
    }
    ### RIBO
    ribo.symbols <- if (!is.null(crb.genes.file)) readRDS(crb.genes.file) else if (!is.null(sobj@misc$params$QC$ribo.symbols)) sobj@misc$params$QC$ribo.symbols else NULL
    if (!is.null(ribo.symbols)) {
      ## Manual percent
      sobj@misc$params$QC$ribo.symbols = ribo.symbols
      sobj@misc$params$QC$pcribo.range = pcribo.range
      inribo <- rownames(sobj@assays$RNA@counts) %in% ribo.symbols
      sobj$percent_rb <- as.vector(Matrix::colSums(sobj@assays$RNA@counts[inribo,]) / sobj$nCount_RNA)
      message('% Ribosomal expression :')
      print(summary(sobj$percent_rb))
      pcribo_leftbound <- sobj$percent_rb >= pcribo.range[1]
      pcribo_rightbound <- sobj$percent_rb <= pcribo.range[2]
      sobj$pcribo_inrange <- pcribo_leftbound & pcribo_rightbound
      message(paste0(pcribo.range[1], ' <= % ribo <= ', pcribo.range[2], ' :'))
      pcribo_factor <- factor(as.numeric(pcribo_leftbound) + abs(as.numeric(pcribo_rightbound)-1))
      levels(pcribo_factor)[levels(pcribo_factor) == 0] <- "out.left"
      levels(pcribo_factor)[levels(pcribo_factor) == 1] <- "in"
      levels(pcribo_factor)[levels(pcribo_factor) == 2] <- "out.right"
      print(table(pcribo_factor))
      ## Seurat AddModuleScore
      # sobj@meta.data$RBscore <- Seurat::AddModuleScore(object = sobj, features = list(ribo.symbols), nbin = nbin, seed = my.seed, assay = 'RNA', name = 'SCORE')@meta.data$SCORE1
      if(is.null(sobj@misc$excel$Cells_Quality$ribo_summary)) sobj@misc$excel$Cells_Quality$ribo_summary <- summary(sobj$percent_rb) else sobj@misc$excel$Final_Cells_Quality$ribo_summary <- summary(sobj$percent_rb)
    }
    ### STRESS
    stress.symbols <- if (!is.null(str.genes.file)) readRDS(str.genes.file) else if (!is.null(sobj@misc$params$QC$stress.symbols)) sobj@misc$params$QC$stress.symbols else NULL
    if (!is.null(stress.symbols)) {
      ## Manual percent
      sobj@misc$params$QC$stress.symbols = stress.symbols
      instress <- rownames(sobj@assays$RNA@counts) %in% stress.symbols
      sobj$percent_st <- as.vector(Matrix::colSums(sobj@assays$RNA@counts[instress,]) / sobj$nCount_RNA)
      message('% Stress response expression :')
      print(summary(sobj$percent_st))
      ## Seurat AddModuleScore
      # sobj@meta.data$STscore <- Seurat::AddModuleScore(object = sobj, features = list(stress.symbols), nbin = nbin, seed = my.seed, assay = 'RNA', name = 'SCORE')@meta.data$SCORE1
      if(is.null(sobj@misc$excel$Cells_Quality$stress_summary)) sobj@misc$excel$Cells_Quality$stress_summary <- summary(sobj$percent_st) else sobj@misc$excel$Final_Cells_Quality$stress_summary <- summary(sobj$percent_st)
    }

    ### NB FEATURES
    sobj@misc$params$QC$min.features = min.features
    sobj$min_features <- sobj$nFeature_RNA >= min.features
    message(paste0('Number of cells with >= ', min.features, ' features :'))
    print(table(sobj$min_features))
    if(is.null(sobj@misc$excel$Cells_Quality$filter_cells_genes)) sobj@misc$excel$Cells_Quality$filter_cells_genes <- sum(sobj$min_features)
    if(is.null(sobj@misc$excel$Cells_Quality$genes_per_cell_summary)) sobj@misc$excel$Cells_Quality$genes_per_cell_summary <- summary(Matrix::colSums(sobj@assays$RNA@counts != 0)) else sobj@misc$excel$Final_Cells_Quality$genes_per_cell_summary <- summary(Matrix::colSums(sobj@assays$RNA@counts != 0))

    ### NB COUNTS
    sobj@misc$params$QC$min.counts = min.counts
    sobj$min_counts <- sobj$nCount_RNA >= min.counts
    message(paste0('Number of cells with >= ', min.counts, ' of total counts:'))
    print(table(sobj$min_counts))
    if(is.null(sobj@misc$excel$Cells_Quality$filter_cells_counts)) sobj@misc$excel$Cells_Quality$filter_cells_counts <- sum(sobj$min_counts)
    if(is.null(sobj@misc$excel$Cells_Quality$UMI_per_cell_summary)) sobj@misc$excel$Cells_Quality$UMI_per_cell_summary <- summary(Matrix::colSums(sobj@assays$RNA@counts)) else sobj@misc$excel$Final_Cells_Quality$UMI_per_cell_summary <- summary(Matrix::colSums(sobj@assays$RNA@counts))

    ## Save parameters
    sobj@misc$params$QC$Rsession <- utils::capture.output(devtools::session_info())
    ## Save packages versions
    sobj@misc$technical_info$scater <- utils::packageVersion('scater')
    
  }
  return(sobj)
}

## Plot QC histograms (after using QC.metrics()  or cells.QC.filters())
QC.hist <- function(sobj = NULL, assay = 'RNA', out.dir = NULL) {
  if ("Seurat" %in% is(sobj) && !is.null(out.dir) && assay %in% names(sobj@assays)) {
    require(patchwork)
    sample.name <- Seurat::Project(sobj)
    qcplots <-list(
      histFEAT <- ggplot2::qplot(sobj[[paste(c('nFeature', assay), collapse = '_'), drop = TRUE]], geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main=paste0(paste(c('nFeature', assay), collapse = '_'), ' (>= ', sobj@misc$params$QC$min.features, " : ", length(which(sobj$min_features)), " cells)"), xlab = paste(c('nFeature', assay), collapse = '_')) + ggplot2::geom_vline(xintercept = sobj@misc$params$QC$min.features, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme(),
      histNC <- ggplot2::qplot(sobj[[paste(c('nCount', assay), collapse = '_'), drop = TRUE]], geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main=paste0(paste(c('nCount', assay), collapse = '_'), ' (>= ', sobj@misc$params$QC$min.counts, " : ", length(which(sobj$min_counts)), " cells)"), xlab = paste(c('nCount', assay), collapse = '_')) + ggplot2::geom_vline(xintercept = sobj@misc$params$QC$min.counts, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme(),
      histMT <- ggplot2::qplot(sobj$percent_mt, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main=paste0("%mito (in [", sobj@misc$params$QC$pcmito.range[1], ';', sobj@misc$params$QC$pcmito.range[2], "] : ", length(which(sobj$pcmito_inrange)), " cells)"), xlab = "%mito") + ggplot2::geom_vline(xintercept=c(0.05,0.10,0.15,0.20), col = "cyan4", linetype = "dashed", size = 1) + ggplot2::geom_vline(xintercept=sobj@misc$params$QC$pcmito.range, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme(),
      histRB <- ggplot2::qplot(sobj$percent_rb, geom = "histogram", bins = 101, fill = I("white"), col = I("black"), main=paste0("%ribo (in [", sobj@misc$params$QC$pcribo.range[1], ';', sobj@misc$params$QC$pcribo.range[2], "] : ", length(which(sobj$pcribo_inrange)), " cells)"), xlab = "%ribo") + ggplot2::geom_vline(xintercept=sobj@misc$params$QC$pcribo.range, col = "red", linetype = "dashed", size = 2) + Seurat::DarkTheme()
    )
    png(paste0(out.dir, '/', sample.name, '_QChist.png'), width = 1000, height = 1000)
    print(patchwork::wrap_plots(qcplots))
    dev.off()
  }
  ## Save packages versions
  sobj@misc$technical_info$ggplot2 <- utils::packageVersion('ggplot2')
}

## Filter cells based on QC metrics
cells.QC.filter <- function(sobj = NULL, min.features = 200, min.counts = 1000, pcmito.range = c(0, .1), pcribo.range = c(0, 1)) {
  if (!is.null(sobj) & all(c('nFeature_RNA', 'nCount_RNA', 'percent_mt', 'percent_rb') %in% colnames(sobj@meta.data))) {
    ## Save command
    sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("cells.QC.filter(sobj = sobj, min.features = ", min.features, ", min.counts = ", min.counts, ", pcmito.range = c(", pcmito.range[1], ",", pcmito.range[2], ") , pcribo.range = c(", pcribo.range[1], ",", pcribo.range[2], "))"))
    ## Fitrers processing
    message('Cell expression matrix dimensions (unfiltered) :')
    print(dim(sobj))
    cellskeep <- sobj$nFeature_RNA >= min.features & sobj$nCount_RNA >= min.counts & sobj$percent_mt >= pcmito.range[1] & sobj$percent_mt <= pcmito.range[2] & sobj$percent_rb >= pcribo.range[1] & sobj$percent_rb <= pcribo.range[2]
    sobj <- sobj[, cellskeep]
    message('Cell expression matrix dimensions (filtered) :')
    print(dim(sobj))
    ## Save parameters
    sobj@misc$params$QC$pcmito.range = pcmito.range
    sobj@misc$params$QC$pcribo.range = pcribo.range
    sobj@misc$params$QC$min.features = min.features
    sobj@misc$params$QC$min.counts = min.counts
    sobj@misc$params$QC$Rsession <- utils::capture.output(devtools::session_info())
  }
  return(sobj)
}

## Predict cell cycle phase
cell.cycle.predict <- function(sobj = NULL, assay = 'RNA', cc.cyclone.file = NULL, cc.seurat.file = NULL, BPPARAM = BiocParallel::SerialParam(), nbin = 24) {
  if(!is.null(sobj)) {
    if (!file.exists(cc.cyclone.file)) stop('cc.cyclone.file not found !')
    if (!file.exists(cc.seurat.file)) stop('cc.seurat.file not found !')
    if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))

    ## Save command
    sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("cell.cycle.predict(sobj = sobj, assay = ", assay, ", cc.cyclone.file = ", cc.cyclone.file, ", cc.seurat.file = ", cc.seurat.file, ", BPPARAM = BiocParallel::SerialParam(), nbin = ", nbin, ")"))

    ## Restoring seed
    my.seed <- sobj@misc$params$seed

    ## Load files
    cc_cyclone <- readRDS(cc.cyclone.file)
    cc_seurat <- readRDS(cc.seurat.file)

    ### Seurat
    sobj <- Seurat::CellCycleScoring(object = sobj, s.features = cc_seurat$s.genes, g2m.features = cc_seurat$g2m.genes, assay = assay, nbin = nbin, seed = my.seed)
    sobj$Seurat.S.Score <- sobj$S.Score
    sobj$Seurat.G2M.Score <- sobj$G2M.Score
    sobj$Seurat.SmG2M.Score <- sobj$S.Score - sobj$G2M.Score
    sobj$Seurat.Phase <- sobj$Phase
    sobj$S.Score <- sobj$G2M.Score <- sobj$Phase <- NULL
    message("Cell cycle phases according to Seurat: ")
    print(table(sobj$Seurat.Phase))

    ### Scran: Cyclone
    set.seed(my.seed)
    cycres <- scran::cyclone(Seurat::as.SingleCellExperiment(sobj, assay = assay), pairs=cc_cyclone, BPPARAM = BPPARAM, verbose = TRUE)
    sobj$Cyclone.Phase <- as.factor(cycres$phases)
    sobj$Cyclone.nG1.Score <- cycres$normalized.scores$G1
    sobj$Cyclone.nS.Score <- cycres$normalized.scores$S
    sobj$Cyclone.nG2M.Score <- cycres$normalized.scores$G2M
    sobj$Cyclone.nSmG2M.Score <- cycres$normalized.scores$S - cycres$normalized.scores$G2M
    sobj$Cyclone.G1.Score <- cycres$scores$G1
    sobj$Cyclone.S.Score <- cycres$scores$S
    sobj$Cyclone.G2M.Score <- cycres$scores$G2M
    sobj$Cyclone.SmG2M.Score <- cycres$scores$S - cycres$scores$G2M
    message("Cell cycle phases according to Cyclone: ")
    print(table(sobj$Cyclone.Phase))
    sobj@misc$excel$After_QC_cells_filtering$estim_cells_G1 <- sum(sobj$Cyclone.Phase == "G1")
    sobj@misc$excel$After_QC_cells_filtering$estim_cells_G2M <- sum(sobj$Cyclone.Phase == "G2M")
    sobj@misc$excel$After_QC_cells_filtering$estim_cells_S <- sum(sobj$Cyclone.Phase == "S")

    ## Save parameters
    sobj@misc$params$QC$cell.cycle$cyclone.cell.cycle.genes = cc_cyclone;
    sobj@misc$params$QC$cell.cycle$seurat.cell.cycle.genes = cc_seurat;
    sobj@misc$params$QC$cell.cycle$Rsession <- utils::capture.output(devtools::session_info())
    ## Save packages versions
    sobj@misc$technical_info$scran <- utils::packageVersion('scran')
  }
  return(sobj)
}

## Features filtering on fixed criteria
features.filter <- function(sobj = NULL, min.cells = 5)  {
  if(!is.null(sobj)) {
    ## Save command
    sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("features.filter(sobj = sobj, min.cells = ", min.cells, ")"))

    if(all(Seurat::Assays(sobj) == 'RNA'))  {
      ## Filer processing
      message('Cell expression matrix dimensions (unfiltered) :')
      print(dim(sobj))
      featureskeep <- which(apply(sobj@assays$RNA@counts,1,function(x){length(which(x>0))}) >= min.cells)
      sobj <- sobj[featureskeep,]
      message('Cell expression matrix dimensions (filtered) :')
      print(dim(sobj))
      sobj@misc$excel$After_QC_cells_filtering$Genes_covered <- dim(sobj)[1]
      ## Save parameters
      sobj@misc$params$QC$min.cells <- min.cells
      sobj@misc$Rsession <- sessionInfo()
    } else {
      message("WARNING : Seurat object contained other assays than 'RNA', no filtering performed.")
    }
  }
  return(sobj)
}

## Finding cell doublets
find.doublets <- function(sobj = NULL, assay = 'RNA') {
  if (!is.null(sobj)) {
    if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
    ## Restoring seed
    my.seed <- sobj@misc$params$seed

    ##Save command
    sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("find.doublets(sobj = sobj, assay = ", assay, ")"))

    ## scDblFinder
    set.seed(my.seed)
    sobj$scDblFinder.class <- Seurat::as.Seurat(scDblFinder::scDblFinder(Seurat::as.SingleCellExperiment(sobj, assay = assay)))$scDblFinder.class ## Do not try using BPPARAM, it only works when computing on multiple samples
    sobj$scDblFinder.class <- unname(sobj$scDblFinder.class == "doublet")
    message('scDblFinder doublets :')
    print(table(sobj$scDblFinder.class))

    # scds
    set.seed(my.seed)
    sobj$hybrid_score <- scds::cxds_bcds_hybrid(Seurat::as.SingleCellExperiment(sobj, assay = assay))$hybrid_score
    sobj$hybrid_score.class <- unname(sobj$hybrid_score > 1)
    message('scds-hybrid doublets :')
    print(table(sobj$hybrid_score.class))

    ## union scDblFinder & scds
    sobj$doublets_consensus.class <- sobj$scDblFinder.class | sobj$hybrid_score.class
    message('Consensus doublets :')
    print(table(sobj$doublets_consensus.class))
    sobj@misc$excel$After_QC_cells_feature_filtering$estim_doublets <- sum(sobj$doublets_consensus.class)
    sobj@misc$doublets <- list(scDblFinder = length(which(sobj$scDblFinder.class)), scds_hybrid = length(which(sobj$hybrid_score.class)), union = length(which(sobj$doublets_consensus.class)))

    ## Save parameters
    sobj@misc$params$doublets$Rsession <- utils::capture.output(devtools::session_info())
    ## Save packages versions
    sobj@misc$technical_info$scDblFinder <- utils::packageVersion('scDblFinder')
    sobj@misc$technical_info$scds <- utils::packageVersion('scds')
  }
  return(sobj)
}

## Filtering doublets
filter.doublets <- function(sobj = NULL, method = "all") { ## Method can be 'both', 'scDblFinder', 'scds'
  if (!is.null(sobj) && method %in% c('all', 'scDblFinder', 'scds')) {
    ## Save command
    sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("filter.doublets(sobj = sobj, method = ", method, ")"))
    ## Filter processing
    message('Cell expression matrix dimensions (unfiltered) :')
    print(dim(sobj))
    if (method == 'all') sobj <- sobj[, !sobj$doublets_consensus.class]
    if (method == 'scDblFinder') sobj <- sobj[, !sobj$scDblFinder.class]
    if (method == 'scds') sobj <- sobj[, !sobj$hybrid_score.class]
    message('Cell expression matrix dimensions (filtered) :')
    print(dim(sobj))
    ## Save parameters
    sobj@misc$params$doublets$method_filtering <- method

  }
  return(sobj)
}

## Control genes (if any)
tag.ctrl.genes <- function(sobj = NULL, ctrl.genes = c("GAPDH"), ctrl.min.counts = 3, assay = "RNA") {
  if (!is.null(ctrl.genes) && !is.null(sobj)) {
    if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
    ##Save command
    sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("tag.ctrl.genes(sobj = sobj, ctrl.genes = c(", paste(ctrl.genes, collapse = ","), ") , ctrl.min.counts = ", ctrl.min.counts, ", assay = ", assay, ")"))

    if (length(ctrl.genes) > 0) {
      ok.ctrl <- ctrl.genes %in% rownames(sobj@assays[[assay]]@counts)
      if (any(ok.ctrl)) for (ctrlg in ctrl.genes[ctrl.genes %in% rownames(sobj@assays[[assay]]@counts)]) suppressMessages(sobj[[paste0('ctrl_', assay, '_', ctrlg)]] <- sobj@assays[[assay]]@data[rownames(sobj@assays[[assay]]@data) == ctrlg] >= ctrl.min.counts)

      ## Save parameters
      sobj@misc$params$ctrl.genes <- ctrl.genes
      sobj@misc$params$ctrl.min.counts <- ctrl.min.counts
    }
  }
  return(sobj)
}

## Normalization (Seurat: SCTransform ou LogNormalize ou CLR)
sc.normalization <- function(sobj = NULL, assay = 'RNA', normalization.method = "SCTransform", features.n = 3000, vtr = NULL) {
  if (!is.null(sobj)) {
    if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
    if(!is.null(vtr)) vtr <- sort(vtr)

    ## Save command
    sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("sc.normalization(sobj = sobj, assay = ", assay, ", normalization.method = ", normalization.method, ", features.n = ", features.n, ", vtr = ", if(is.null(vtr)) "NULL" else paste0("c(", paste(vtr, collapse = ","),")"), ")"))

    ## Restoring seed and set assay.ori
    my.seed <- sobj@misc$params$seed
    assay.ori <- assay

    if (toupper(normalization.method) == toupper("SCTransform")) {
      sobj <- Seurat::SCTransform(object = sobj, assay = assay, seed.use = my.seed, variable.features.n = features.n, vars.to.regress = vtr, return.only.var.genes = TRUE)
      assay <- 'SCT'
    } else if (toupper(normalization.method) == toupper("LogNormalize")) {
      sobj <- Seurat::NormalizeData(sobj, normalization.method = 'LogNormalize', assay = assay)
      sobj <- Seurat::FindVariableFeatures(sobj, assay = assay, nfeatures = features.n)
    } else if (toupper(normalization.method) == toupper("CLR")) {
      sobj <- Seurat::NormalizeData(sobj, normalization.method = 'CLR', assay = assay)
      sobj <- Seurat::FindVariableFeatures(sobj, assay = assay, nfeatures = features.n)
    } else stop('Unknown or unsupported normalization method !')

    ## Save parameters
    sobj@assays[[assay]]@misc$params$normalization <- list(normalization.method = normalization.method, assay.ori = assay.ori, assay.out = assay, features.used = features.n)
    sobj@assays[[assay]]@misc$scaling = list(vtr = if(is.null(vtr)) NA else vtr)
    sobj@misc$params$normalization <- c(sobj@assays[[assay]]@misc$params$normalization, sobj@assays[[assay]]@misc$scaling)
    sobj@misc$params$normalization$Rsession <- utils::capture.output(devtools::session_info())
  }
  return(sobj)
}

## Dimensions reduction (now, character or factors are converted to a model.matrix)
dimensions.reduction <- function(sobj = NULL, reduction.method = 'pca', assay = 'RNA', max.dims = 100L, vtr = NULL, vtr.scale = TRUE, red.name = NULL) {
  raw.methods <- c('scbfa', 'bpca')
  all.methods <- c(raw.methods, 'pca', 'mds', 'ica')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (!reduction.method %in% all.methods) stop(paste0('Unknown reduction method ! Expecting any of : ', paste(all.methods, collapse = ', ')))
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if(is.null(red.name)) red.name <- paste0(assay, '_', reduction.method)
  if(!is.null(vtr)) vtr <- sort(vtr)

  ## Save command
  sobj@misc$pipeline_commands = c(sobj@misc$pipeline_commands, paste0("dimensions.reduction(sobj = sobj, reduction.method = ", reduction.method, ", assay = ", assay, ", max.dims = ", max.dims, ", vtr = ", if(is.null(vtr)) "NULL" else paste0("c(", paste(vtr, collapse = ","),")"), ", vtr.scale =", vtr.scale, ")"))

  ## Restoring seed
  my.seed <- sobj@misc$params$seed

  if(reduction.method %in% raw.methods) {
    ## Formatting covariates for scbfa
    if (!is.null(vtr)) {
      if (!all(vtr %in% colnames(sobj@meta.data))) stop('Not all vtr names found in Seurat object meta data !')
      minimeta <- sobj@meta.data[,colnames(sobj@meta.data) %in% vtr, drop = FALSE]
      X <- matrix(ncol = 0, nrow = nrow(sobj@meta.data))
      for(v in vtr) {
        if(is.character(minimeta[[v]])) minimeta[[v]] <- as.factor(minimeta[[v]])
        if(any(is.na(minimeta[[v]]))) stop(paste0("Covariate '", v, "' contains NA value(s) !"))
        if(is.factor(minimeta[[v]])) {
          message(paste0("Converting '", v, "' factor into model matrix and adding to the regression..."))
          mm <- model.matrix(~minimeta[[v]])[,-1, drop = FALSE]
          # mm <- model.matrix(~minimeta[[v]])
          X <- cbind(X, mm)
        } else {
          message(paste0("Adding '", v, "' covariate", if(vtr.scale) ' (scaled)' else NULL,  ' to the regression ...'))
          X <- cbind(X, if(vtr.scale) scale(minimeta[[v]]) else minimeta[[v]])
        }
      }
    } else {
      X <- NULL
    }
    # ## Formatting covariates if needed
    # if (!is.null(vtr)) {
    #   if(is.matrix(vtr)) {
    #     message(paste0("Regressing matrix named '", vtr.matrix.name, "' ..."))
    #     X <- vtr
    #   } else {
    #     if (!all(vtr %in% colnames(sobj@meta.data))) stop('Not all vtr names found in Seurat object meta data !')
    #     X <- sobj@meta.data[, vtr, drop = FALSE]
    #     for(x in seq_len(ncol(X))) {
    #       if (!is.numeric(X[,x])) X[,x] <- as.numeric(as.factor(X[,x]))
    #     }
    #     X <- as.matrix(X)
    #   }
    #   if(vtr.scale) X <- scale(X)
    # } else X <- vtr
  } else {
    ## If vtr is called, stop !
    if(!is.null(vtr)) stop("vtr can only be called when using any reduction method out of : ", paste(raw.methods, collapse = ', '))
    ## Scaling if necessary
    if (sum(dim(sobj@assays[[assay]]@scale.data)) < 3) {
      scale.vtr.all <- NULL
      if(!any(is.na(sobj@assays[[assay]]@misc$scaling$vtr))) {
        message(paste0("Found scaling coveriate(s) '", paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = "', '"), "' to regress from normalization ..."))
        scale.vtr.all <- c(scale.vtr.all, sobj@assays[[assay]]@misc$scaling$vtr)
      }
      assay.ori <- Seurat::DefaultAssay(sobj)
      Seurat::DefaultAssay(sobj) <- assay
      if(assay == 'SCT') {
        sobj <- Seurat::ScaleData(object = sobj,
                                  vars.to.regress = scale.vtr.all,
                                  do.scale = FALSE, scale.max = Inf, block.size = 750)
      } else {
        sobj <- Seurat::ScaleData(object = sobj,
                                  vars.to.regress = scale.vtr.all,
                                  do.scale = TRUE, scale.max = 10, block.size = 1000)
      }
      Seurat::DefaultAssay(sobj) <- assay.ori
    }
  }

  if (reduction.method == 'pca') {
    sobj <- Seurat::RunPCA(object = sobj, assay = assay, verbose = FALSE, npcs = max.dims, reduction.name = red.name, reduction.key = paste0(red.name, '_'), seed.use = my.seed)
    sobj@assays[[assay]]@misc$params$reductions$vtr <- sobj@reductions[[red.name]]@misc$vtr <- NA
  } else if (reduction.method == 'ica') {
    sobj <- Seurat::RunICA(object = sobj, assay = assay, verbose = FALSE, nics = max.dims, reduction.name = red.name, reduction.key = paste0(red.name, '_'), seed.use = my.seed)
    sobj@assays[[assay]]@misc$params$reductions$vtr <- sobj@reductions[[red.name]]@misc$vtr <- NA
  } else if (reduction.method == 'mds') {
    set.seed(my.seed)
    mds.res <- scater::calculateMDS(x = sobj@assays[[assay]]@scale.data, ncomponents = max.dims)
    sobj@reductions[[red.name]] <- Seurat::CreateDimReducObject(embeddings = mds.res, loadings = matrix(nrow = 0, ncol = 0), assay = assay, stdev = matrixStats::colSds(mds.res), key = paste0(red.name, '_'), misc = list())
    sobj@assays[[assay]]@misc$params$reductions$vtr <- sobj@reductions[[red.name]]@misc$vtr <- NA
    ## Save packages versions
    sobj@misc$technical_info$scater <- utils::packageVersion('scater')
  } else if (reduction.method == 'bpca') {
    set.seed(my.seed)
    bpca.res <- scBFA::BinaryPCA(scData = as.matrix(sobj@assays[[assay]]@counts[sobj@assays[[assay]]@var.features,]), X = X)
    colnames(bpca.res$x) <- paste0('BPCA_', 1L:features.n)
    colnames(bpca.res$rotation) <- paste0('BPCA_', 1L:features.n)
    sobj@reductions[[red.name]] <- Seurat::CreateDimReducObject(embeddings = bpca.res$x, loadings = bpca.res$rotation, assay = assay, stdev = bpca.res$sdev, key = paste0(red.name, '_'), misc = list(center = bpca.res$center, scale = bpca.res$scale))
    ## Save parameters
    sobj@assays[[assay]]@misc$params$reductions <- sobj@reductions[[red.name]]@misc <- list(numFactors = bpca.res$numFactors,
                                                                                                       X = if(is.null(X)) NA else X,
                                                                                                       vtr = if(is.null(vtr)) NA else vtr,
                                                                                                       vtr.scale = vtr.scale)
    ## Save packages versions
    sobj@misc$technical_info$scBFA <- utils::packageVersion('scBFA')
    rm(bpca.res)

  } else if (reduction.method == 'scbfa') {
    set.seed(my.seed)
    bfa.res <- scBFA::scBFA(scData = as.matrix(sobj@assays[[assay]]@counts[sobj@assays[[assay]]@var.features,]), numFactors = max.dims, X = X)
    dimnames(bfa.res$ZZ) <- list(colnames(sobj@assays[[assay]]@counts), paste0('SCBFA_', 1L:max.dims))
    dimnames(bfa.res$AA) <- list(sobj@assays[[assay]]@var.features, paste0('SCBFA_', 1L:max.dims))
    sobj@reductions[[red.name]] <- Seurat::CreateDimReducObject(embeddings = bfa.res$ZZ, loadings = bfa.res$AA, assay = assay, stdev = matrixStats::colSds(bfa.res$ZZ), key = paste0(red.name, '_'), misc = list())
    ## Save parameters
    sobj@assays[[assay]]@misc$params$reductions <- sobj@reductions[[red.name]]@misc <- list(binary.matrix = bfa.res$BB,
                                                                                                       numFactors = bfa.res$numFactors,
                                                                                                       X = if(is.null(X)) NA else X,
                                                                                                       vtr = if(is.null(vtr)) NA else vtr,
                                                                                                       vtr.scale = vtr.scale)
    ## Save packages versions
    sobj@misc$technical_info$scBFA <- utils::packageVersion('scBFA')
    rm(bfa.res)

  } else stop("Unknown reduction method !")

  ## Deleting @scale.data
  sobj@assays[[assay]]@scale.data <- matrix(nrow = 0, ncol = 0)

  ## Filling @misc
  # sobj@reductions[[red.name]]@misc$vtr <- if (!is.null(vtr)) vtr else NA
  sobj@reductions[[red.name]]@misc$from.assay <- assay
  ## Save parameters
  sobj@assays[[assay]]@misc$params$reductions <- c(sobj@assays[[assay]]@misc$params$reductions, list(method = tolower(reduction.method), assay = assay, max.dims = max.dims))
  sobj@misc$params$reductions <- sobj@assays[[assay]]@misc$params$reductions
  sobj@misc$params$reductions$Rsession <- utils::capture.output(devtools::session_info())
  return(sobj)
}

## Harmonize dataset (remove sample effect) for integration
harmonize <- function(sobj = NULL, reduction = NULL, vtr = NULL, harmony.max.iter = 100, ncomp = 20, out.dir = NULL) {
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(reduction)) stop('No reduction name provided !')
  if (!reduction %in% names(sobj@reductions)) stop(paste0("Reduction '", reduction, "' not found !"))
  if (is.null(vtr)) stop("No variable name to regress provided ! At least one (ex : 'orig.ident') is required.")
  if (!all(vtr %in% colnames(sobj@meta.data))) stop('At least one of the variable names to regress was not found Seurat object metadata !')

  ## Save command
  sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("harmonize(sobj = sobj, reduction = ", reduction, " vtr = ", if(is.null(vtr)) "NULL" else paste0("c(", paste(vtr, collapse = ","),")"), ", harmony.max.iter =", harmony.max.iter, ", ncomp =", ncomp, ", out.dir =", out.dir, ")"))

  message(paste0('Harmonizing for : ', paste(vtr, collapse = ', ')))
  set.seed(sobj@misc$params$seed)
  png(paste0(out.dir, '/harmony_convergence_plot.png'), width = 1000, height = 1000)
  harmonized.red <- harmony::HarmonyMatrix(data_mat = sobj@reductions[[reduction]]@cell.embeddings[,1:ncomp], meta_data = sobj@meta.data[, vtr, drop = FALSE], vars_use = vtr, do_pca = FALSE, max.iter.harmony = harmony.max.iter, max.iter.cluster = 1000, plot_convergence = TRUE)
  dev.off()
  harm.red.name <- paste(c(reduction, paste(c(vtr, 'harmonized'), collapse = '.')), collapse = '_')
  sobj@reductions[[harm.red.name]] <- sobj@reductions[[reduction]]
  sobj@reductions[[harm.red.name]]@cell.embeddings <- harmonized.red
  sobj@reductions[[harm.red.name]]@feature.loadings <- sobj@reductions[[harm.red.name]]@feature.loadings[,1:ncomp]
  sobj@reductions[[harm.red.name]]@stdev <- sobj@reductions[[harm.red.name]]@stdev[1:ncomp]
  sobj@reductions[[harm.red.name]]@misc$harmony <- list(vtr = vtr, max.iter = harmony.max.iter)
  sobj@reductions[[harm.red.name]]@misc$harmony$Rsession <- utils::capture.output(devtools::session_info())
  ## Save packages versions
  sobj@misc$technical_info$harmony <- utils::packageVersion('harmony')
  return(sobj)
}

## DecontX processing
### This requires a seurat object with clusters (see celda::decontX)
### This function returns a new seurat object with corrected raw count matrix
decontx.process <- function(sobj, assay = 'RNA', idents = NULL, ...) {
  if (is.null(sobj)) stop("A Seurat object is required !")
  if (is.null(idents)) stop("An ident name is required ! Use get.idents() to get available idents.")

  ## Save command
  sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("decontx.process(sobj = sobj, assay = ", assay, ", idents = ", idents, ", ...)"))

  ## Processing
  my.counts = as.matrix(sobj@assays[[assay]]@counts)
  mode(my.counts) = "integer"
  my.z = as.factor(as.character(as.numeric(sobj@meta.data[[idents]])))
  dx.res <- celda::decontX(counts = my.counts, z = my.z, seed = sobj@misc$params$seed, ...)
  sobj@assays[[assay]]@counts <- as(dx.res$resList$estNativeCounts, "dgCMatrix")

  ## Cleaning
  rm(my.counts, dx.res, my.z)
  ## Save packages versions
  sobj@misc$technical_info$celda <- utils::packageVersion('celda')
  return(sobj)
}

## Reduction dims correlation with bias sources
### If future gives a failure do to RAM object size :
### 1) Kill current future instance :
### future:::ClusterRegistry("stop")
### 2) Increase max transfer size using :
### options(future.globals.maxSize=2*1024^3)
### 3) Re-run dimensions.eval()
dimensions.eval <- function(sobj = NULL, reduction = 'RNA_scbfa', cor.method = 'spearman', meta.names = c('nCount_RNA', 'nFeature_RNA', 'percent_mt', 'MTscore', 'percent_rb', 'RBscore', 'percent_st', 'STscore', "Cyclone.S.Score", "Cyclone.G1.Score", "Cyclone.G2M.Score", "Cyclone.SmG2M.Score"), eval.markers = c('GAPDH'), slot = 'data', max.dims = 100L, out.dir = NULL, nthreads = 1) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (!(reduction %in% names(sobj@reductions))) stop(paste0('Reduction "', reduction, '" not present in the provided Seurat object !'))
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  sample.name <- Seurat::Project(sobj)
  ndims <- min(max.dims, ncol(sobj@reductions[[reduction]]@cell.embeddings))
  if (!dir.exists(out.dir)) dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  if(slot == 'scale.data' && length(eval.markers) > 0 && (assay == "integrated")) stop("scale.data matrix can't be recomputed after integration by seurat")
  if(slot == 'scale.data' && length(eval.markers) > 0 && (!is.null(sobj@misc$params$integration$method)) && sobj@misc$params$integration$method == 'Liger') stop("scale.data matrix can't be recomputed after integration by Liger")
  if(slot == 'scale.data' && length(eval.markers) > 0 && (!is.null(sobj@misc$params$group$Keep.Norm)) && sobj@misc$params$group$Keep.Norm == 'TRUE') stop("scale.data matrix can't be recomputed after grouping data if normalization is kept.")
  
  ## Regenerate full scaled matrix if eval.markers are available
  if (!is.null(eval.markers)) {
    eval.markers <- sort(eval.markers[eval.markers %in% rownames(sobj@assays[[assay]]@data)])
    if ((slot == 'scale.data') && (length(eval.markers) > 0)) {
      ## Scaling if necessary
      if (sum(dim(sobj@assays[[assay]]@scale.data)) < 3) {
        scale.vtr <- NULL
        if(!is.na(sobj@reductions[[reduction]]@misc$vtr)) scale.vtr <- c(scale.vtr, sobj@reductions[[reduction]]@misc$vtr)
        if("harmony" %in% names(sobj@reductions[[reduction]]@misc)) scale.vtr <- c(scale.vtr, sobj@reductions[[reduction]]@misc$harmony$vtr)
  
        if (length(scale.vtr) > 0) future::plan("multiprocess", workers = nthreads, gc = TRUE)
        if (assay == 'SCT') {
          sobj <- Seurat::ScaleData(object = sobj,
                                    vars.to.regress = scale.vtr,
                                    do.scale = FALSE, scale.max = Inf, block.size = 750)
        }
        else {
          sobj <- Seurat::ScaleData(object = sobj,
                                    vars.to.regress = scale.vtr,
                                    do.scale = TRUE, scale.max = 10, block.size = 1000)
        }
        if (length(scale.vtr) > 0) future:::ClusterRegistry("stop")
      }
    }
  }

  ## Empty plot + Parameters for meta.names plotting
  png(paste0(out.dir, '/', sample.name, '_', reduction, '_dims.bias.cor.png'), width = 1400, height = 800)
  suppressWarnings(plot(0, 0, log = 'x', xlim = c(1L, ndims), ylim = c(0,1), type = "n", xlab = paste0(reduction, " dimension"), ylab = paste0(cor.method, " correlation"), xaxs = "i", yaxs = "i"))
  meta.names <- meta.names[meta.names %in% colnames(sobj@meta.data)]
  meta.cols <- if (length(meta.names) > 12) scales::hue_pal()(length(meta.names)) else RColorBrewer::brewer.pal(length(meta.names), name = "Paired")
  ## Calculation and plottinf of meta.names correlations
  for(myf in seq_along(meta.names)) {
    corvec <- vapply(seq_len(ndims), function(k) { (abs(cor(sobj@reductions[[reduction]]@cell.embeddings[,k], sobj[[meta.names[myf]]], method = cor.method))) }, .1)
    lines(corvec, type = "l", ylim = c(0,1), col = meta.cols[myf], lwd = 4)
  }
  ## Calculation and plottinf of eval.markers correlations
  if (!is.null(eval.markers)) {
    mydata <- slot(sobj@assays[[assay]], slot)
    for(mym in seq_along(eval.markers)) {
      corvec <- vapply(seq_len(ndims), function(k) { (abs(cor(sobj@reductions[[reduction]]@cell.embeddings[,k], mydata[rownames(mydata) == eval.markers[mym],], method = cor.method))) }, .1)
      lines(corvec, type = "l", ylim = c(0,1), col = 1, lwd = 5, lty = mym+1L)
    }
  }
  ## Adding legend
  legend(x = ndims, y = 1, legend = c(meta.names, eval.markers), col = c(meta.cols, rep(1, length(eval.markers))), lty = c(rep(1, length(meta.names)), seq_along(eval.markers)), lwd = 5, xjust = 1, yjust = 1)
  dev.off()
}

## Evaluating Louvain clusters (pca dims / resolution) (multithreading version)
clustering.eval.mt <- function(sobj = NULL, reduction = 'RNA_scbfa', dimsvec = seq.int(3,99, 2), resvec = seq(.1,1.2,.1), out.dir = NULL, solo.pt.size = 2, BPPARAM = BiocParallel::SerialParam()) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (!"seed" %in% names(sobj@misc$params)) stop("No seed found in Seurat object !")
  if (!(reduction %in% names(sobj@reductions))) stop(paste0('Reduction "', redution, '" not present in the provided Seurat object !'))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if (max(dimsvec) > ncol(sobj@reductions[[reduction]]@cell.embeddings)) stop(paste0('Max dimsvec requested is ', max(dimsvec), ' whereas max dimension in "', reduction, '" reduction is ', ncol(sobj@reductions[[reduction]]@cell.embeddings)))
  if (1 %in% dimsvec) stop("One should not request a reduction to a single dimension !")

  ## Restoring seed and sample name from within the object
  my.seed <- sobj@misc$params$seed
  sample.name <- Seurat::Project(sobj)

  ## Create output dir
  clustree.dir <- paste0(out.dir, "/clustree_", reduction, '/')
  umaps.clustree.dir <- paste0(clustree.dir, "/uMAPs/")
  pca.clustree.dir <- paste0(clustree.dir, "/dimensions/")
  res.clustree.dir <- paste0(clustree.dir, "/louvain_resolution/")
  dir.create(umaps.clustree.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(pca.clustree.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(res.clustree.dir, recursive = TRUE, showWarnings = FALSE)
  require(ggplot2)
  require(clustree)

  ## Parameters
  plot.pix <- 600
  # `%dopar%` <- foreach::"%dopar%"
  `%mydo%` <- if (BiocParallel::bpworkers(BPPARAM) > 1) foreach::"%dopar%" else foreach::"%do%"
  `%do%` <- foreach::"%do%"

  ## Builiding minimal Seurat object (for memory sake)
  miniobj <- sobj
  ### Removing other assays
  Seurat::DefaultAssay(miniobj) <- assay
  other.assays <- names(miniobj@assays)[names(miniobj@assays) != assay]
  if (length(other.assays) > 0) miniobj@assays[other.assays] <- NULL
  ### Removing scale.data slot
  miniobj@assays[[assay]]@scale.data <- matrix(nrow = 0, ncol = 0)
  ### Removing misc
  miniobj@assays[[assay]]@misc <- list()
  ### Removing other reductions
  other.reductions <- names(miniobj@reductions)[names(miniobj@reductions) != reduction]
  if (length(other.reductions) > 0) miniobj@reductions[other.reductions] <- NULL
  ### Slimming other slots
  miniobj@commands <- miniobj@misc <- list()
  require(Matrix)
  miniobj@assays[[assay]]@counts <- new("dgCMatrix")

  ## Purging putative residues from a former run
  miniobj@meta.data[, grep(colnames(miniobj@meta.data), pattern = "seurat_clusters_LE_")] <- NULL
  rm(sobj)
  gc()

  ## Reordering dimsvec if needed (as higher dims are longer to compute)
  if (which.max(dimsvec) != 1) dimsvec <- sort(dimsvec, decreasing = TRUE)

  resclust.all <- foreach::foreach(my.dims = dimsvec, .combine = "cbind", .inorder = FALSE, .errorhandling = "stop", .packages = c("Seurat", "ggplot2")) %mydo% {
    
    miniobj@graphs <- list()
    message(paste0("Dimensions 1 to ", my.dims))
    suppressMessages(miniobj <- Seurat::FindNeighbors(object = miniobj, assay = assay, dims = 1L:my.dims, reduction = reduction))

    resloop = list()
    resloop <- foreach::foreach(my.res = resvec, .inorder = FALSE, .errorhandling = "stop", .noexport = objects()) %do% {
      # for (my.res in resvec) {

      message(paste0("Testing resolution ", format(my.res, digits=2, nsmall=1, decimal.mark="."), " ..."))

      miniobj <- Seurat::FindClusters(object = miniobj, assay = assay, random.seed = my.seed, resolution = my.res, graph.name = paste0(assay, '_snn'))
      miniobj <- Seurat::RunUMAP(object = miniobj, assay = assay, dims = 1L:my.dims, reduction = reduction, seed.use = my.seed, reduction.name = paste(c(assay, reduction, my.dims, 'umap'), collapse = '_'))
      png(paste0(umaps.clustree.dir, '/', sample.name,'_uMAP_', reduction, my.dims, "_res", format(my.res, digits=2, nsmall=1, decimal.mark="."), '.png'), width = 1100, height = 1000)
      resdim.plot <- Seurat::LabelClusters(plot = Seurat::DimPlot(object = miniobj, reduction = paste(c(assay, reduction, my.dims, 'umap'), collapse = '_'), pt.size = solo.pt.size) + ggplot2::ggtitle(paste0(toupper(reduction), " dims =  ", my.dims, " ; resolution = ", format(my.res, digits=2, nsmall=1, decimal.mark="."))) + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
      print(resdim.plot)
      dev.off()

      return(list(resplot = resdim.plot, clusters = miniobj$seurat_clusters))
    }

    ## Decomposing resloop
    resplots <- unlist(resloop, recursive = FALSE)[seq.int(1, length(resloop)*2, 2)]
    resclust <- as.data.frame(unlist(resloop, recursive = FALSE)[seq.int(2, length(resloop)*2, 2)])
    resclust <- cbind(resclust, resclust)
    rm(resloop)
    colnames(resclust) <- c(paste0(paste0("seurat_clusters_LE_", reduction, my.dims, "_res", format(resvec, digits=2, nsmall=1, decimal.mark="."))), paste0(paste0("seurat_clusters_LE_res", format(resvec, digits=2, nsmall=1, decimal.mark="."), "_", reduction, my.dims)))
    miniobj@meta.data <- cbind(miniobj@meta.data, resclust)

    grid.xy <- grid.scalers(length(resplots))
    png(paste0(umaps.clustree.dir, '/', sample.name, '_uMAPs_', reduction, my.dims, '_ALLres.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    print(patchwork::wrap_plots(resplots) + patchwork::plot_layout(ncol = grid.xy[1]))
    dev.off()

    rm(resplots)


    if (length(resvec) > 1) {
      cres <- clustree::clustree(subset(miniobj@meta.data, select = grepl(paste0("seurat_clusters_LE_", reduction, my.dims, "_res"), names(miniobj@meta.data))), prefix = paste0("seurat_clusters_LE_", reduction, my.dims, "_res"))
      png(paste0(pca.clustree.dir, '/', sample.name, '_', reduction, my.dims, '.png'), width = 800, height = 1000)
      print(cres + ggplot2::ggtitle(paste0(sample.name, ', ', toupper(reduction), ' = ', my.dims)))
      dev.off()
    }

    gc(verbose = FALSE)

    return(resclust)
  }
  miniobj@meta.data[, grep(colnames(miniobj@meta.data), pattern = "seurat_clusters_LE_")] <- NULL
  miniobj@meta.data <- cbind(miniobj@meta.data, resclust.all)
  rm(resclust.all)

  for (my.res in resvec) {
    if (length(dimsvec) > 1) {
      cres <- clustree::clustree(miniobj, prefix = paste0("seurat_clusters_LE_res", format(my.res, digits=2, nsmall=1, decimal.mark="."), "_", reduction))
      png(paste0(res.clustree.dir, '/', sample.name, '_', assay, '_res', format(my.res, digits=2, nsmall=1, decimal.mark="."), '.png'), width = 800, height = 1000)
      print(cres + ggplot2::ggtitle(paste0(sample.name, ", res = ", format(my.res, digits=2, nsmall=1, decimal.mark="."))))
      dev.off()
    }
  }
  rm(miniobj)
  gc(verbose = FALSE)

}


## Louvain clustering + UMAP
louvain.cluster <- function(sobj = NULL, reduction = 'RNA_scbfa', max.dim = 100L, algorithm = 1, resolution = .8, out.dir = NULL, solo.pt.size = 2) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (!"seed" %in% names(sobj@misc$params)) stop("No seed found in Seurat object !")
  if (!(reduction %in% names(sobj@reductions))) stop(paste0('Reduction "', reduction, '" not present in the provided Seurat object !'))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  if (max.dim > ncol(sobj@reductions[[reduction]]@cell.embeddings)) stop(paste0('Max dimension requested is ', max.dim, ' whereas max dimension in "', reduction, '" reduction is ', ncol(sobj@reductions[[reduction]]@cell.embeddings)))

  ## Save command
  sobj@misc$pipeline_commands = c(sobj@misc$pipeline_commands, paste0("louvain.cluster(sobj = sobj, reduction = ", reduction, ", max.dim = ", max.dim, ", algorithm = ", algorithm, ", resolution = ", resolution, ", out.dir = ", out.dir, ", solo.pt.size = ", solo.pt.size, ")"))

  ## Restoring seed and sample name from within the object, and deleting graphs slot
  my.seed <- sobj@misc$params$seed
  sample.name <- Seurat::Project(sobj)
  sobj@graphs <- list()
  set.seed(my.seed)

  ## Find Neighbors
  sobj <- Seurat::FindNeighbors(object = sobj, assay = assay, reduction = reduction, compute.SNN = TRUE, dims = 1L:max.dim)
  names(sobj@graphs)[names(sobj@graphs) == paste0(assay, '_nn')] <- paste0(paste(c(reduction, max.dim), collapse = '.'), '_nn')
  names(sobj@graphs)[names(sobj@graphs) == paste0(assay, '_snn')] <- paste0(paste(c(reduction, max.dim), collapse = '.'))
  ori.assay <- Seurat::DefaultAssay(sobj)
  Seurat::DefaultAssay(sobj) <- assay

  ## Find clusters
  sobj <- Seurat::FindClusters(object = sobj, resolution = resolution, random.seed = my.seed, algorithm = algorithm, graph.name = paste0(paste(c(reduction, max.dim), collapse = '.')))
  Seurat::DefaultAssay(sobj) <- ori.assay

  ## Run UMAP
  sobj <- Seurat::RunUMAP(object = sobj, assay = assay, dims = 1L:max.dim, reduction = reduction, graph.name = paste0(paste(c(reduction, max.dim), collapse = '.')), reduction.name = paste(c(reduction, max.dim, 'umap'), collapse = "_"), reduction.key = tolower(paste0(reduction, max.dim, 'umap_')), seed.use = my.seed)
  sobj <- Seurat::RunUMAP(object = sobj, assay = assay, dims = 1L:max.dim, reduction = reduction, n.components = 3L, graph.name = paste0(reduction, '_snn'), reduction.name = paste(c(reduction, max.dim, 'umap3d'), collapse = "_"), reduction.key = tolower(paste0(reduction, max.dim, 'umap3d_')), seed.use = my.seed)

  if (!is.null(out.dir)) {
    png(paste0(out.dir, '/', paste0(c(sample.name, reduction, 'uMAP', 'dim'), collapse = "_"), max.dim, "_res", resolution, '.png'), width = 1000, height = 1000)
    print(Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = paste(c(reduction, max.dim, 'umap'), collapse = "_"), pt.size = solo.pt.size) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold"))
    dev.off()
    u3dlist <- list(
      Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = paste(c(reduction, max.dim, 'umap3d'), collapse = "_"), pt.size = solo.pt.size, dims = c(1,2)) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold"),
      Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = paste(c(reduction, max.dim, 'umap3d'), collapse = "_"), pt.size = solo.pt.size, dims = c(1,3)) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold"),
      Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = paste(c(reduction, max.dim, 'umap3d'), collapse = "_"), pt.size = solo.pt.size, dims = c(2,3)) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = solo.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
    )
    png(paste0(out.dir, '/', paste0(c(sample.name, reduction, 'uMAP3d', 'dim'), collapse = "_"), max.dim, "_res", resolution, '.png'), width = 2000, height = 2000)
    print(patchwork::wrap_plots(u3dlist) + patchwork::plot_layout(ncol = 2))
    dev.off()
  }

  ## Save parameters
  # ident.name <- paste0(sub(pattern = '_', replacement = '.', x = reduction), '.', max.dim, '_res.', resolution)
  ident.name <- paste0(reduction, '.', max.dim, '_res.', resolution)
  umap.name <- paste(c(reduction, max.dim, 'umap'), collapse = '_')
  umap3d.name <- paste(c(reduction, max.dim, 'umap3d'), collapse = '_')
  sobj@reductions[[umap3d.name]]@misc <- list(from.reduction = reduction,
                                              dimensions = max.dim)
  sobj@reductions[[umap.name]]@misc <- list(from.reduction = reduction,
                                            dimensions = max.dim,
                                            clustering = list(algorithm = algorithm,
                                                              resolution = resolution,
                                                              ident = ident.name))
  sobj@assays[[assay]]@misc$params$clustering[[ident.name]] = list(method = "louvain",
                                                                  algorithm = algorithm,
                                                                  resolution = resolution)
  sobj@misc$params$clustering <- c(sobj@assays[[assay]]@misc$params$clustering[[ident.name]], ident = ident.name, umap=paste(c(reduction, max.dim, 'umap'), collapse = "_"), dimensions = max.dim)


  if (!'ident2reduction' %in% names(sobj@misc)) sobj@misc$ident2reduction <- list()
  sobj@misc$ident2reduction[[ident.name]] <- umap.name
  if (!'reduction2ident' %in% names(sobj@misc)) sobj@misc$reduction2ident <- list()
  sobj@misc$reduction2ident[[umap.name]] <- ident.name
  sobj@misc$params$clustering$Rsession <- utils::capture.output(devtools::session_info())

  return(sobj)
}

## Automatic annotation of cells/clusters
### NOTE : Already tried to use @scale.data (with regression) instead of @data : IT IS NOT A GOOD IDEA,  VERY FEW CELLS SCORED!
cells.annot <- function(sobj = NULL, ident = NULL, slot = 'data', singler.setnames = NULL, clustifyr.setnames = NULL, sr.minscore = .25, cfr.minscore = .35, out.dir = NULL, solo.pt.size = 2, BPPARAM = NULL) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(ident)) stop("No ident name provided !")
  if (!ident %in% colnames(sobj@meta.data)) stop(paste0("Ident '", ident, "' not fond in Seurat object !"))
  umap.name <- sobj@misc$ident2reduction[[ident]]
  if (!umap.name %in% names(sobj@reductions)) stop(paste0("Could not find expected uMAP reduction '", umap.name, "' from ident !"))
  reduction <- sobj@reductions[[umap.name]]@misc$from.reduction
  if (!reduction %in% names(sobj@reductions)) stop(paste0("Could not find reduction '", reduction, "' from ident !"))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist in the provided Seurat object !'))
  if (!"seed" %in% names(sobj@misc$params)) stop("No seed found in Seurat object !")

  ## Save command
  sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("cells.annot(sobj = sobj, ident = ", ident, ", slot = ", slot, ", singler.setnames = ", paste0("c(", paste(singler.setnames, collapse = ","),")"), ", clustifyr.setnames = ", paste0("c(", paste(clustifyr.setnames, collapse = ","),")"), ", sr.minscore = ", sr.minscore, ", cfr.minscore = ", cfr.minscore, ", out.dir = ", out.dir, ", solo.pt.size = ", solo.pt.size, ")"))

  ## NEVER CHANGE THIS !!
  # slot <- 'data'

  ## Restoring seed and sample name from within the object
  my.seed <- sobj@misc$params$seed
  sample.name <- Seurat::Project(sobj)

  ## Building output structure
  cellannot.dir <- paste0(out.dir, "/cells_annotation/")
  sr.cellannot.dir <- paste0(cellannot.dir, "singler")
  cfr.cellannot.dir <- paste0(cellannot.dir, "clustifyr")
  dir.create(sr.cellannot.dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(cfr.cellannot.dir, recursive = TRUE, showWarnings = FALSE)

  ## Retrieving data
  mydata <- methods::slot(sobj@assays[[assay]], slot)

  ### SingleR
  suppressMessages(require(SingleR))
  suppressMessages(require(celldex))
  suppressMessages(require(scRNAseq))
  for (setname in singler.setnames) {
    print(setname)
    singler.ref <- suppressMessages(do.call(match.fun(paste0(setname)), args = list()))
    ## per cell
    entryname <- paste0("SR_", setname, "_cells")
    tryCatch( { # singler.res <- SingleR::SingleR(test = sobj@assays[[assay]]@data, ref = singler.ref, labels = singler.ref$label.main)
                #singler.res <- SingleR::SingleR(sc_data = mydata, ref_data = singler.ref@assays@data$logcounts, types = singler.ref$label.main, numCores = nthreads)
                singler.res <- SingleR::SingleR(test = mydata, ref = singler.ref, labels = singler.ref$label.main)
                singler.res$pruned.labels[singler.res@listData$tuning.scores$first < sr.minscore] <- NA
                sobj@meta.data[[entryname]] <- as.factor(unname(singler.res$pruned.labels))
                png(paste0(sr.cellannot.dir, '/', sample.name, '_', assay, '_uMAP_', entryname, '.png'), width = 1500, height = 1500)
                if (all(is.na(sobj@meta.data[[entryname]]))) {
                  suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, cols = "grey50", group.by = entryname) + ggplot2::ggtitle(paste0("SingleR predicted cell types (", setname, ") for ", ncol(sobj), " cells")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
                } else {
                  suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, group.by = entryname) + ggplot2::ggtitle(paste0("SingleR predicted cell types (", setname, ") for ", ncol(sobj), " cells")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
                }
                dev.off()
    },  error=function(error_message) { message(paste0("Error in singleR function for cells with database ", setname))} )

    ## per cluster
    entryname <- paste0("SR_", setname, "_clust")
    tryCatch( { sobj@meta.data[[entryname]] <- sobj@meta.data[[ident]]
                # singler.res <- SingleR::SingleR(test = sobj@assays[[assay]]@data, ref = singler.ref, labels = singler.ref$label.main, method = "cluster", clusters = sobj[[entryname, drop = TRUE]])
                singler.res <- SingleR::SingleR(test = mydata, ref = singler.ref, labels = singler.ref$label.main, method = "cluster", clusters = sobj[[entryname, drop = TRUE]])
                singler.res$pruned.labels[singler.res@listData$tuning.scores$first < sr.minscore] <- NA
                levels(sobj@meta.data[[entryname]]) <- singler.res$pruned.labels

                png(paste0(sr.cellannot.dir, '/', sample.name, '_', assay, '_uMAP_', entryname, '.png'), width = 1500, height = 1500)
                if (all(is.na(sobj@meta.data[[entryname]]))) {
                  suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, cols = "grey50", group.by = entryname) + ggplot2::ggtitle(paste0("SingleR predicted cell types (", setname, ") for ", length(unique(sobj$seurat_clusters)), " cluster(s)")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
                } else suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, group.by = entryname) + ggplot2::ggtitle(paste0("SingleR predicted cell types (", setname, ") for ", nlevels(sobj@meta.data[[entryname]]), " cluster(s)")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
                dev.off()
    },  error=function(error_message) { message(paste0("Error in singleR function for clusters with database ", setname))} )
  }

  ### CLUSTIFYR
  suppressPackageStartupMessages(require(clustifyrdata))
  for (setname in clustifyr.setnames) {
    print(setname)

    non.var.genes <- rownames(sobj@assays[[assay]]@data)[!(rownames(sobj@assays[[assay]]@data) %in% sobj@assays[[assay]]@var.features)]
    ## per cell
    entryname <- paste0("CFR_", setname, "_cells")
    tryCatch( { sobj@meta.data[[entryname]] <- sobj@meta.data[[ident]]
                # cfyr.res <- suppressMessages(clustifyr::clustify(input = sobj@assays[[assay]]@data, metadata = sobj[[entryname, drop = TRUE]], ref_mat = get0(setname), per_cell = TRUE, dr = umap.name, exclude_genes = non.var.genes))
                cfyr.res <- suppressMessages(clustifyr::clustify(input = mydata, metadata = sobj[[entryname, drop = TRUE]], ref_mat = get0(setname), per_cell = TRUE, dr = umap.name, exclude_genes = non.var.genes))
                cfyr.res.tbl <- clustifyr::cor_to_call(cfyr.res)
                cfyr.res.tbl <- cfyr.res.tbl[order(cfyr.res.tbl$cluster),]
                sobj@meta.data[[entryname]] <- cfyr.res.tbl$type
                sobj@meta.data[[entryname]][cfyr.res.tbl$r < cfr.minscore] <- NA

                png(paste0(cfr.cellannot.dir, '/', sample.name, '_', assay, '_uMAP_', entryname, '.png'), width = 1500, height = 1500)
                if (all(is.na(sobj@meta.data[[entryname]]))) {
                  suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, cols = "grey50", group.by = entryname) + ggplot2::ggtitle(paste0("clustifyr predicted cell types (", setname, ") for ", ncol(sobj), " cells")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
                } else {
                  suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, group.by = entryname) + ggplot2::ggtitle(paste0("clustifyr predicted cell types (", setname, ") for ", ncol(sobj), " cells")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
                }
                dev.off()
    },  error=function(error_message) { message(paste0("Error in clustifyr function for cells with database ", setname))} )

    ## per cluster
    entryname <- paste0("CFR_", setname, "_clust")
    tryCatch( { sobj@meta.data[[entryname]] <- sobj@meta.data[[ident]]
                # cfyr.res <- suppressMessages(clustifyr::clustify(input = sobj@assays[[assay]]@data, metadata = sobj[[entryname, drop = TRUE]], ref_mat = get0(setname), per_cell = FALSE, dr = umap.name, exclude_genes = non.var.genes))
                cfyr.res <- suppressMessages(clustifyr::clustify(input = mydata, metadata = sobj[[entryname, drop = TRUE]], ref_mat = get0(setname), per_cell = FALSE, dr = umap.name, exclude_genes = non.var.genes))
                cfyr.res.tbl <- clustifyr::cor_to_call(cfyr.res)
                cfyr.res.tbl <- cfyr.res.tbl[order(as.numeric(cfyr.res.tbl$cluster)),]
                levels(sobj@meta.data[[entryname]]) <- cfyr.res.tbl$type
                sobj@meta.data[[entryname]][cfyr.res.tbl$r < cfr.minscore] <- NA

                png(paste0(cfr.cellannot.dir, '/', sample.name, '_', assay, '_uMAP_', entryname, '.png'), width = 1500, height = 1500)
                if (all(is.na(sobj@meta.data[[entryname]]))) {
                  suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, cols = "grey50", group.by = entryname) + ggplot2::ggtitle(paste0("clustifyr predicted cell types (", setname, ") for ", nlevels(sobj@meta.data[[entryname]]), " cluster(s)")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
                } else {
                  suppressMessages(print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = solo.pt.size, group.by = entryname) + ggplot2::ggtitle(paste0("clustifyr predicted cell types (", setname, ") for ", nlevels(sobj@meta.data[[entryname]]), " cluster(s)")) + ggplot2::theme(legend.position = "bottom") + Seurat::DarkTheme()))
                }
                dev.off()
    },  error=function(error_message) { message(paste0("Error in clustifyr function for clusters with database ", setname))} )
  }

  ## Save parameters
  sobj@misc$params$annot <- list(ident = ident,
                                 slot = slot,
                                 singler.setnames = singler.setnames,
                                 clustifyr.setnames = clustifyr.setnames,
                                 sr.minscore = sr.minscore,
                                 cfr.minscore = cfr.minscore)
  sobj@misc$params$annot$Rsession <- utils::capture.output(devtools::session_info())
  ## Save packages versions
  sobj@misc$technical_info$SingleR <- utils::packageVersion('SingleR')
  sobj@misc$technical_info$clustifyrdata <- utils::packageVersion('clustifyrdata')
  sobj@misc$technical_info$clustifyr <- utils::packageVersion('clustifyr')

  return(sobj)
}


## Find markers
find.markers.quick <- function(sobj = NULL, ident = NULL, slot = 'data', test.use = 'wilcox', min.pct = .75, logfc.threshold = .5, only.pos = TRUE, adjp.p.max = 5E-02, topn = 10, heatmap.cols = c("gold", "blue"), out.dir = NULL) {
  if (is.null(out.dir)) stop('No output directory provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(ident)) stop("No ident name provided !")
  if (!ident %in% colnames(sobj@meta.data)) stop(paste0("Ident '", ident, "' not fond in Seurat object !"))
  umap.name <- sobj@misc$ident2reduction[[ident]]
  if (!umap.name %in% names(sobj@reductions)) stop(paste0("Could not find expected uMAP reduction '", umap.name, "' from ident !"))
  reduction <- sobj@reductions[[umap.name]]@misc$from.reduction
  if (!reduction %in% names(sobj@reductions)) stop(paste0("Could not find reduction '", reduction, "' from ident !"))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist in the provided Seurat object !'))
  if (!"seed" %in% names(sobj@misc$params)) stop("No seed found in Seurat object !")
  suppressPackageStartupMessages(require(UpSetR))
  suppressPackageStartupMessages(require(dplyr))

  #Save command
  sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("find.markers.quick(sobj = sobj, ident = ", ident, ", slot = ", slot, ", test.use = ", test.use, ", min.pct = ", min.pct, ", logfc.threshold = ", logfc.threshold, ", only.pos = ", only.pos, ", adjp.p.max = ", adjp.p.max, ", topn = ", topn, ", heatmap.cols = ", paste0("c(", paste(heatmap.cols, collapse = ","),")"), ", out.dir = ", out.dir, ")"))

  ## NEVER CHANGE THIS
  # slot <- 'data'

  ## Restoring seed and sample name from within the object
  my.seed <- sobj@misc$params$seed
  sample.name <- Seurat::Project(sobj)

  ## Scaling (for heatmap)
  if (sum(dim(sobj@assays[[assay]]@scale.data)) < 3) {
    assay.ori <- Seurat::DefaultAssay(sobj)
    Seurat::DefaultAssay(sobj) <- assay
    if(assay == 'SCT') {
      sobj <- Seurat::ScaleData(object = sobj,
                                vars.to.regress = NULL,
                                do.scale = FALSE, scale.max = Inf, block.size = 750)
    } else {
      sobj <- Seurat::ScaleData(object = sobj,
                                vars.to.regress = NULL,
                                do.scale = TRUE, scale.max = 10, block.size = 1000)
    }
    Seurat::DefaultAssay(sobj) <- assay.ori
  }

  ## Keeping track of original ident (to restore it after our computations)
  ori.ident <- Seurat::Idents(sobj)
  Seurat::Idents(sobj) <- ident

  ## Computing differentials
  fmark <- Seurat::FindAllMarkers(sobj, assay = assay, slot = slot, test.use = test.use, min.pct = min.pct, logfc.threshold = logfc.threshold, only.pos = only.pos, random.seed = my.seed)
  if (dim(fmark)[1]!=0){
    ## Filtering markers
    fmark <- fmark[fmark$p_val_adj < adjp.p.max, ]
    if (dim(fmark)[1]!=0){
      fmark <- fmark[order(fmark$cluster, fmark$p_val_adj),]

      if(slot == 'data') {
        mytop <- fmark %>% group_by(cluster) %>% top_n(n = topn, wt = avg_log2FC)
        avg_name <- 'avg_log2FC'
      } else if(slot == 'scale.data') {
        mytop <- fmark %>% group_by(cluster) %>% top_n(n = topn, wt = avg_diff)
        avg_name <- 'avg_diff'
      }

      ## Creating output dir
      fmark.dir <- paste0(out.dir, '/found_markers/')
      dir.create(fmark.dir, recursive = TRUE, showWarnings = FALSE)

      ## Heatmap
      suppressPackageStartupMessages(require(ggplot2))
      png(paste0(fmark.dir, '/', sample.name, '_findmarkers_top', topn, '_heatmap.png'), width = 1600, height = 1000)
      print(Seurat::DoHeatmap(sobj, slot = 'scale.data', features = mytop$gene, angle = 0, hjust = .5, assay = assay) + ggplot2::scale_fill_gradientn(colors = heatmap.cols))
      dev.off()

      ## Upset plots
      g.list <- sapply(unique(fmark$cluster), function(k) { return(fmark$gene[fmark$cluster == k]) }, simplify = FALSE)
      names(g.list) <- unique(fmark$cluster)
      if(length(names(g.list))>1){
        png(paste0(fmark.dir, '/', sample.name, '_findmarkers_upset_all.png'), width = 1600, height = 1000)
        up1 <- upset(fromList(g.list), sets = rev(names(g.list)), nsets = length(g.list), text.scale = 3, keep.order = TRUE)
        print(up1)
        grid::grid.text("All markers",x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
        dev.off()
      }
      g.top.list <- sapply(unique(fmark$cluster), function(k) { return(mytop$gene[mytop$cluster == k]) }, simplify = FALSE)
      names(g.top.list) <- unique(fmark$cluster)
      if(length(names(g.top.list))>1){
        png(paste0(fmark.dir, '/', sample.name, '_findmarkers_upset_top', topn, '.png'), width = 1600, height = 1000)
        up2 <- upset(fromList(g.top.list), sets = rev(names(g.top.list)), nsets = length(g.top.list), text.scale = 3, keep.order = TRUE)
        print(up2)
        grid::grid.text(paste0("Top ", topn, " markers"),x = 0.65, y=0.95, gp=grid::gpar(fontsize=30))
        dev.off()
      }

      ## Violinplots
      for(k in unique(fmark$cluster)) {
        mytop.k <- mytop[mytop$cluster == k,]
        if (nrow(mytop.k) > 0) {
          fm.list <- sapply(seq_len(nrow(mytop.k)), function(g) { Seurat::VlnPlot(sobj, assay = assay, features = mytop.k$gene[g]) + Seurat::NoLegend() + ggplot2::ggtitle(mytop.k$gene[g], subtitle = paste0(avg_name, ' = ', format(mytop.k[[avg_name]][g], digits= 3), ' ; adj.p = ', format(mytop.k$p_val_adj[g], digits = 2, scientific = TRUE))) }, simplify = FALSE)
          png(paste0(fmark.dir, '/', sample.name, '_findmarkers_top', topn, '_cluster', k, '_vln.png'), width = 1800, height = 1000)
          print(patchwork::wrap_plots(fm.list) + patchwork::plot_layout(ncol = 5))
          dev.off()
        }
      }
    }else message("No markers pass the pvalue adjusted threshold!")
  }else message("No markers identified!")

  ## Restoring original ident
  Seurat::Idents(sobj) <- ori.ident

  ## Save table
  write.table(fmark, file = paste0(fmark.dir, '/', sample.name, '_', ident, '_findmarkers_all.txt'), sep = "\t", row.names = FALSE, quote = FALSE)

  ## Cleaning
  sobj@assays[[assay]]@scale.data <- matrix(nrow = 0, ncol = 0)

  ## Save results and parameters
  sobj@misc$find.markers.quick <- fmark
  sobj@misc$params$find.markers.quick <- list(ident = ident,
                                              method = test.use,
                                              min.pct = min.pct,
                                              logfc.threshold = logfc.threshold,
                                              adjp.p.max = adjp.p.max,
                                              only.pos = only.pos,
                                              topn = topn)
  sobj@misc$params$find.markers.quick$Rsession <-  utils::capture.output(devtools::session_info())
  
  return(sobj)
}

## Technical uMAPs
technical.plot <- function(sobj = NULL, ident = NULL, out.dir = NULL, multi.pt.size = 1, gradient.cols = c("gold", "blue")) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(ident)) stop("No ident name provided !")
  if (!ident %in% colnames(sobj@meta.data)) stop(paste0("Ident '", ident, "' not found in Seurat object !"))
  umap.name <- sobj@misc$ident2reduction[[ident]]
  if (!umap.name %in% names(sobj@reductions)) stop(paste0("Could not find expected uMAP reduction '", umap.name, "' from ident !"))
  reduction <- sobj@reductions[[umap.name]]@misc$from.reduction
  if (!reduction %in% names(sobj@reductions)) stop(paste0("Could not find reduction '", reduction, "' from ident !"))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist in the provided Seurat object !'))
  require(patchwork)

  ## Creating ouput directory
  tech.dir <- paste0(out.dir, '/technical/')
  dir.create(tech.dir, recursive = TRUE, showWarnings = FALSE)

  ## Restoring sample name from within the object
  sample.name <- Seurat::Project(sobj)

  ## Setting plot parameters
  plot.num.add <- if(!is.null(ident)) 1 else 0
  plot.pix <- 600

  ## CLUSTERS
  if (!is.null(ident)) {
    # ori.ident <- Seurat::Idents(sobj)
    Seurat::Idents(sobj) <- ident
    upCLUST <- Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = multi.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
    # Seurat::Idents(sobj) <- ori.ident
  }

  ## METRICS
  metrics.plotlist <- list(
    'percent_mt' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "percent_mt", max.cutoff = .2, cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
    'percent_rb' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "percent_rb", max.cutoff = .25, cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
    'percent_st' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "percent_st", max.cutoff = .25, cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
    'log_nCount_RNA' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "log_nCount_RNA", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
    'nFeature_RNA' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "nFeature_RNA", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme()
  )

  for (p in seq_along(metrics.plotlist)) {
    png(paste0(tech.dir, '/', sample.name, '_technical_SINGLE_METRICS_', names(metrics.plotlist)[p], '_uMAP.png'), width = plot.pix, height = plot.pix)
    print(metrics.plotlist[[p]])
    dev.off()
  }
  grid.xy <- grid.scalers(length(metrics.plotlist) + plot.num.add)
  png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_METRICS_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
  if(!is.null(ident)) print(wrap_plots(c(list(upCLUST), metrics.plotlist))) else print(wrap_plots(metrics.plotlist))
  dev.off()

  ## CELL CYCLE (cyclone)
  if ('Cyclone.Phase' %in% colnames(sobj@meta.data)) {
    c_cycle.plotlist <- list(
      'Cyclone_Phase' = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "Cyclone.Phase")  + ggplot2::ggtitle("Cell Phase (cyclone)") + Seurat::DarkTheme(),
      'Cyclone_S' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Cyclone.S.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'Cyclone_G2M' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Cyclone.G2M.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'Cyclone_SmG2M' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Cyclone.SmG2M.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'Cyclone_G1' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Cyclone.G1.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme()
    )
    for (p in seq_along(c_cycle.plotlist)) {
      png(paste0(tech.dir, '/', sample.name, '_technical_SINGLE_CYCLE_', names(c_cycle.plotlist)[p], '_uMAP.png'), width = 1000, height = 1000)
      print(c_cycle.plotlist[[p]])
      dev.off()
    }
    grid.xy <- grid.scalers(length(c_cycle.plotlist) + plot.num.add)
    png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_CYCLE_Cyclone_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(wrap_plots(c(list(upCLUST), c_cycle.plotlist))) else print(wrap_plots(c_cycle.plotlist))
    dev.off()
  }else{
    c_cycle.plotlist <- list(plot_spacer())
  }

  ## CELL CYCLE (Seurat)
  if ('Seurat.Phase' %in% colnames(sobj@meta.data)) {
    s_cycle.plotlist <- list(
      'Seurat_Phase' = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "Seurat.Phase")  + ggplot2::ggtitle("Cell Phase (Seurat)") + Seurat::DarkTheme(),
      'Seurat_S' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Seurat.S.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'Seurat_G2M' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Seurat.G2M.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme(),
      'Seurat_SmG2M' = Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "Seurat.SmG2M.Score", cols = gradient.cols, order = TRUE) + Seurat::DarkTheme()
    )
    for (p in seq_along(s_cycle.plotlist)) {
      png(paste0(tech.dir, '/', sample.name, '_technical_SINGLE_CYCLE_', names(s_cycle.plotlist)[p], '_uMAP.png'), width = 1000, height = 1000)
      print(s_cycle.plotlist[[p]])
      dev.off()
    }
    grid.xy <- grid.scalers(length(s_cycle.plotlist) + plot.num.add)
    png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_CYCLE_Seurat_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(wrap_plots(c(list(upCLUST), s_cycle.plotlist))) else print(wrap_plots(s_cycle.plotlist))
    dev.off()
  }else{
    s_cycle.plotlist <- list(plot_spacer())
  }

  ## CELL CYCLE (all)
  if (('Cyclone.Phase' %in% colnames(sobj@meta.data)) & ('Seurat.Phase' %in% colnames(sobj@meta.data))) {
    grid.xy <- grid.scalers(length(c_cycle.plotlist) + length(s_cycle.plotlist) + plot.num.add)
    png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_CYCLE_ALL_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(wrap_plots(c(list(upCLUST), c_cycle.plotlist, s_cycle.plotlist))) else print(wrap_plots(c(c_cycle.plotlist, s_cycle.plotlist)))
    dev.off()
  }

  ## DOUBLETS
  doublets.plotlist <- list()
  if("scDblFinder.class" %in% colnames(sobj@meta.data)) doublets.plotlist[['ScDblFinder']] <- Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "scDblFinder.class") + ggplot2::ggtitle("Cell doublets (scDblFinder)") + Seurat::DarkTheme()
  if("hybrid_score.class" %in% colnames(sobj@meta.data)) doublets.plotlist[['scds_hybrid']] <- Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "hybrid_score.class") + ggplot2::ggtitle("Cell doublets (scds-hybrid)") + Seurat::DarkTheme()
  if("doublets_consensus.class" %in% colnames(sobj@meta.data)) doublets.plotlist[['Doublets_union']] <- Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "doublets_consensus.class") + ggplot2::ggtitle("Cell doublets (union)") + Seurat::DarkTheme()
  if("log_scran.doubletscore" %in% colnames(sobj@meta.data)) doublets.plotlist[['Doublets_scran']] <- Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = "log_scran.doubletscore", cols = gradient.cols, order = TRUE)  + ggplot2::ggtitle("Cell doublets (scran (log))") + Seurat::DarkTheme()
  if (length(doublets.plotlist) > 0) {
    for (p in seq_along(doublets.plotlist)) {
      png(paste0(tech.dir, '/', sample.name, '_technical_SINGLE_', names(doublets.plotlist)[p], '_uMAP.png'), width = 1000, height = 1000)
      print(doublets.plotlist[[p]])
      dev.off()
    }
    grid.xy <- grid.scalers(length(doublets.plotlist) + plot.num.add)
    png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_DOUBLETS_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    if(!is.null(ident)) print(wrap_plots(c(list(upCLUST), doublets.plotlist))) else print(wrap_plots(doublets.plotlist))
    dev.off()
  }else{
    doublets.plotlist <- list(plot_spacer())
  }

  ## MULTI:ALL
  grid.xy <- grid.scalers(length(metrics.plotlist) + length(doublets.plotlist) + length(c_cycle.plotlist) + length(s_cycle.plotlist) + plot.num.add)
  png(paste0(tech.dir, '/', sample.name, '_technical_MULTI_ALL_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
  if(!is.null(ident)) print(wrap_plots(c(list(upCLUST), metrics.plotlist, doublets.plotlist, c_cycle.plotlist, s_cycle.plotlist))) else print(wrap_plots(c(metrics.plotlist,doublets.plotlist,c_cycle.plotlist,s_cycle.plotlist)))
  dev.off()
}

orig.ident.plot <- function(sobj = NULL, ident = NULL, out.dir = NULL, multi.pt.size = 1) {
  if (is.null(out.dir)) stop('No output dir provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(ident)) stop("No ident name provided !")
  if (!ident %in% colnames(sobj@meta.data)) stop(paste0("Ident '", ident, "' not found in Seurat object !"))
  umap.name <- sobj@misc$ident2reduction[[ident]]
  if (!umap.name %in% names(sobj@reductions)) stop(paste0("Could not find expected uMAP reduction '", umap.name, "' from ident !"))

  ## Restoring sample name from within the object
  sample.name <- Seurat::Project(sobj)

  ## Plot
  png(paste0(out.dir, '/', sample.name, '_orig.ident_uMAP.png'), width = 1000, height = 1000)
  print(Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = "orig.ident")  + ggplot2::ggtitle("Origin identity (samples)") + Seurat::DarkTheme())
  dev.off()
}

## Markers
markers.umap.plot <- function(sobj = NULL, markers = NULL, ident = NULL, out.dir = NULL, dimplot.cols = c("gold", "blue"), multi.pt.size = 1) {
  if (is.null(out.dir)) stop('No output directory provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(ident)) stop("No ident name provided !")
  if (!ident %in% colnames(sobj@meta.data)) stop(paste0("Ident '", ident, "' not fond in Seurat object !"))
  umap.name <- sobj@misc$ident2reduction[[ident]]
  if (!umap.name %in% names(sobj@reductions)) stop(paste0("Could not find expected uMAP reduction '", umap.name, "' from ident !"))
  reduction <- sobj@reductions[[umap.name]]@misc$from.reduction
  if (!reduction %in% names(sobj@reductions)) stop(paste0("Could not find reduction '", reduction, "' from ident !"))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist in the provided Seurat object !'))
  if (is.null(markers)) stop('No marker name provided !')
  require(patchwork)

  ## Save command
  sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("markers.umap.plot(sobj = sobj, markers = ", paste0("c(", paste(markers, collapse = ","),")"), ", ident = ", ident, ", out.dir = ", out.dir, ", dimplot.cols = ", paste0("c(", paste(dimplot.cols, collapse = ","),")"), ", multi.pt.size = ", multi.pt.size, ")"))

  ## Restoring sample name from within the object and cheking markers
  sample.name <- Seurat::Project(sobj)
  markers <- markers[markers %in% rownames(sobj@assays[[assay]]@data)]

  #Creating output directory
  mark.dir <- paste0(out.dir, '/markers/')
  dir.create(mark.dir, recursive = TRUE, showWarnings = FALSE)

  ## Setting plot parameters
  plot.num.add <- if(!is.null(ident)) 1 else 0
  plot.pix <- 600

  ## CLUSTERS
  if (!is.null(ident)) {
    ori.ident <- Seurat::Idents(sobj)
    Seurat::Idents(sobj) <- ident
    upCLUST <- Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = multi.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
    Seurat::Idents(sobj) <- ori.ident
  }

  marklist <- list()
  marknames <- unique(names(markers))
  if (length(marknames) > 0) {
    for (mn in marknames) {
      mini.markers <- markers[names(markers) == mn]
      mn.plotlist <- sapply(mini.markers, function(x) {Seurat::FeaturePlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, features = x, cols = dimplot.cols) + Seurat::DarkTheme() }, simplify = FALSE)
      grid.xy <- grid.scalers(length(mn.plotlist) + plot.num.add)
      png(paste0(mark.dir, '/', sample.name, '_markers_TYPE_', gsub(pattern = " ", replacement = "_", mn), '_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
      # if(!is.null(ident)) print(upCLUST + mn.plotlist) else print(mn.plotlist)
      if(!is.null(ident)) print(patchwork::wrap_plots(append(list(upCLUST), mn.plotlist))) + patchwork::plot_layout(ncol = grid.xy[1]) else print(patchwork::wrap_plots(mn.plotlist)) + patchwork::plot_layout(ncol = grid.xy[1])
      dev.off()
      marklist <- c(marklist, list(mn.plotlist))
    }
  }
  marklist <- unlist(marklist, recursive = FALSE)

  ## SINGLES
  for(x in seq_along(marklist)) {
    png(paste0(mark.dir, '/', sample.name, '_markers_SINGLE_', markers[x], '_uMAP.png'), width = 1000, height = 1000)
    print(marklist[[x]])
    dev.off()
  }

  ## ALL
  if (length(markers) <= 15) {
    grid.xy <- grid.scalers(length(marklist) + plot.num.add)
    png(paste0(mark.dir, '/', sample.name, '_markers_ALL_uMAPs.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    # if(!is.null(ident)) print(upCLUST + patchwork::wrap_plots(marklist)) + patchwork::plot_layout(ncol = grid.xy[1]) else print(patchwork::wrap_plots(marklist))+ patchwork::plot_layout(ncol = grid.xy[1])
    if(!is.null(ident)) print(patchwork::wrap_plots(append(list(upCLUST), marklist))) + patchwork::plot_layout(ncol = grid.xy[1]) else print(patchwork::wrap_plots(marklist)) + patchwork::plot_layout(ncol = grid.xy[1])
    dev.off()
  }

  sobj@misc$markers.in <- markers

  return(sobj)
}

norm.red.plot.quick <- function(sobj = NULL, assay = 'RNA', sample.name.GE = NULL, pre.out.dir = NULL, norm.method = 'LogNormalize', features.n = 3000, normalization.vtr = NULL, dimred.method = 'pca', reduction.vtr = NULL, vtr.scale = FALSE, max.dims = 50, keep.dims = 15, keep.res = 0.8, solo.pt.size = 3, multi.pt.size = 2, file.name = NULL){
  if (is.null(pre.out.dir)) stop('No output dir provided !')
  if (!dir.exists(pre.out.dir)) stop('Output directory does not exist !')

  ### Basic normalization and dimension reduction
  sobj <- sc.normalization(sobj = sobj, assay = assay, normalization.method = norm.method, features.n = 3000, vtr = normalization.vtr)
  if(tolower(norm.method == 'sctransform')) assay <- 'SCT'
  sobj <- dimensions.reduction(sobj = sobj, reduction.method = dimred.method, assay = assay, max.dims = keep.dims, vtr = reduction.vtr, vtr.scale = vtr.scale)

  ### Building reduced normalized output dir
  norm_vtr = paste0(norm.method, if(!is.na(sobj@assays[[assay]]@misc$scaling$vtr)) paste(sobj@assays[[assay]]@misc$scaling$vtr, collapse = '_') else NULL)
  dimred_vtr = paste0(dimred.method, if(!is.na(sobj@reductions[[paste(c(assay, dimred.method), collapse = '_')]]@misc$vtr)) paste(sobj@reductions[[paste(c(assay, reduction.method), collapse = '_')]]@misc$vtr, collapse = '_') else NULL)
  norm.dim.red.dir = paste(c(pre.out.dir, norm_vtr, dimred_vtr), collapse = '/')
  dir.create(path = norm.dim.red.dir, recursive = TRUE, showWarnings = FALSE)

  ### Saving reduced normalized object
  save(sobj, file = paste0(norm.dim.red.dir, '/', paste(c(sample.name.GE, file.name, norm.method, dimred.method), collapse = '_'), '.rda'), compress = "bzip2")

  ### Building clustered output dir
  clust.dir <- paste(norm.dim.red.dir, paste(c(dimred.method, keep.dims, keep.res), collapse = '_'), sep = '/')
  dir.create(path = clust.dir, recursive = TRUE, showWarnings = FALSE)

  ## Clustering + uMAP
  sobj <- louvain.cluster(sobj = sobj, reduction = paste0(assay, "_", dimred.method), max.dim = keep.dims, resolution = keep.res, out.dir = clust.dir, solo.pt.size = solo.pt.size)

  ## Setting ident name
  ident.name <- paste0(paste0(assay, "_", dimred.method, ".", keep.dims), '_res.', stringr::str_replace(keep.res, pattern = ",", replacement = "."))

  ### Technical plots
  technical.plot(sobj = sobj, ident = ident.name, out.dir = clust.dir, multi.pt.size = multi.pt.size)

  ### Saving reduced object
  save(sobj, file = paste0(clust.dir, '/', paste(c(sample.name.GE, file.name, norm_vtr, dimred_vtr, keep.dims, keep.res), collapse = '_'), '.rda'), compress = "bzip2")
}

## Control genes
ctrl.umap.plot <- function(sobj = NULL, ctrl.genes = NULL, ident = NULL, out.dir = NULL, dimplot.cols = c("gold", "blue"), multi.pt.size = 2) {
  if (is.null(out.dir)) stop('No output directory provided !')
  if (!dir.exists(out.dir)) stop('Output directory does not exist !')
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(ident)) stop("No ident name provided !")
  if (!ident %in% colnames(sobj@meta.data)) stop(paste0("Ident '", ident, "' not fond in Seurat object !"))
  umap.name <- sobj@misc$ident2reduction[[ident]]
  if (!umap.name %in% names(sobj@reductions)) stop(paste0("Could not find expected uMAP reduction '", umap.name, "' from ident !"))
  reduction <- sobj@reductions[[umap.name]]@misc$from.reduction
  if (!reduction %in% names(sobj@reductions)) stop(paste0("Could not find reduction '", reduction, "' from ident !"))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist in the provided Seurat object !'))
  if(is.null(ctrl.genes)) stop('No control genes provided !')

  ## Save command
  sobj@misc$pipeline_commands=c(sobj@misc$pipeline_commands, paste0("ctrl.umap.plot(sobj = sobj, ctrl.genes = ",  paste0("c(", paste(ctrl.genes, collapse = ","),")"), ident, ", out.dir = ", out.dir, ", dimplot.cols = ", paste0("c(", paste(dimplot.cols, collapse = ","),")"), ", multi.pt.size = ", multi.pt.size, ")"))

  ## Restoring sample name from within the object and cheking ctrl.genes
  sample.name <- Seurat::Project(sobj)
  ctrl.genes <- unique(ctrl.genes[ctrl.genes %in% rownames(sobj@assays[[assay]]@data)])

  ## Setting plot parameters
  plot.num.add <- if(!is.null(ident)) 1 else 0
  plot.pix <- 600

  ## CLUSTERS
  if (!is.null(ident)) {
    ori.ident <- Seurat::Idents(sobj)
    Seurat::Idents(sobj) <- ident
    upCLUST <- Seurat::LabelClusters(plot = Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size) + ggplot2::ggtitle("Louvain clusters (Seurat)") + Seurat::DarkTheme(), id = "ident", size = multi.pt.size*3, repel = FALSE, color = "white", fontface = "bold")
    Seurat::Idents(sobj) <- ori.ident
  }

  ## Creating output directory
  ctrl.dir <- paste0(out.dir, '/control_genes/')
  dir.create(ctrl.dir, recursive = TRUE, showWarnings = FALSE)

  ## Plotting
  require(patchwork)
  for (x in ctrl.genes) {
    exp.umap.counts <- Seurat::FeaturePlot(object = sobj, reduction = umap.name, slot = 'counts', pt.size = multi.pt.size, features = x, cols = dimplot.cols) + ggplot2::ggtitle(paste0(x, ' expression level (counts)')) + Seurat::DarkTheme()
    exp.umap.log <- Seurat::FeaturePlot(object = sobj, reduction = umap.name, slot = 'data', pt.size = multi.pt.size, features = x, cols = dimplot.cols) + ggplot2::ggtitle(paste0(x, ' expression level (logcounts)')) + Seurat::DarkTheme()
    ctrl.list <- list(exp.umap.counts, exp.umap.log)
    if(paste0('ctrl_', assay, '_', x) %in% colnames (sobj@meta.data)) {
      tag.umap <- Seurat::DimPlot(object = sobj, reduction = umap.name, pt.size = multi.pt.size, group.by = paste0('ctrl_', assay, '_', x)) + ggplot2::ggtitle(paste0(x, ' is expressed')) + Seurat::DarkTheme()
      ctrl.list <- append(ctrl.list, list(tag.umap))
    }
    grid.xy <- grid.scalers(length(ctrl.list) + plot.num.add)
    png(paste0(ctrl.dir, '/', sample.name, '_control.genes_', x, '_uMAP.png'), width = grid.xy[1]*plot.pix, height = grid.xy[2]*plot.pix)
    # if(!is.null(ident)) print(upCLUST + ctrl.list) else print(ctrl.list)
    if(!is.null(ident)) print(patchwork::wrap_plots(append(list(upCLUST), ctrl.list))) + patchwork::plot_layout(ncol = grid.xy[1]) else print(patchwork::wrap_plots(ctrl.list)) + patchwork::plot_layout(ncol = grid.xy[1])
    dev.off()
  }

  ## Save results
  sobj@misc$ctrl.genes <- ctrl.genes

  return(sobj)
}

## Build annotated Cerebro binary
### Additional parameters (...) correspond to the cerebroApp::getMarkerGenes() function, itself forwarding parameters to the Seurat::FindAllMarkers() function
### Those recommended parameters are : only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 1E-02, test = 'wilcox'
### To filter out mito, ribo and/or stress response genes from the MostExpressed results, one can pass the corresponding RDS files to parameters mt.genes.file, crb.genes.file and/or str.genes.file, respectively.
### The 'other.reductions' filter is a boolean to exclude other reductions than the one which the 'ident.name' ident came from (with the exception of '3d' version kept).
### Same goes for the 'other.idents'
seurat2cerebro <- function(sobj = NULL, ident = NULL, clusters.colnames = NULL, sample.colname = 'orig.ident', force.assay = NULL, species = 'hg', remove.other.reductions = TRUE, remove.other.idents = TRUE, gmt.file = NULL, remove.mt.genes = FALSE, remove.crb.genes = FALSE, remove.str.genes = FALSE, file = NULL, nthreads = 1, remove.custom.DE = FALSE, only_pos_DE = FALSE, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 5E-02, test = 'wilcox') {
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(ident)) stop("No ident name provided !")
  if (!ident %in% colnames(sobj@meta.data)) stop(paste0("Ident '", ident, "' not fond in Seurat object !"))
  umap.name <- sobj@misc$ident2reduction[[ident]]
  if (!umap.name %in% names(sobj@reductions)) stop(paste0("Could not find expected uMAP reduction '", umap.name, "' from ident !"))
  reduction <- sobj@reductions[[umap.name]]@misc$from.reduction
  if (!reduction %in% names(sobj@reductions)) stop(paste0("Could not find reduction '", reduction, "' from ident !"))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist in the provided Seurat object !'))
  if (!species %in% c('hg', 'mm', 'rn')) stop('species must be "hg", "mm" or "rn" !')
  if(!sample.colname %in% colnames(sobj@meta.data)) stop(paste0("Sample column '", sample.colname, '" does not exist !'))
  if (remove.mt.genes && is.null(sobj@misc$params$QC$mito.symbols)){
    message("Warning: no mitochondrial genes symbols in seurat object! Mitochondrial genes not filtred!")
    remove.mt.genes <- FALSE
  }else if (remove.mt.genes){
    mito.symbols <- sobj@misc$params$QC$mito.symbols
  }
  if (remove.crb.genes && is.null(sobj@misc$params$QC$ribo.symbols)){
    message("Warning: no ribosomal genes symbols in seurat object! Ribosomal genes not filtred!")
    remove.crb.genes <- FALSE
  }else if (remove.crb.genes){
    ribo.symbols <- sobj@misc$params$QC$ribo.symbols
  }
  if (remove.str.genes && is.null(sobj@misc$params$QC$stress.symbols)){
    message("Warning: no stress genes symbols in seurat object! Stress genes not filtred!")
    remove.str.genes <- FALSE
  }else  if (remove.str.genes){
    stress.symbols <- sobj@misc$params$QC$stress.symbols
  }
  if (!is.null(force.assay)) {
    if (!force.assay %in% names(sobj@assays)) stop(paste0("Force assay '", force.assay, "' was not found !"))
    message(paste0("Assay forced to : '", force.assay, "'."))
  }

  ## Restriction to the provided reduction, if requested
  if(remove.other.reductions) {
    cat("\nRemoving other reductions...\n")
    keep.red <- grep(pattern = umap.name, x = names(sobj@reductions), value = TRUE)
    sobj@reductions <- sobj@reductions[keep.red]
  }
  ## Restriction to the provided ident, if requested
  if(remove.other.idents) {
    cat("\nRemoving other idents...\n")
    idents.idx <- grep(pattern = '\\_res\\.', x = colnames(sobj@meta.data))
    cur.ident.idx <- grep(pattern = paste0("^",ident), x = colnames(sobj@meta.data))
    out.idents <- idents.idx[idents.idx != cur.ident.idx]
    if(length(out.idents) > 0) sobj@meta.data <- sobj@meta.data[,-c(out.idents)]
  }
  ## Translation
  cat("\nTranslation names...\n")
  sample.name <- Seurat::Project(sobj)
  sobj@meta.data$sample <- as.factor(sobj@meta.data[[sample.colname]])
  sobj@meta.data$cluster <- if(is.null(clusters.colnames)) sobj@meta.data[[ident]] else sobj@meta.data[[clusters.colnames]]
  Seurat::Idents(sobj) <- sobj@meta.data$cluster
  sobj@meta.data$nUMI <- if(paste0('nCount_', assay) %in% colnames(sobj@meta.data)) sobj@meta.data[[paste0('nCount_', assay)]] else sobj@meta.data$nCount_RNA
  sobj@meta.data$nGene <- if(paste0('nFeature_', assay) %in% colnames(sobj@meta.data)) sobj@meta.data[[paste0('nFeature_', assay)]] else sobj@meta.data$nFeature_RNA
  if('percent_rb' %in% colnames(sobj@meta.data)) sobj@meta.data$percent_ribo <- sobj$percent_rb
  if('Seurat.Phase' %in% colnames(sobj@meta.data)) sobj@meta.data$cell_cycle_seurat <- as.factor(sobj$Seurat.Phase)
  if('Cyclone.Phase' %in% colnames(sobj@meta.data)) sobj@meta.data$cell_cycle_cyclone <- as.factor(sobj$Cyclone.Phase)
  projections_available <- names(sobj@reductions)
  if(length(projections_available) >= 1){
    projections_available_pca <- projections_available[grep(projections_available, pattern = 'pca', ignore.case = TRUE, invert = FALSE)]
    if (length(projections_available_pca) >= 1){
      for (projection in projections_available_pca){
        new_projection=stringr::str_replace(projection, "pca", "PrincipalComponentAnalysis")
        sobj@reductions[new_projection]=sobj@reductions[projection]
        sobj@reductions[projection]=NULL
      }
    }
  }
  ## Clean meta.data
  cat("\nClean meta.data...\n")
  # sobj@meta.data[[sample.colname]] <- NULL
  # if(is.null(clusters.colnames) && remove.other.idents) sobj[[ident]] <- NULL
  if(paste0('nCount_', assay) %in% colnames(sobj@meta.data)) sobj@meta.data[[paste0('nCount_', assay)]] <- NULL else sobj@meta.data$nCount_RNA <- NULL
  if(paste0('nFeature_', assay) %in% colnames(sobj@meta.data)) sobj@meta.data[[paste0('nFeature_', assay)]] <- NULL else sobj@meta.data$nFeature_RNA <- NULL
  if('percent_rb' %in% colnames(sobj@meta.data)) sobj@meta.data$percent_rb <- NULL
  if('Seurat.Phase' %in% colnames(sobj@meta.data)) sobj@meta.data$Seurat.Phase <- NULL
  if('Cyclone.Phase' %in% colnames(sobj@meta.data)) sobj@meta.data$Cyclone.Phase <- NULL
  if('pcmito_inrange' %in% colnames(sobj@meta.data)) sobj@meta.data$pcmito_inrange <- NULL
  if('pcribo_inrange' %in% colnames(sobj@meta.data)) sobj@meta.data$pcribo_inrange <- NULL
  if('min_features' %in% colnames(sobj@meta.data)) sobj@meta.data$min_features <- NULL
  if('min_counts' %in% colnames(sobj@meta.data)) sobj@meta.data$min_counts <- NULL
  if('seurat_clusters' %in% colnames(sobj@meta.data)) sobj@meta.data$seurat_clusters <- NULL
  percent_top_col = grep("percent_top_", colnames(sobj@meta.data), value = TRUE)
  if(length(percent_top_col) > 0) sobj@meta.data[percent_top_col] <- NULL

  ## Get most expressed genes
  cat("\nGet most expressed genes...\n")
  try(sobj <- cerebroApp::getMostExpressedGenes(object = sobj, assay = assay), silent = TRUE)
  if('most_expressed_genes' %in% names(sobj@misc)) {
    ## Removing MT genes (if requested)
    if (remove.mt.genes) {
      message('Removing MT genes from most expressed genes results ...')
      if('by_sample' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_sample <- sobj@misc$most_expressed_genes$by_sample[!sobj@misc$most_expressed_genes$by_sample$gene %in% mito.symbols,]
      if('by_cluster' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_cluster <- sobj@misc$most_expressed_genes$by_cluster[!sobj@misc$most_expressed_genes$by_cluster$gene %in% mito.symbols,]
    }
    ## Removing RB genes (if requested)
    if (remove.crb.genes) {
      message('Removing RB genes from most expressed genes results ...')
      if('by_sample' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_sample <- sobj@misc$most_expressed_genes$by_sample[!sobj@misc$most_expressed_genes$by_sample$gene %in% ribo.symbols,]
      if('by_cluster' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_cluster <- sobj@misc$most_expressed_genes$by_cluster[!sobj@misc$most_expressed_genes$by_cluster$gene %in% ribo.symbols,]
    }
    ## Removing STR genes (if requested)
    if (remove.str.genes) {
      message('Removing STR genes from most expressed genes results ...')
      if('by_sample' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_sample <- sobj@misc$most_expressed_genes$by_sample[!sobj@misc$most_expressed_genes$by_sample$gene %in% stress.symbols,]
      if('by_cluster' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_cluster <- sobj@misc$most_expressed_genes$by_cluster[!sobj@misc$most_expressed_genes$by_cluster$gene %in% stress.symbols,]
    }
  }

  ## Get Marker Genes
  cat("\nGet Marker Genes...\n")
  sobj <- cerebroApp::getMarkerGenes(object = sobj, assay = assay, organism = species, only_pos = only_pos, min_pct = min_pct, thresh_logFC = thresh_logFC, thresh_p_val = thresh_p_val, test = test)

  ## Checking if all samples/clusters have marker genes (to bypass a cerebroApp bug that goes on error if at least one sample has no marker genes)
  ### For samples
  if('by_sample' %in% names(sobj@misc$marker_genes)) {
    if(is.data.frame(sobj@misc$marker_genes$by_sample)) {
      samps.with.marks <- unique(sort(sobj@misc$marker_genes$by_sample$sample))
      missing.samps <- unique(sobj@meta.data$sample[!sobj@meta.data$sample %in% samps.with.marks])
      if (length(missing.samps) > 0) {
        message(paste0("No marker found for sample(s) '", paste(missing.samps, collapse = "', '"), "' ! Skipping this step ..."))
        sobj@misc$marker_genes$by_sample <- 'no_markers_found'
      }
      # for(s in missing.samps) sobj@misc$marker_genes$by_sample <- rbind(sobj@misc$marker_genes$by_sample, c(s, 'DUMMY', 1, 0, 0, 0, 1, FALSE))
      # sobj@misc$marker_genes$by_sample$sample <- as.factor(sobj@misc$marker_genes$by_sample$sample)
    }
  }
  ### For clusters
  # if('by_cluster' %in% names(sobj@misc$marker_genes)) {
  #   if(is.data.frame(sobj@misc$marker_genes$by_cluster)) {
  #     clusts.with.marks <- levels(sobj@misc$marker_genes$by_cluster$cluster)
  #     missing.clusts <- unique(sobj@meta.data[[ident]][!sobj@meta.data[[ident]] %in% clusts.with.marks])
  #     for(s in missing.clusts) sobj@misc$marker_genes$by_cluster <- rbind(sobj@misc$marker_genes$by_cluster, c(s, 'DUMMY', 1, 0, 0, 0, 1, FALSE))
  #     sobj@misc$marker_genes$by_cluster$cluster <- as.factor(sobj@misc$marker_genes$by_cluster$cluster)
  #   }
  # }

  ## Add test DEG results
  if(!remove.custom.DE){
    #get test DEG results names
    names_misc_DE = names(sobj@misc)[stringr::str_detect(names(sobj@misc),  "vs|conditions")]
    if(length(names_misc_DE)>0){
      #get test names
      test.used=unique(sort(stringr::str_replace_all(names_misc_DE, "1vsAll_|1vs1_|SvsS_|conditions_", "")))

      df_all_clust=NULL
      df_all_cond=NULL
      for (test in test.used){
        sub_names_misc_DE=names_misc_DE[stringr::str_detect(names_misc_DE, test)]

        for (i in sub_names_misc_DE){
          #get DE_type
          DE_type=stringr::str_split(i, "_")[[1]][1]
          #get data
          a=sobj@misc[[i]]
          #filtering by logFC et pval
          if (only_pos_DE) a=subset(a, logFC>=thresh_logFC & adj.P.Val<=thresh_p_val) else a=subset(a, abs(logFC)>=thresh_logFC & adj.P.Val<=thresh_p_val)
          if(dim(a)[1]==0){
            print(paste0(i, " No genes after filtering"))
            next
          }
          names(a)[names(a) == "FDR"] <- "p_val" #EdgeR
          names(a)[names(a) == "P.Value"] <- "p_val" #Limma
          #création du data.frame: NB: le slot "cluster"/"sample" correspond au test/contrast:element1,element2,min.pct,batch
          if(DE_type == "conditions"){
            df_cond = data.frame(sample=paste0(test,"_cluster_",a$selected_clusters,"_contrast_",a$contrast,"_min.pct_", a$min.pct, "_Batch_", a$batch), gene=a$genes, p_val=a$p_val, avg_logFC=a$logFC, pct.1=a$pct.1, pct.2=a$pct.2, p_val_adj=a$adj.P.Val)
            df_all_cond = rbind(df_all_cond, df_cond)
          }else{
            df_clust = data.frame(cluster=paste0(test,"_contrast_",a$tested_cluster,"-vs-", a$control_cluster,"_min.pct_", a$min.pct, "_Batch_", a$batch), gene=a$genes, p_val=a$p_val, avg_logFC=a$logFC, pct.1=a$pct.1, pct.2=a$pct.2, p_val_adj=a$adj.P.Val)
            df_all_clust = rbind(df_all_clust, df_clust)
          }
        }
      }
      ## Download bbd for "genes on cell surface" column for get markers
      if (species == 'hg' || species == 'human'){
        temp_attributes <- 'hgnc_symbol'
        temp_dataset <- 'hsapiens_gene_ensembl'
      } else if (species == 'mm' || species == 'mouse'){
        temp_attributes <- 'external_gene_name'
        temp_dataset <- 'mmusculus_gene_ensembl'
      }
      genes_on_cell_surface <- biomaRt::getBM(attributes = temp_attributes, filters = 'go', values = 'GO:0009986', mart = biomaRt::useMart('ensembl', dataset = temp_dataset))[,1]
      ## add cell surface column
      require(dplyr)
      if(!is.null(df_all_cond)) df_all_cond = df_all_cond %>% dplyr::mutate(on_cell_surface = df_all_cond$gene %in% genes_on_cell_surface)
      if(!is.null(df_all_clust)) df_all_clust = df_all_clust %>% dplyr::mutate(on_cell_surface = df_all_clust$gene %in% genes_on_cell_surface)

      #fusion avec quick.markers
      if(!is.null(df_all_cond)){
        sobj@misc$marker_genes$by_sample <- if(sobj@misc$marker_genes$by_sample %in% 'no_markers_found') df_all_cond else rbind(sobj@misc$marker_genes$by_sample, df_all_cond)
      }
      if(!is.null(df_all_clust)){
        sobj@misc$marker_genes$by_cluster <- if(sobj@misc$marker_genes$by_cluster %in% 'no_markers_found') df_all_clust else rbind(sobj@misc$marker_genes$by_cluster, df_all_clust)
      }
    }
  }

  ## Perform Enchichment by Cerebro
  cat("\nPerform Enchichment...\n")
   #save real clusters and sample name
  true_clusters <- sobj$cluster
  true_sample <- sobj$sample
   #replace by customized clusters and samples, caused by customized genes expression analysis
  if(!'no_markers_found'  %in% sobj@misc$marker_genes$by_cluster){
    false_cluster = noquote(unique(sort(sobj@misc$marker_genes$by_cluster$cluster)))
    false_cluster <- c(rep(as.character(false_cluster[1]), length(sobj$cluster)-length(false_cluster)), as.character(false_cluster))
    names(false_cluster) <- names(sobj$cluster)
    sobj$cluster <- false_cluster
  }
  if(!'no_markers_found'  %in% sobj@misc$marker_genes$by_sample){
    false_sample = noquote(unique(sort(sobj@misc$marker_genes$by_sample$sample)))
    false_sample <- c(rep(as.character(false_sample[1]), length(sobj$sample)-length(false_sample)), as.character(false_sample))
    names(false_sample) <- names(sobj$sample)
    sobj$sample <- false_sample
  }
  tryCatch( { 
    sobj <- cerebroApp::getEnrichedPathways(object = sobj)
  },
  error=function(error_message) { message("Error in cerebroApp::getEnrichedPathways function.")} )

   #restaure real clusters and sample name
  sobj$cluster <- true_clusters
  sobj$sample <- true_sample
  
  ## Perform GSEA by Cerebro
  cat("\nPerform GSEA...\n")
  if(!is.null(gmt.file)) sobj <- cerebroApp::performGeneSetEnrichmentAnalysis(object = sobj, assay = assay, GMT_file = gmt.file, parallel.sz = nthreads)

  ## Add some experiment information
  cat("\nAdd experiment information...\n")
  organism <- if (species == 'hg') "Homo sapiens" else if(species == 'mm') "Mus musculus" else if(species == 'rn') "Rattus norvegicus"
  if (!is.null(sobj@misc$params$analysis_type)) sobj@misc$parameters$analysis_type <- sobj@misc$params$analysis_type
  sobj@misc$parameters$seed <- sobj@misc$params$seed
  if (!is.null(sobj@assays$RNA@misc$params$sobj_creation)) sobj@misc$parameters$sobj_creation <- sobj@assays$RNA@misc$params$sobj_creation[1:5]
  if (!is.null(sobj@misc$params$analysis_type) && sobj@misc$params$analysis_type == "Individual analysis"){
    if(!is.null(sobj@misc$params$QC)) sobj@misc$parameters$QC <- list()
    if(!is.null(sobj@misc$params$QC$pcmito.range[1])) sobj@misc$parameters$QC$filter.droplets_percentage.mito.min <- paste0(sobj@misc$params$QC$pcmito.range[1]*100, " %")
    if(!is.null(sobj@misc$params$QC$pcmito.range[2])) sobj@misc$parameters$QC$filter.droplets_percentage.mito.max <- paste0(sobj@misc$params$QC$pcmito.range[2]*100, " %")
    if(!is.null(sobj@misc$params$QC$pcribo.range[1])) sobj@misc$parameters$QC$filter.droplets_percentage.ribo.min <- paste0(sobj@misc$params$QC$pcribo.range[1]*100, " %")
    if(!is.null(sobj@misc$params$QC$pcribo.range[2])) sobj@misc$parameters$QC$filter.droplets_percentage.ribo.max <- paste0(sobj@misc$params$QC$pcribo.range[2]*100, " %")
    if(!is.null(sobj@misc$params$QC$min.features)) sobj@misc$parameters$QC$filter.droplets_features.min <- paste0(sobj@misc$params$QC$min.features, " genes")
    if(!is.null(sobj@misc$params$QC$min.counts)) sobj@misc$parameters$QC$filter.droplets_transcript.min <- paste0(sobj@misc$params$QC$min.counts, " UMIs")
    if(!is.null(sobj@misc$params$QC$min.cells)) sobj@misc$parameters$QC$filter.genes_cell.min <- paste0(sobj@misc$params$QC$min.cells, " cells")
    if(!is.null(sobj@misc$params$doublets$method_filtering)) sobj@misc$parameters$QC$filter.doublets_method <- if (sobj@misc$params$doublets$method_filtering == "all") "scDblFinder, scds" else sobj@misc$params$doublets$method_filtering
    #sobj@misc$gene_lists$cyclone.cell.cycle.genes <- sobj@misc$params$QC$cell.cycle$cyclone.cell.cycle.genes
    #sobj@misc$gene_lists$seurat.cell.cycle.genes <- sobj@misc$params$QC$cell.cycle$seurat.cell.cycle.genes
  }
  if (!is.null(sobj@assays[[assay]]@misc$params$normalization)){
    sobj@misc$parameters$normalization <- c(list(method = sobj@assays[[assay]]@misc$params$normalization$normalization.method),
                                            sobj@assays[[assay]]@misc$params$normalization[2:3],
                                            list(High_Variable_Genes_used = sobj@assays[[assay]]@misc$params$normalization$features.used,
                                                 regressed_variables = paste0(sobj@assays[[assay]]@misc$scaling$vtr, collapse = ", ")))
  }
  if (!is.null(sobj@assays[[assay]]@misc$params$reductions$method)){
    sobj@misc$parameters$reductions <- list( method = sobj@assays[[assay]]@misc$params$reductions$method,
                                             dims.max = sobj@assays[[assay]]@misc$params$reductions$max.dims,
                                             regressed_variables = sobj@reductions[[reduction]]@misc$vtr,
                                             scale_regressed_variables = sobj@reductions[[reduction]]@misc$vtr.scale )
  }
  if (!is.null(sobj@misc$params$clustering$method)){
    sobj@misc$parameters$clustering <- list( method = sobj@misc$params$clustering$method,
                                             algorithm = sobj@misc$params$clustering$algorithm,
                                             dimensions = sobj@misc$params$clustering$dimensions,
                                             resolution = sobj@misc$params$clustering$resolution,
                                             umap = sobj@misc$params$clustering$umap,
                                             name = sobj@misc$params$clustering$ident )
  }
  if (!is.null(sobj@misc$params$annot)){
    sobj@misc$parameters$annotation <- list( clustering = sobj@misc$params$annot$ident,
                                             singleR_references = paste0(sobj@misc$params$annot$singler.setnames, collapse = ", "),
                                             clustifyr_references = paste0(sobj@misc$params$annot$clustifyr.setnames, collapse = ", "),
                                             singleR_threshold.minscore = sobj@misc$params$annot$sr.minscore,
                                             clustifyr_threshold.minscore = sobj@misc$params$annot$cfr.minscore )
  }
  if (!is.null(sobj@misc$params$find.markers.quick)){
    sobj@misc$parameters$classical_markers <- list( clustering = sobj@misc$params$find.markers.quick$ident,
                                                    method = sobj@misc$params$find.markers.quick$method,
                                                    min.pct = sobj@misc$params$find.markers.quick$min.pct,
                                                    logfs.threshold = sobj@misc$params$find.markers.quick$logfc.threshold,
                                                    adj.pval.threshold = sobj@misc$params$find.markers.quick$adjp.p.max )
  }
  if(!remove.custom.DE && length(names_misc_DE)>0){
    sobj@misc$parameters$custom_markers <- list( logfs.threshold = sobj@misc$params$find.markers.quick$logfc.threshold,
                                                 adj.pval.threshold = sobj@misc$params$find.markers.quick$adjp.p.max,
                                                 only_positive_logFC = only_pos_DE )
  }
  if(!is.null(sobj@misc$params$ADT)) sobj@misc$parameters$ADT <- sobj@misc$params$ADT
  sobj@misc$parameters$cerebro <- list( groups = paste0(groups, collapse = ", "),
                                        remove.other.reductions = remove.other.reductions,
                                        remove.other.idents = remove.other.idents,
                                        remove.mt.genes = remove.mt.genes,
                                        remove.crb.genes = remove.crb.genes,
                                        remove.str.genes = remove.str.genes,
                                        gmt.file = gmt.file )
  sobj@misc$technical_info$R <- sobj@misc$params$sobj_creation$Rsession
  sobj@misc$technical_info$cerebroApp <- utils::packageVersion('cerebroApp')
  #if (!is.null(sobj@misc$parameters$Materials_and_Methods$packages_references)) sobj@misc$parameters$Materials_and_Methods$References_packages <- sobj@misc$parameters$Materials_and_Methods$packages_references

  ## Conversion in cerebro objet
  cat("\nConversion in cerebro objet...\n")
  file = paste(c(file, if(remove.mt.genes) 'noMT' else NULL, if(remove.crb.genes) 'noRB' else NULL, if(remove.str.genes) 'noSTR' else NULL, if(!is.null(clusters.colnames)) paste0('clusterIs_', clusters.colnames) else NULL), collapse = '_')
  cerebroApp::exportFromSeurat(object = sobj, assay = assay, experiment_name = sample.name, organism = species, file = paste0(file, '.crb'), add_all_meta_data = TRUE)
  message("Cerebro file done!")
}


seurat2cerebro_1.3 <- function(sobj = NULL, ident = NULL, groups = NULL, sample.colname = 'orig.ident', force.assay = NULL, species = 'hg', remove.other.reductions = TRUE, remove.other.idents = TRUE, gmt.file = NULL, remove.mt.genes = FALSE, remove.crb.genes = FALSE, remove.str.genes = FALSE, file = NULL, nthreads = 1, remove.custom.DE = FALSE, only_pos_DE = FALSE, only_pos = TRUE, min_pct = .75, thresh_logFC = .5, thresh_p_val = 5E-02, test = 'wilcox') {
  if (is.null(sobj)) stop('No Seurat object provided !')
  if (is.null(ident)) stop("No ident name provided !")
  if (!ident %in% colnames(sobj@meta.data)) stop(paste0("Ident '", ident, "' not fond in Seurat object !"))
  umap.name <- sobj@misc$ident2reduction[[ident]]
  if (!umap.name %in% names(sobj@reductions)) stop(paste0("Could not find expected uMAP reduction '", umap.name, "' from ident !"))
  reduction <- sobj@reductions[[umap.name]]@misc$from.reduction
  if (!reduction %in% names(sobj@reductions)) stop(paste0("Could not find reduction '", reduction, "' from ident !"))
  assay <- sobj@reductions[[reduction]]@misc$from.assay
  if (!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist in the provided Seurat object !'))
  if (!species %in% c('hg', 'mm', 'rn')) stop('species must be "hg", "mm" or "rn" !')
  if(!sample.colname %in% colnames(sobj@meta.data)) stop(paste0("Sample column '", sample.colname, '" does not exist !'))
  if (remove.mt.genes && is.null(sobj@misc$params$QC$mito.symbols)){
    message("Warning: no mitochondrial genes symbols in seurat object! Mitochondrial genes not filtred!")
    remove.mt.genes <- FALSE
  }else if (remove.mt.genes){
    mito.symbols <- sobj@misc$params$QC$mito.symbols
  }
  if (remove.crb.genes && is.null(sobj@misc$params$QC$ribo.symbols)){
    message("Warning: no ribosomal genes symbols in seurat object! Ribosomal genes not filtred!")
    remove.crb.genes <- FALSE
  }else if (remove.crb.genes){
    ribo.symbols <- sobj@misc$params$QC$ribo.symbols
  }
  if (remove.str.genes && is.null(sobj@misc$params$QC$stress.symbols)){
    message("Warning: no stress genes symbols in seurat object! Stress genes not filtred!")
    remove.str.genes <- FALSE
  }else  if (remove.str.genes){
    stress.symbols <- sobj@misc$params$QC$stress.symbols
  }
  if (!is.null(force.assay)) {
    if (!force.assay %in% names(sobj@assays)) stop(paste0("Force assay '", force.assay, "' was not found !"))
    message(paste0("Assay forced to : '", force.assay, "'."))
  }
  groups = c('sample', ident, groups)
  
  ## Translation
  cat("\nTranslation names...\n")
  sample.name <- Seurat::Project(sobj)
  sobj@meta.data$sample <- as.factor(sobj@meta.data[[sample.colname]])
  Seurat::Idents(sobj) <- sobj@meta.data[[ident]]
  sobj@meta.data$nUMI <- if(paste0('nCount_', assay) %in% colnames(sobj@meta.data)) sobj@meta.data[[paste0('nCount_', assay)]] else sobj@meta.data$nCount_RNA
  sobj@meta.data$nGene <- if(paste0('nFeature_', assay) %in% colnames(sobj@meta.data)) sobj@meta.data[[paste0('nFeature_', assay)]] else sobj@meta.data$nFeature_RNA
  if('percent_rb' %in% colnames(sobj@meta.data)) sobj@meta.data$percent_ribo <- sobj$percent_rb
  projections_available <- names(sobj@reductions)
  if(length(projections_available) >= 1){
    projections_available_pca <- projections_available[grep(projections_available, pattern = 'pca', ignore.case = TRUE, invert = FALSE)]
    if (length(projections_available_pca) >= 1){
      for (projection in projections_available_pca){
        new_projection <- stringr::str_replace(projection, "pca", "PrincipalComponentAnalysis")
        sobj@reductions[new_projection]=sobj@reductions[projection]
        sobj@reductions[projection]=NULL
      }
    }
  }
  reduction <- stringr::str_replace(reduction, "pca", "PrincipalComponentAnalysis")
  umap.name <- stringr::str_replace(umap.name, "pca", "PrincipalComponentAnalysis")
  
  ## Add relationship tree based on dimension reduction
  cat("\nAdd relationship tree...\n")
  for (group in groups){
    if(length(unique(sort(sobj@meta.data[[group]]))) >1 ){
      Seurat::Idents(sobj) <- sobj@meta.data[[group]]
      sobj@reductions$pca <- sobj@reductions[[reduction]]
      dim_choosen = as.numeric(gsub("\\D", "", umap.name))
      sobj <- Seurat::BuildClusterTree(sobj, dim=1:dim_choosen, reorder = FALSE, reorder.numeric = FALSE)
      sobj@misc$trees[[group]] <- sobj@tools$'Seurat::BuildClusterTree'
      sobj@reductions$pca <- NULL
    }
  }
  sobj@tools <- list()
  Seurat::Idents(sobj) <- sobj@meta.data[[ident]]
  
  ## Clean meta.data
  cat("\nClean meta.data...\n")
  sobj@meta.data[[sample.colname]] <- NULL
  if(paste0('nCount_', assay) %in% colnames(sobj@meta.data)) sobj@meta.data[[paste0('nCount_', assay)]] <- NULL else sobj@meta.data$nCount_RNA <- NULL
  if(paste0('nFeature_', assay) %in% colnames(sobj@meta.data)) sobj@meta.data[[paste0('nFeature_', assay)]] <- NULL else sobj@meta.data$nFeature_RNA <- NULL
  if('percent_rb' %in% colnames(sobj@meta.data)) sobj@meta.data$percent_rb <- NULL
  if('pcmito_inrange' %in% colnames(sobj@meta.data)) sobj@meta.data$pcmito_inrange <- NULL
  if('pcribo_inrange' %in% colnames(sobj@meta.data)) sobj@meta.data$pcribo_inrange <- NULL
  if('min_features' %in% colnames(sobj@meta.data)) sobj@meta.data$min_features <- NULL
  if('min_counts' %in% colnames(sobj@meta.data)) sobj@meta.data$min_counts <- NULL
  if('seurat_clusters' %in% colnames(sobj@meta.data)) sobj@meta.data$seurat_clusters <- NULL
  percent_top_col = grep("percent_top_", colnames(sobj@meta.data), value = TRUE)
  if(length(percent_top_col) > 0) sobj@meta.data[percent_top_col] <- NULL
  
  ## Get most expressed genes
  cat("\nGet most expressed genes...\n")
  try(sobj <- cerebroApp::getMostExpressedGenes(object = sobj, assay = assay, groups = groups), silent = TRUE)
  if('most_expressed_genes' %in% names(sobj@misc)) {
    ## Removing MT genes (if requested)
    if (remove.mt.genes) {
      message('Removing MT genes from most expressed genes results ...')
      if('by_sample' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_sample <- sobj@misc$most_expressed_genes$by_sample[!sobj@misc$most_expressed_genes$by_sample$gene %in% mito.symbols,]
      if('by_cluster' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_cluster <- sobj@misc$most_expressed_genes$by_cluster[!sobj@misc$most_expressed_genes$by_cluster$gene %in% mito.symbols,]
    }
    ## Removing RB genes (if requested)
    if (remove.crb.genes) {
      message('Removing RB genes from most expressed genes results ...')
      if('by_sample' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_sample <- sobj@misc$most_expressed_genes$by_sample[!sobj@misc$most_expressed_genes$by_sample$gene %in% ribo.symbols,]
      if('by_cluster' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_cluster <- sobj@misc$most_expressed_genes$by_cluster[!sobj@misc$most_expressed_genes$by_cluster$gene %in% ribo.symbols,]
    }
    ## Removing STR genes (if requested)
    if (remove.str.genes) {
      message('Removing STR genes from most expressed genes results ...')
      if('by_sample' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_sample <- sobj@misc$most_expressed_genes$by_sample[!sobj@misc$most_expressed_genes$by_sample$gene %in% stress.symbols,]
      if('by_cluster' %in% names(sobj@misc$most_expressed_genes)) sobj@misc$most_expressed_genes$by_cluster <- sobj@misc$most_expressed_genes$by_cluster[!sobj@misc$most_expressed_genes$by_cluster$gene %in% stress.symbols,]
    }
  }
  
  ## Get Marker Genes
  cat("\nGet Marker Genes...\n")
  sobj <- cerebroApp::getMarkerGenes(object = sobj, assay = assay, organism = species, groups = groups, name = 'classical_markers', only_pos = only_pos, min_pct = min_pct, thresh_logFC = thresh_logFC, thresh_p_val = thresh_p_val, test = test)
  
  ## Add test DEG results
  if(!remove.custom.DE){
    #get test DEG results names
    names_misc_DE = names(sobj@misc)[stringr::str_detect(names(sobj@misc),  "vs|conditions")]
    if(length(names_misc_DE)>0){
      #get test names
      test.used=unique(sort(stringr::str_replace_all(names_misc_DE, "1vsAll_|1vs1_|SvsS_|conditions_", "")))
      
      df_all_clust=NULL
      df_all_cond=NULL
      for (test in test.used){
        sub_names_misc_DE=names_misc_DE[stringr::str_detect(names_misc_DE, test)]
        
        for (i in sub_names_misc_DE){
          #get DE_type
          DE_type=stringr::str_split(i, "_")[[1]][1]
          #get data
          a=sobj@misc[[i]]
          #filtering by logFC et pval
          if (only_pos_DE) a=subset(a, logFC>=thresh_logFC & adj.P.Val<=thresh_p_val) else a=subset(a, abs(logFC)>=thresh_logFC & adj.P.Val<=thresh_p_val)
          if(dim(a)[1]==0){
            print(paste0(i, " No genes after filtering"))
            next
          }
          names(a)[names(a) == "FDR"] <- "p_val" #EdgeR
          names(a)[names(a) == "P.Value"] <- "p_val" #Limma
          #création du data.frame: NB: le slot "cluster"/"sample" correspond au test/contrast:element1,element2,min.pct,batch
          if(DE_type == "conditions"){
            df_cond = data.frame(conditions_custom=paste0(test,"_cluster_",a$selected_clusters,"_contrast_",a$contrast,"_min.pct_", a$min.pct, "_Batch_", a$batch), gene=a$genes, p_val=a$p_val, avg_logFC=a$logFC, pct.1=a$pct.1, pct.2=a$pct.2, p_val_adj=a$adj.P.Val)
            df_all_cond = rbind(df_all_cond, df_cond)
          }else{
            df_clust = data.frame(clusters_custom=paste0(test,"_contrast_",a$tested_cluster,"-vs-", a$control_cluster,"_min.pct_", a$min.pct, "_Batch_", a$batch), gene=a$genes, p_val=a$p_val, avg_logFC=a$logFC, pct.1=a$pct.1, pct.2=a$pct.2, p_val_adj=a$adj.P.Val)
            df_all_clust = rbind(df_all_clust, df_clust)
          }
        }
      }
      ## Download bbd for "genes on cell surface" column for get markers
      if (species == 'hg'){
        temp_attributes <- 'hgnc_symbol'
        temp_dataset <- 'hsapiens_gene_ensembl'
      } else if (species == 'mm'){
        temp_attributes <- 'external_gene_name'
        temp_dataset <- 'mmusculus_gene_ensembl'
      }
      genes_on_cell_surface <- biomaRt::getBM(attributes = temp_attributes, filters = 'go', values = 'GO:0009986', mart = biomaRt::useMart('ensembl', dataset = temp_dataset))[,1]
      ## add cell surface column
      require(dplyr)
      if(!is.null(df_all_cond)) df_all_cond = df_all_cond %>% dplyr::mutate(on_cell_surface = df_all_cond$gene %in% genes_on_cell_surface)
      if(!is.null(df_all_clust)) df_all_clust = df_all_clust %>% dplyr::mutate(on_cell_surface = df_all_clust$gene %in% genes_on_cell_surface)
      
      #fusion avec quick.markers
      if(!is.null(df_all_cond)){
        sobj@misc$marker_genes[["differential_custom"]]$conditions_custom <- if(is.null(sobj@misc$marker_genes[["differential_custom"]]$conditions_custom)) df_all_cond else rbind(sobj@misc$marker_genes[["differential_custom"]]$conditions_custom, df_all_cond)
      }
      if(!is.null(df_all_clust)){
        sobj@misc$marker_genes[["differential_custom"]]$clusters_custom <- if(is.null(sobj@misc$marker_genes[["differential_custom"]]$clusters_custom)) df_all_clust else rbind(sobj@misc$marker_genes[["differential_custom"]]$clusters_custom, df_all_clust)
      }
    }
  }
  
  ## Perform Enchichment by Cerebro
  cat("\nPerform Enchichment...\n")
  tryCatch( { sobj <- cerebroApp::getEnrichedPathways(object = sobj, marker_genes_input = 'classical_markers') },
            error=function(error_message) { message("Error in cerebroApp::getEnrichedPathways function.")} )
  if(!is.null(sobj@misc$marker_genes[["differential_custom"]])) tryCatch( { sobj <- cerebroApp::getEnrichedPathways(object = sobj, marker_genes_input = 'differential_custom') },
            error=function(error_message) { message("Error in cerebroApp::getEnrichedPathways function.")} )
  
  ## Perform GSEA by Cerebro
  cat("\nPerform GSEA...\n")
  if(!is.null(gmt.file)) for (group in groups){ try(sobj <- cerebroApp::performGeneSetEnrichmentAnalysis(object = sobj, assay = assay, groups = group, name = 'GSVA', GMT_file = gmt.file, parallel.sz = nthreads)) }
  
  ## Add some experiment information
  cat("\nAdd experiment information...\n")
  organism <- if (species == 'hg') "Homo sapiens" else if(species == 'mm') "Mus musculus" else if(species == 'rn') "Rattus norvegicus"
  if (!is.null(sobj@misc$params$analysis_type)) sobj@misc$parameters$analysis_type <- sobj@misc$params$analysis_type
  sobj@misc$parameters$seed <- sobj@misc$params$seed
  if (!is.null(sobj@assays$RNA@misc$params$sobj_creation)) sobj@misc$parameters$sobj_creation <- sobj@assays$RNA@misc$params$sobj_creation[1:5]
  if (!is.null(sobj@misc$params$analysis_type) && sobj@misc$params$analysis_type == "Individual analysis"){
    if(!is.null(sobj@misc$params$QC)) sobj@misc$parameters$QC <- list()
    if(!is.null(sobj@misc$params$QC$pcmito.range[1])) sobj@misc$parameters$QC$filter.droplets_percentage.mito.min <- paste0(sobj@misc$params$QC$pcmito.range[1]*100, " %")
    if(!is.null(sobj@misc$params$QC$pcmito.range[2])) sobj@misc$parameters$QC$filter.droplets_percentage.mito.max <- paste0(sobj@misc$params$QC$pcmito.range[2]*100, " %")
    if(!is.null(sobj@misc$params$QC$pcribo.range[1])) sobj@misc$parameters$QC$filter.droplets_percentage.ribo.min <- paste0(sobj@misc$params$QC$pcribo.range[1]*100, " %")
    if(!is.null(sobj@misc$params$QC$pcribo.range[2])) sobj@misc$parameters$QC$filter.droplets_percentage.ribo.max <- paste0(sobj@misc$params$QC$pcribo.range[2]*100, " %")
    if(!is.null(sobj@misc$params$QC$min.features)) sobj@misc$parameters$QC$filter.droplets_features.min <- paste0(sobj@misc$params$QC$min.features, " genes")
    if(!is.null(sobj@misc$params$QC$min.counts)) sobj@misc$parameters$QC$filter.droplets_transcript.min <- paste0(sobj@misc$params$QC$min.counts, " UMIs")
    if(!is.null(sobj@misc$params$QC$min.cells)) sobj@misc$parameters$QC$filter.genes_cell.min <- paste0(sobj@misc$params$QC$min.cells, " cells")
    if(!is.null(sobj@misc$params$doublets$method_filtering)) sobj@misc$parameters$QC$filter.doublets_method <- if (sobj@misc$params$doublets$method_filtering == "all") "scDblFinder, scds" else sobj@misc$params$doublets$method_filtering
    #sobj@misc$gene_lists$cyclone.cell.cycle.genes <- sobj@misc$params$QC$cell.cycle$cyclone.cell.cycle.genes
    #sobj@misc$gene_lists$seurat.cell.cycle.genes <- sobj@misc$params$QC$cell.cycle$seurat.cell.cycle.genes
  }
  if (!is.null(sobj@assays[[assay]]@misc$params$normalization)){
    sobj@misc$parameters$normalization <- c(list(method = sobj@assays[[assay]]@misc$params$normalization$normalization.method),
                                            sobj@assays[[assay]]@misc$params$normalization[2:3],
                                            list(High_Variable_Genes_used = sobj@assays[[assay]]@misc$params$normalization$features.used,
                                                 regressed_variables = paste0(sobj@assays[[assay]]@misc$scaling$vtr, collapse = ", ")))
  }
  if (!is.null(sobj@assays[[assay]]@misc$params$reductions$method)){
    sobj@misc$parameters$reductions <- list( method = sobj@assays[[assay]]@misc$params$reductions$method,
                                       dims.max = sobj@assays[[assay]]@misc$params$reductions$max.dims,
                                       regressed_variables = sobj@reductions[[reduction]]@misc$vtr,
                                       scale_regressed_variables = sobj@reductions[[reduction]]@misc$vtr.scale )
  }
  if (!is.null(sobj@misc$params$clustering$method)){
    sobj@misc$parameters$clustering <- list( method = sobj@misc$params$clustering$method,
                                          algorithm = sobj@misc$params$clustering$algorithm,
                                          dimensions = sobj@misc$params$clustering$dimensions,
                                          resolution = sobj@misc$params$clustering$resolution,
                                          umap = sobj@misc$params$clustering$umap,
                                          name = sobj@misc$params$clustering$ident )
  }
  if (!is.null(sobj@misc$params$annot)){
    sobj@misc$parameters$annotation <- list( clustering = sobj@misc$params$annot$ident,
                                           singleR_references = paste0(sobj@misc$params$annot$singler.setnames, collapse = ", "),
                                           clustifyr_references = paste0(sobj@misc$params$annot$clustifyr.setnames, collapse = ", "),
                                           singleR_threshold.minscore = sobj@misc$params$annot$sr.minscore,
                                           clustifyr_threshold.minscore = sobj@misc$params$annot$cfr.minscore )
  }
  if (!is.null(sobj@misc$params$find.markers.quick)){
    sobj@misc$parameters$classical_markers <- list( clustering = sobj@misc$params$find.markers.quick$ident,
                                                  method = sobj@misc$params$find.markers.quick$method,
                                                  min.pct = sobj@misc$params$find.markers.quick$min.pct,
                                                  logfs.threshold = sobj@misc$params$find.markers.quick$logfc.threshold,
                                                  adj.pval.threshold = sobj@misc$params$find.markers.quick$adjp.p.max )
  }
  if(!remove.custom.DE && length(names_misc_DE)>0){
    sobj@misc$parameters$custom_markers <- list( logfs.threshold = sobj@misc$params$find.markers.quick$logfc.threshold,
                                                 adj.pval.threshold = sobj@misc$params$find.markers.quick$adjp.p.max,
                                                 only_positive_logFC = only_pos_DE )
  }
  if(!is.null(sobj@misc$params$ADT)) sobj@misc$parameters$ADT <- sobj@misc$params$ADT
  sobj@misc$parameters$cerebro <- list( groups = paste0(groups, collapse = ", "),
                                        remove.other.reductions = remove.other.reductions,
                                        remove.other.idents = remove.other.idents,
                                        remove.mt.genes = remove.mt.genes,
                                        remove.crb.genes = remove.crb.genes,
                                        remove.str.genes = remove.str.genes,
                                        gmt.file = gmt.file )
  sobj@misc$technical_info$R <- sobj@misc$params$sobj_creation$Rsession
  sobj@misc$technical_info$cerebroApp <- utils::packageVersion('cerebroApp')
  #if (!is.null(sobj@misc$parameters$Materials_and_Methods$packages_references)) sobj@misc$parameters$Materials_and_Methods$References_packages <- sobj@misc$parameters$Materials_and_Methods$packages_references

  ## Restriction to the provided reduction, if requested
  if(remove.other.reductions) {
    cat("\nRemoving other reductions ...\n")
    keep.red <- grep(pattern = umap.name, x = names(sobj@reductions), value = TRUE)
    sobj@reductions <- sobj@reductions[keep.red]
  }

  ## Restriction to the provided ident, if requested
  if(remove.other.idents) {
    cat("\nRemoving other idents ...\n")
    idents.idx <- grep(pattern = '\\_res\\.', x = colnames(sobj@meta.data))
    cur.ident.idx <- grep(pattern = paste0("^",ident), x = colnames(sobj@meta.data))
    out.idents <- idents.idx[idents.idx != cur.ident.idx]
    if(length(out.idents) > 0) sobj@meta.data <- sobj@meta.data[,-c(out.idents)]
  }
  
  ## Conversion in cerebro objet
  cat("\nConversion in cerebro objet...\n")
  file = paste(c(file, if(remove.mt.genes) 'noMT' else NULL, if(remove.crb.genes) 'noRB' else NULL, if(remove.str.genes) 'noSTR' else NULL), collapse = '_')
  cerebroApp::exportFromSeurat(object = sobj, assay = assay, groups = groups, cell_cycle = c("Seurat.Phase","Cyclone.Phase"), experiment_name = sample.name, organism = organism, file = paste0(file, '.crb'), add_all_meta_data = TRUE)
  message("Cerebro file done!")
}


## Computing correlations between two assays with corresponding features (ex: RNA and ADT)
### assay features must be in the same order !
feature.cor <- function(sobj = NULL, assay1 = 'RNA', assay2 = 'ADT', assay1.features = NULL, assay2.features = NULL, slot = 'data', cor.method = 'spearman', zero.filter = TRUE, gene.names = NULL) {
  #selecting cells with no zero expression
  cell.idx.list <- sapply(
    seq_along(gene.names),
    function(k) {
      if (zero.filter) zerocells.get(sobj = sobj, assay1 = assay1, assay2 = assay2, assay1.feature = assay1.features[k], assay2.feature = assay2.features[k], slot = slot) else rep(TRUE, ncol(sobj@assays[[assay1]]))
    },
    simplify = FALSE
  )
  #correlation
  corvec <- vapply(
    seq_along(gene.names),
    function(k) {
      cell.idx <- cell.idx.list[[k]]
      if(length(which(cell.idx)) < 2) return(NA) # must to get at least 2 cells with RNA and protein expressions to compute correlation
      return(
        tryCatch({
          cor.test(slot(sobj@assays[[assay1]], slot)[rownames(slot(sobj@assays[[assay1]], slot)) == assay1.features[k],cell.idx],
                   slot(sobj@assays[[assay2]], slot)[rownames(slot(sobj@assays[[assay2]], slot)) == assay2.features[k],cell.idx],
                   method = cor.method)$estimate
          }, error=function(cond) {
             return(NA)
          })
      )
    },
    .1)
  out.df <- if (zero.filter) {
    data.frame(vapply(cell.idx.list, function(x) { length(which(x)) }, 1L),  corvec, stringsAsFactors = FALSE)
  } else {
    data.frame(corvec, stringsAsFactors = FALSE)
  }
  colnames(out.df) <- if (zero.filter) {
    c(paste(c(slot, 'non0'), collapse = '_'), paste(c('cor', slot, cor.method, '0filt'), collapse = '_'))
  } else {
    c(paste(c('cor', slot, cor.method), collapse = '_'))
  }
  return(out.df)
  # if (zero.filter) {
  #   out.df <- data.frame(vapply(cell.idx.list, function(x) { length(which(x)) }, 1L),  corvec, stringsAsFactors = FALSE)
  #   colnames(out.df) <- c(paste(c('cor', slot, 'non0'), collapse = '_'), paste(c('cor', slot, cor.method, '0filt'), collapse = '_'))
  # } else {
  #   out.df <-data.frame(corvec, stringsAsFactors = FALSE)
  #   colnames(out.df) <- c(paste(c('cor', slot, cor.method, '0filt'), collapse = '_'))
  # }
}

## Make FeaturePlot for a list of genes with exception capture
feature_plots <- function(sobj, assay, features, slot, reduction, min.cutoff, max.cutoff){
  if(!(assay %in% names(sobj@assays))) stop(paste0('Assay "', assay, '" does not exist !'))
  Seurat::DefaultAssay(sobj) <- assay

  plots=patchwork::wrap_plots(lapply(seq_len(length(features)), function(k){
    tryCatch({
      Seurat::FeaturePlot(sobj, features = features[k], slot = slot, reduction = reduction, ncol = 1, pt.size = 2, order = TRUE, min.cutoff = min.cutoff[k], max.cutoff = max.cutoff[k])
    },
    warning=function(cond) {
      if(grepl("All cells have the same value", as.character(cond)) ) tmp_plot = Seurat::FeaturePlot(sobj, features = features[k], slot = slot, reduction = reduction, ncol = 1, pt.size = 2, order = TRUE, min.cutoff = min.cutoff[k], max.cutoff = max.cutoff[k])
      if(grepl("The following requested variables were not found", as.character(cond)) ) tmp_plot = patchwork::wrap_elements(patchwork::plot_spacer() + patchwork::plot_annotation(title = features[k], theme = ggplot2::theme(plot.title = ggplot2::element_text(size=18, hjust=0.5, face="bold"))))
      if(grepl("Could not find", as.character(cond)) ) tmp_plot = patchwork::wrap_elements(patchwork::plot_spacer() + patchwork::plot_annotation(title = features[k], theme = ggplot2::theme(plot.title = ggplot2::element_text(size=18, hjust=0.5, face="bold"))))
      return(tmp_plot)
    })
  }),
  ncol=1) + plot_annotation(title = assay, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=50, hjust=0.5, face="bold")))

  return(plots)
}

## Get cells with non-null expression in assay1 and assay2, with corresponding features between assays
zerocells.get <- function(sobj = NULL, assay1 = 'RNA', assay2 = 'ADT', assay1.feature = NULL, assay2.feature = NULL, slot = 'data') {
  cell.idx1 <- slot(sobj@assays[[assay1]], slot)[rownames(slot(sobj@assays[[assay1]], slot)) == assay1.feature,] > 0
  cell.idx2 <- slot(sobj@assays[[assay2]], slot)[rownames(slot(sobj@assays[[assay2]], slot)) == assay2.feature,] > 0
  if(length(cell.idx1)==0) cell.idx1 = rep(FALSE, ncol(sobj@assays[[assay1]]))
  if(length(cell.idx2)==0) cell.idx2 = rep(FALSE, ncol(sobj@assays[[assay2]]))
  return(cell.idx1 & cell.idx2)
}

## Get the list of assays of a Seurat object
get.assays <- function(sobj = NULL) {
  if (is.null(sobj)) stop('No Seurat object provided !')
  print(Seurat::Assays(sobj))
}

## Get the list of reductions of a Seurat object
get.reductions <- function(sobj = NULL) {
  if (is.null(sobj)) stop('No Seurat object provided !')
  print(Seurat::Reductions(sobj))
}

## Get the list of idents that can be used out of a Seurat object
get.idents <- function(sobj = NULL, pattern = '_res') {
  if (is.null(sobj)) stop('No Seurat object provided !')
  idents.list <- grep(pattern = pattern, x = colnames(sobj@meta.data))
  for (id in idents.list) {
    message(paste0("\n", colnames(sobj@meta.data)[id]))
    print(table(sobj@meta.data[,id]))
  }
}

## Get the list of graphs of a Seurat object
get.graphs <- function(sobj = NULL) {
  if (is.null(sobj)) stop('No Seurat object provided !')
  # print(grep("_res", x = colnames(sobj@meta.data), value = T))
  print(names(sobj@graphs))
}

## Function to get the multiplot scalers
grid.scalers <- function(n = 1) {
  x <- ceiling(sqrt(n))
  y <- ceiling(n / x)
  # print(paste0("X : ", x))
  # print(paste0("Y : ", y))
  return(c(x, y))
}
# save statistics from sobj (nb cells, genes, mito etc.)
save_stat <- function(sobj = NULL, sample.names = NULL, title = NULL, out.dir = NULL){
  if(length(sobj) != length(sample.names)) stop(paste0("sample.names does not correspond to sobj provided! (not the same length: ", length(sobj), "sobj  for ", length(sample.names), "name(s) )"))
  if(class(sobj) != "list") sobj.list = list(sobj) else sobj.list = sobj
  rm(sobj)
  gc()
  names(sobj.list) = sample.names
  stat_tot=data.frame(matrix(vector(), 24, 0),stringsAsFactors=F)
  for (i in sample.names){
    stat=format(as.data.frame(unlist(sobj.list[[i]]@misc$excel)), trim = TRUE, drop0trailing = TRUE, decimal.mark = ",", scientific = FALSE)
    names_row_stat=c("Droplet_Quality.captured_droplet", "Droplet_Quality.total_number_UMI", "Droplet_Quality.estimated_cells",
                     "Droplet_Quality.estimated_UMI", "Droplet_Quality.fraction_read_in_cells", "Cells_Quality.mito_summary.Median",
                     "Cells_Quality.ribo_summary.Median", "Cells_Quality.stress_summary.Median", "Cells_Quality.filter_cells_genes",
                     "Cells_Quality.filter_cells_counts", "Cells_Quality.genes_per_cell_summary.Median", "Cells_Quality.UMI_per_cell_summary.Median",
                     "After_QC_cells_filtering.estim_cells_G1","After_QC_cells_filtering.estim_cells_G2M", "After_QC_cells_filtering.estim_cells_S",
                     "After_QC_cells_filtering.Genes_covered", "After_QC_cells_feature_filtering.estim_doublets","Final_Cells_Quality.mito_summary.Median",
                     "Final_Cells_Quality.ribo_summary.Median","Final_Cells_Quality.stress_summary.Median","Final_Cells_Quality.genes_per_cell_summary.Median",
                     "Final_Cells_Quality.UMI_per_cell_summary.Median","Final_Cells_Quality.nb_genes", "Final_Cells_Quality.nb_cells")
    stat=as.data.frame(stat[names_row_stat,], row.names = names_row_stat)
    colnames(stat)=i
    stat_tot=cbind(stat_tot, stat)
  }
  write.table(stat_tot, file = paste0(out.dir, title, '_stat.txt'), sep = ";", row.names = TRUE, col.names = TRUE, quote = FALSE)
}

## Loading TCR data and filter barcode according to sobj
load.sc.tcr.bcr <- function(x=1, sobj=NULL, vdj.input.file, sample.name=NULL){
  message(paste0("Loading '", vdj.input.file[x], "' ..."))
  df <- read.table(file = vdj.input.file[x], sep = ',', header = TRUE, stringsAsFactors = FALSE)
  df$barcode <- gsub(pattern = '-1', replacement = '', x = df$barcode)
  if (!is.null(sample.name)) df$barcode=paste0(sample.name[x],"_", df$barcode)
  df <- df[df$barcode %in% colnames(sobj),]
  if (!is.null(sample.name)) df$barcode=gsub(pattern =paste0(sample.name[x],"_"), replacement = '', x = df$barcode)
  return(df)
}

## Basic QC metrics for TCR/BCR
QC.tcr.bcr <- function(cr_res=NULL, out.dir=global_output, type = 'TCR'){
  require(patchwork)
  require(ggplot2)
  require(dplyr)
  ## Quantification analysis
  cr_res$productive[cr_res$productive == "True"] <- "Productive"
  cr_res$productive[cr_res$productive == "False"] <- "No productive"
  ### Quatification total Productive TCR: (sometimes 2 TCR for 1 cell)
  plot_productive = ggplot(cr_res, aes(x=chain, fill=productive)) +
    geom_bar(stat="count", position=position_dodge()) +
    xlab("") + ylab("") + ggtitle("Productivity") +
    labs(caption = paste("(Don't forget: sometimes there are more than 1 ", type," by cell)"))
  ### Quantification nb TCR/cell/productivity
  table_receptor_number = cr_res %>% select("barcode","chain", "productive") %>% group_by(productive, barcode, chain) %>% summarise(nb = n()) %>% as.data.frame()
  table_receptor_number$nb = as.character(table_receptor_number$nb)
  plot_receptor_number = ggplot(table_receptor_number, aes(x=chain, fill=nb)) +
    geom_bar(stat="count", position=position_dodge2(preserve = "single")) + facet_grid(. ~ productive) +
    xlab("") + ylab("log10(counts)") + ggtitle(type," quantification by cell") + scale_y_log10()
  ### Save
  nb_cell_without_receptor=length(colnames(sobj)[!(colnames(sobj) %in% cr_res$barcode)]) #nb of cells without TCR/BCR
  png(paste0(out.dir,'/QC_quantif.png'), width = 1200, height = 400)
  print(((plot_productive | plot_receptor_number ) +
      plot_annotation(title = sample.name, subtitle = paste0("(",nb_cell_without_receptor, " cells without ", type," on ",dim(sobj)[2]," cells)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) +
      plot_layout(width = c(1, 2))))
  dev.off()
}

## Quantification of global unique contig analysis
Quantif.unique.g <- function(combined = NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL){
  require(patchwork)
  require(dplyr)
  for(x in list_type_clT) {
    ### Tables
    scR_res = scRepertoire::quantContig(combined, cloneCall=x, scale = T, exportTable = T)
    tmp = scR_res %>% select("total_contigs_nb"="total","unique_contig_nb"="contigs","unique_contig_pct"="scaled") %>% round(2)
    assign(paste0("table_quantContig_",sub("\\+","_",x)),data.frame(lapply(tmp, as.character), row.names = scR_res$values))
    ### Plots
    assign(paste0("plot_quantUniqueContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::quantContig(combined, cloneCall=x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold"))) + Seurat::NoLegend()))
  }
  ### Save
  png(paste0(out.dir,'/quantUniqueContig.png'), width = 1800, height = 600)
  print((plot_quantUniqueContig_gene_nt | plot_quantUniqueContig_gene | plot_quantUniqueContig_nt | plot_quantUniqueContig_aa ) / (gridExtra::grid.arrange( gridExtra::tableGrob(table_quantContig_gene_nt, theme = gridExtra::ttheme_default(base_size = 10)) , gridExtra::tableGrob(table_quantContig_gene, theme = gridExtra::ttheme_default(base_size = 10)), gridExtra::tableGrob(table_quantContig_nt, theme = gridExtra::ttheme_default(base_size = 10)), gridExtra::tableGrob(table_quantContig_aa, theme = gridExtra::ttheme_default(base_size = 10)), nrow=1) ) +
      plot_annotation(title = sample.name, subtitle = paste0("(",dim(sobj)[2]," cells)"), caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) +
      plot_layout(heights = c(2, 1)))
  dev.off()
}

## Global clonal Homeostasis analysis
Homeo.g <-function(combined = NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL){
  require(patchwork)
  require(dplyr)
  for(x in list_type_clT) {
    ### Tables
    scR_res = scRepertoire::clonalHomeostasis(combined, cloneCall = x, exportTable = TRUE) %>% as.data.frame() %>% round(4)
    colnames(scR_res)=c("Rare","Small","Medium", "Large", "Hyperexpanded")
    assign(paste0("table_clhomeo_",sub("\\+","_",x)),data.frame(lapply(scR_res, as.character), row.names = rownames(scR_res)))
    ### Plots
    assign(paste0("plot_clhomeo_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::clonalHomeostasis(combined, cloneCall = x) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.65, face="bold"))) ))
  }
  ### Save
  png(paste0(out.dir,'/clhomeo.png'), width = 2100, height = 800)
  print( (plot_clhomeo_gene_nt | plot_clhomeo_gene | plot_clhomeo_nt | plot_clhomeo_aa )  / (gridExtra::grid.arrange( gridExtra::tableGrob(table_clhomeo_gene_nt, theme = gridExtra::ttheme_default(base_size = 10)),gridExtra::tableGrob(table_clhomeo_gene, theme = gridExtra::ttheme_default(base_size = 10)),gridExtra::tableGrob(table_clhomeo_nt, theme = gridExtra::ttheme_default(base_size = 10)),gridExtra::tableGrob(table_clhomeo_aa, theme = gridExtra::ttheme_default(base_size = 10)), nrow=1) ) +
           plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) +
           plot_layout(heights = c(2, 1)))
  dev.off()
}

## Global Clonal Proportions analysis
Prop.g <-function(combined = NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL){
  require(patchwork)
  require(dplyr)
  for(x in list_type_clT) {
    ### Tables
    scR_res = scRepertoire::clonalProportion(combined, cloneCall = x, exportTable = TRUE) %>% as.data.frame() %>% round(4)
    tmp=data.frame(lapply(scR_res, as.character), row.names = rownames(scR_res))
    colnames(tmp)=colnames(scR_res)
    assign(paste0("table_clprop_",sub("\\+","_",x)),tmp)
    ### Plots
    assign(paste0("plot_clprop_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::clonalProportion(combined, cloneCall = x) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.65, face="bold"))) ))
  }
  ### Save
  png(paste0(out.dir,'/clprop.png'), width = 2100, height = 800)
  print( (plot_clprop_gene_nt | plot_clprop_gene | plot_clprop_nt | plot_clprop_aa )  / (gridExtra::grid.arrange( gridExtra::tableGrob(table_clprop_gene_nt, theme = gridExtra::ttheme_default(base_size = 10)),gridExtra::tableGrob(table_clprop_gene, theme = gridExtra::ttheme_default(base_size = 10)),gridExtra::tableGrob(table_clprop_nt, theme = gridExtra::ttheme_default(base_size = 10)),gridExtra::tableGrob(table_clprop_aa, theme = gridExtra::ttheme_default(base_size = 10)), nrow=1) ) +
           plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) +
           plot_layout(heights = c(2, 1)))
  dev.off()
}

## Diversity analysis
Div.g <-function(combined = NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL){
  require(patchwork)
  require(dplyr)
  for(x in list_type_clT) {
    ### Tables
    scR_res = scRepertoire::clonalDiversity(combined, cloneCall = x, group = "samples", exportTable = TRUE) %>% select("Shannon", "Inv.Simpson", "Chao", "ACE") %>% round(2)
    assign(paste0("table_cldiv_",sub("\\+","_",x)),data.frame(lapply(scR_res, as.character), row.names = rownames(scR_res)))
    ### Plots
    assign(paste0("plot_cldiv_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::clonalDiversity(combined, cloneCall = x, group = "samples") + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.5, face="bold"))) ))
  }
  ### Save
  png(paste0(out.dir,'/cldiv.png'), width = 2100, height = 800)
  print( (plot_cldiv_gene_nt | plot_cldiv_gene | plot_cldiv_nt | plot_cldiv_aa )  / (gridExtra::grid.arrange( gridExtra::tableGrob(table_cldiv_gene_nt, theme = gridExtra::ttheme_default(base_size = 10)),gridExtra::tableGrob(table_cldiv_gene, theme = gridExtra::ttheme_default(base_size = 10)),gridExtra::tableGrob(table_cldiv_nt, theme = gridExtra::ttheme_default(base_size = 10)),gridExtra::tableGrob(table_cldiv_aa, theme = gridExtra::ttheme_default(base_size = 10)), nrow=1) ) +
           plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) +
           plot_layout(heights = c(2, 1)))
  dev.off()
}

### Spliting CTstrict (into separate columns for TRA-V/J/C, TRB-V/J/C and corresponding sequences, with 2 possible clonotypes) and save as metadata
### and Adding length of TR sequence to meta.data
split.CTstrict.tcr <- function(sobj=NULL){
  ctsplit <- t(vapply(sobj@meta.data$CTstrict, function(x) {
    outvec <- rep('NA', 18)
    if (!is.na(x)) {
      pass1 <- unlist(strsplit(x = x, split = "\\_"))
      pass2a <- unlist(strsplit(x = pass1[1], split = "(\\;|\\.)"))
      if (length(pass2a) %% 3 == 0) {
        outvec[1:length(pass2a)] <- pass2a
      } else if (length(pass2a) > 1) {
        outvec[4:6] <- pass2a[2:4]
      }
      pass2b <- unlist(strsplit(x = pass1[2], split = "\\;"))
      outvec[7:(6+length(pass2b))] <- pass2b
      pass2c <- unlist(strsplit(x = pass1[3], split = "(\\;|\\.)"))
      if (length(pass2c) %% 4 == 0) {
        outvec[9:(8 + length(pass2c))] <- pass2c
      }
      pass2d <- unlist(strsplit(x = pass1[4], split = "\\;"))
      outvec[17:(16+length(pass2d))] <- pass2d
    }
    return(outvec)
  }, rep('NA', 18), USE.NAMES = FALSE))
  colnames(ctsplit) <- c(
    'TRAV_1', 'TRAJ_1', 'TRAC_1',
    'TRAV_2', 'TRAJ_2', 'TRAC_2',
    'TRA_nt_1', 'TRA_nt_2',
    'TRBV_1', 'TRBJ_1', 'None_1', 'TRBC_1',
    'TRBV_2', 'TRBJ_2', 'None_2', 'TRBC_2',
    'TRB_nt_1', 'TRB_nt_2')
  ctsplit[ctsplit == "NA"] <- NA # Correcting NAs
  sobj@meta.data <- cbind(sobj@meta.data, ctsplit)
  ### Adding length of TR sequence
  for(x in c("TRA_nt_1", "TRA_nt_2", "TRB_nt_1", "TRB_nt_2")) sobj@meta.data[paste0(x,"_len")] = nchar(sapply(sobj@meta.data[x], as.character))

  return(sobj)
}

### Spliting CTnt and CTgenes (into separate columns for IG-V/L and corresponding sequences) and save as metadata
### and Adding length of IG sequence to meta.data
split.bcr <- function(sobj=NULL){
  sobj@meta.data = tidyr::separate(sobj@meta.data, CTnt, c("cdr3_nt_IGH","cdr3_nt_IGL"), sep = "_", remove = FALSE)
  ctsplit <- t(vapply(sobj@meta.data$CTgene, function(x) {
    outvec <- rep('NA', 7)
    if (!is.na(x)) {
      pass1 <- unlist(strsplit(x = x, split = "\\_"))
      #IGH
      pass2a <- unlist(strsplit(x = pass1[1], split = "(\\.)"))
      outvec[1:length(pass2a)] <- pass2a
      #IGL
      pass2b <- unlist(strsplit(x = pass1[2], split = "(\\.)"))
      outvec[5:(4 + length(pass2b))] <- pass2b
    }
    return(outvec)
  }, rep('NA', 7), USE.NAMES = FALSE))
  colnames(ctsplit) <- c('IGHV', 'IGHJ', 'IGHD', 'Isotype',
                         'IGLV', 'IGLJ', 'IGLC')
  ctsplit[ctsplit == "NA"] <- NA # Correcting NAs
  sobj@meta.data <- cbind(sobj@meta.data, ctsplit)
  ### Adding length of TR sequence
  for(x in c("cdr3_nt_IGH", "cdr3_nt_IGL")) sobj@meta.data[paste0(x,"_len")] = nchar(sapply(sobj@meta.data[x], as.character))

  return(sobj)
}

## Frequency analysis
Freq.g <- function(sobj=NULL, out.dir = NULL, sample.name=NULL, reduction=NULL, freq_col="Frequency"){
  require(patchwork)
  require(dplyr)

  ## Frequency analysis
  ### UMAP of all frequencies
  png(paste0(out.dir,'/Frequency_umap',sample.name,'.png'), width = 800, height = 1000)
  print((Seurat::FeaturePlot(sobj, reduction = reduction, features = freq_col) + Seurat::DarkTheme() + ggplot2::ggtitle("")) /
          (Seurat::DimPlot(sobj, reduction = reduction, group.by = "cloneType") + Seurat::DarkTheme() + ggplot2::ggtitle("") ) +
          plot_annotation(title = "Clonotype frequency", theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))))
  dev.off()
  ### Calculation of clonotypes
  #All
  sequences = sobj@meta.data %>% select(.data[[freq_col]],CTaa) %>% distinct() %>% arrange(desc(.data[[freq_col]])) %>% na.omit()
  sobj <- scRepertoire::highlightClonotypes(sobj, cloneCall = "aa", sequence = sequences$CTaa)
  sobj$highlight_aa_all = as.character(sobj$highlight)
  sobj$highlight=NULL
  #Top 10 frequencies, and top 10 to top 20 frequencies
  top20_freq = sobj@meta.data %>% select(.data[[freq_col]],CTaa,highlight_aa_all) %>% distinct() %>% arrange(desc(.data[[freq_col]])) %>% na.omit() %>% top_n(n = 20, wt = .data[[freq_col]])
  top20_freq = top20_freq[1:20,]
  rownames(top20_freq)=top20_freq$highlight_aa_all
  sobj$highlight_aa_top10_freq <- ifelse(sobj$highlight_aa_all %in% top20_freq$highlight[1:10], sobj$highlight_aa_all, NA)
  sobj$highlight_aa_top11to20_freq <- ifelse(sobj$highlight_aa_all %in% top20_freq$highlight[11:length(top20_freq$highlight)], sobj$highlight_aa_all, NA)
  sobj$highlight_aa_top20_freq <- ifelse(sobj$highlight_aa_all %in% top20_freq$highlight[1:length(top20_freq$highlight)], sobj$highlight_aa_all, NA)
  
  #UMAP of top 10 frequencies
  png(paste0(out.dir,'/Frequency_top_10_umap',sample.name,'.png'), width = 800, height = (400+350))
  print(patchwork::wrap_elements( (Seurat::DimPlot(sobj, reduction = reduction, group.by = "highlight_aa_top10_freq")  + Seurat::DarkTheme()) / gridExtra::tableGrob(top20_freq[1:10,c(freq_col,"CTaa")], theme = gridExtra::ttheme_default(base_size = 10)) +
                                    plot_annotation(title = paste0("Top 10 Clonotypes (by frequencies)"),  subtitle = paste0("(",dim(sobj)[2]," cells)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0.5, face="bold"))) +
                                    plot_layout(heights = c(2, 1))))
  dev.off()
  #UMAP of top 11 to 20 frequencies
  png(paste0(out.dir,'/Frequency_top11to20_umap',sample.name,'.png'), width = 800, height = (400+350))
  print(patchwork::wrap_elements( (Seurat::DimPlot(sobj, reduction = reduction, group.by = "highlight_aa_top11to20_freq")  + Seurat::DarkTheme()) / gridExtra::tableGrob(top20_freq[11:length(top20_freq$CTaa),c(freq_col,"CTaa")], theme = gridExtra::ttheme_default(base_size = 10)) +
                                    plot_annotation(title = paste0("Top 11 to ", dim(top20_freq)[1], " Clonotypes (by frequencies)"),  subtitle = paste0("(",dim(sobj)[2]," cells)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0.5, face="bold"))) +
                                    plot_layout(heights = c(2, 1))))
  dev.off()
  return(sobj)
}

## Physicochemical properties of the CDR3
Physicochemical_properties.g <- function(sobj=NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, sample.name=NULL, type='TCR'){
  require(patchwork)
  require(ggplot2)
  ### Get cdr3 aa
  if(type == 'TCR'){
    CTaa_split <- as.data.frame(t(sapply(sobj@meta.data$CTaa, function(x) {
      outvec <- rep(NA, 4)
      if (!is.na(x)) {
        tmp <- unlist(strsplit(x = x, split = "\\_"))
        tmp1 <- unlist(strsplit(x = tmp[1], split = "(\\;)"))
        tmp2 <- unlist(strsplit(x = tmp[2], split = "(\\;)"))
        outvec[1:2]=tmp1[1:2]
        outvec[3:4]=tmp2[1:2]
      }
      return(outvec)
    }, USE.NAMES = FALSE)))
    colnames(CTaa_split)=c("TRA_cdr3_1", "TRA_cdr3_2", "TRB_cdr3_1", "TRB_cdr3_2")
  }else if(type == 'BCR'){
    CTaa_split <- as.data.frame(t(sapply(sobj@meta.data$CTaa, function(x) {
      outvec <- rep(NA, 2)
      if (!is.na(x)) outvec <- unlist(strsplit(x = x, split = "\\_"))
      return(outvec)
    }, USE.NAMES = FALSE)))
    colnames(CTaa_split)=c("IGH_cdr3", "IGL_cdr3")
  }
  CTaa_split[CTaa_split == "NA"] <- NA # Correcting NAs
  rownames(CTaa_split)=rownames(sobj@meta.data)
  ### Calculation of aa properties for each TRA/TRB or IG
  for (i in colnames(CTaa_split)){
    CTaa_split_global_withoutNA = CTaa_split[(!is.na(CTaa_split[i])),] #suppr NA values
    if(dim(CTaa_split_global_withoutNA)[1]!=0){
      aaProPerties <- alakazam::aminoAcidProperties(CTaa_split_global_withoutNA, seq=i, nt=F, trim=F, label=i)
      ### Merging results
      CTaa_split_global = if(type=='TCR') merge(CTaa_split, aaProPerties[,5:dim(aaProPerties)[2]], by="row.names", all.x=TRUE, all.y=TRUE) else merge(CTaa_split, aaProPerties[,3:dim(aaProPerties)[2]], by="row.names", all.x=TRUE, all.y=TRUE)
      rownames(CTaa_split_global) <- CTaa_split_global$Row.names
      CTaa_split_global$Row.names <- NULL
      summary=NULL
      names_properties = stringr::str_replace(grep("_aa_",colnames(aaProPerties), value=TRUE), ".+_aa_", "")
      for (j in names_properties){
        ### Table
        summary=rbind(summary,summary(CTaa_split_global[[paste0(i,"_aa_", j)]]) %>% round(2))
        ### Plots
        assign(paste0("plot_aaProperties_",j),patchwork::wrap_elements(ggplot(cbind(sobj@meta.data,CTaa_split_global), aes(x=orig.ident, y=.data[[paste0(i,"_aa_", j)]])) +
                                                                         geom_boxplot(aes(fill=orig.ident)) +
                                                                         ggtitle(paste0(ifelse(j=="gravy", "hydrophobicity", j))) + xlab("") + ylab("") + theme(legend.position='none')))
      }
      ### Formatting table
      rownames(summary)=names_properties
      ### Formatting Plot
      assign(paste0("plot_aaProperties_",i),patchwork::wrap_elements((plot_aaProperties_length | plot_aaProperties_gravy | plot_aaProperties_bulk | plot_aaProperties_aliphatic | plot_aaProperties_polarity | plot_aaProperties_charge | plot_aaProperties_basic | plot_aaProperties_acidic | plot_aaProperties_aromatic | gridExtra::tableGrob(summary, theme = gridExtra::ttheme_default(base_size = 10))) +
                                                                       patchwork::plot_annotation(title =  paste(unlist(strsplit(x = i, split = "_cdr3_")), collapse=" ")) +
                                                                       plot_layout(widths = c(1,1,1,1,1,1,1,1,1,2)) ) )
    }else{
      assign(paste0("plot_aaProperties_",i),patchwork::wrap_elements((patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer()) +
                                                                       patchwork::plot_annotation(title =  paste(unlist(strsplit(x = i, split = "_cdr3_")), collapse=" ")) +
                                                                       plot_layout(widths = c(1,1,1,1,1,1,1,1,1,2)) ) )
    }
  }
  ### Save
  png(paste0(out.dir,'/aaProperties.png'), width = 3500, height = 250*length(colnames(CTaa_split)))
  if(type == 'TCR'){
    print( (plot_aaProperties_TRA_cdr3_1 / plot_aaProperties_TRA_cdr3_2 / plot_aaProperties_TRB_cdr3_1 / plot_aaProperties_TRB_cdr3_2 ) +
             plot_annotation(title = sample.name, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) )
  }else if(type =='BCR'){
    print( (plot_aaProperties_IGH_cdr3 / plot_aaProperties_IGL_cdr3 ) +
             plot_annotation(title = sample.name, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) )
  }
  dev.off()
}

## Quantification of unique contig analysis
Quantif.unique.c <- function(sobj = NULL, ident.name=NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL, filtred_metadata_aa=NULL, filtred_metadata_nt=NULL, filtred_metadata_gene=NULL, filtred_metadata_gene_nt=NULL){
  require(patchwork)
  require(dplyr)
  df_nb_cells=as.data.frame(table(sobj[[ident.name]]), row.names = paste0("Cluster ", as.data.frame(table(sobj[[ident.name]]))$Var1))
  df_nb_cells$Var1=NULL
  colnames(df_nb_cells)="nb_cells"
  for(x in list_type_clT){
    ### Tables
    tmp = scRepertoire::quantContig(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall=x, scale = T, exportTable = T)
    tmp = data.frame(lapply((tmp %>% select("total_contigs_numbers"="total","unique_contig_numbers"="contigs","unique_contig_percents"="scaled") %>% round(2)), as.character), row.names = paste0("Cluster ", tmp$values))
    tmp = merge(df_nb_cells, tmp, by="row.names", all=TRUE) %>% arrange(desc(nb_cells))
    rownames(tmp)=tmp$Row.names
    tmp$Row.names=NULL
    colnames(tmp)[1]="total_nb_cells"
    assign(paste0("clust_quantUniqueContig_",sub("\\+","_",x)),tmp)
    sobj@misc$scRepertoire$clonalQuantifContig[[paste0("clust_quantUniqueContig_",sub("\\+","_",x))]]=get(paste0("clust_quantUniqueContig_",sub("\\+","_",x)))
    ### Plots
    assign(paste0("plot_clust_quantUniqueContig_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::quantContig(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall=x, scale = T) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.6, face="bold")))))
  }
  ### Save
  png(paste0(out.dir,'/clust_quantContig_',sample.name,'.png'), width =2100, height = 800)
  print(( (plot_clust_quantUniqueContig_gene_nt | plot_clust_quantUniqueContig_gene | plot_clust_quantUniqueContig_nt | plot_clust_quantUniqueContig_aa ) /
            gridExtra::grid.arrange(gridExtra::tableGrob(clust_quantUniqueContig_gene_nt, theme = gridExtra::ttheme_default(base_size = 10)) , gridExtra::tableGrob(clust_quantUniqueContig_gene, theme = gridExtra::ttheme_default(base_size = 10)) , gridExtra::tableGrob(clust_quantUniqueContig_nt, theme = gridExtra::ttheme_default(base_size = 10)) ,  gridExtra::tableGrob(clust_quantUniqueContig_aa, theme = gridExtra::ttheme_default(base_size = 10)), nrow=1) ) +
          plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) )
  dev.off()
  return(sobj)
}

## Clonal Homeostasis analysis
Homeo.c <- function(sobj = NULL, ident.name=NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL, filtred_metadata_aa=NULL, filtred_metadata_nt=NULL, filtred_metadata_gene=NULL, filtred_metadata_gene_nt=NULL){
  ### Tables
  for(x in list_type_clT){
    ### Tables
    scR_res = scRepertoire::clonalHomeostasis(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, exportTable = TRUE) %>% as.data.frame() %>% round(4)
    colnames(scR_res)=c("Rare","Small","Medium", "Large", "Hyperexpanded")
    assign( paste0("clust_clhomeo_",sub("\\+","_",x)), data.frame(lapply(scR_res, as.character), row.names = paste0("Cluster ", row.names(scR_res))))
    # sobj@misc$scRepertoire$clonalHomeostasis[[paste0("clust_clhomeo_",sub("\\+","_",x))]]=get(paste0("clust_clhomeo_",sub("\\+","_",x)))
    ### Plots
    assign(paste0("plot_clust_clhomeo_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::clonalHomeostasis(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.5, face="bold")))))
  }
  ### Save
  png(paste0(out.dir,'/clust_clhomeo',sample.name,'.png'), width =2000, height = 800)
  print(( (plot_clust_clhomeo_gene_nt | plot_clust_clhomeo_gene | plot_clust_clhomeo_nt | plot_clust_clhomeo_aa ) /
            gridExtra::grid.arrange(gridExtra::tableGrob(clust_clhomeo_gene_nt, theme = gridExtra::ttheme_default(base_size = 10)) , gridExtra::tableGrob(clust_clhomeo_gene, theme = gridExtra::ttheme_default(base_size = 10)) , gridExtra::tableGrob(clust_clhomeo_nt, theme = gridExtra::ttheme_default(base_size = 10)) ,  gridExtra::tableGrob(clust_clhomeo_aa, theme = gridExtra::ttheme_default(base_size = 10)), nrow=1) ) +
          plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) )
  dev.off()

  return(sobj)
}

## Clonal Proportions analysis
Prop.c <- function(sobj = NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL, filtred_metadata_aa=NULL, filtred_metadata_nt=NULL, filtred_metadata_gene=NULL, filtred_metadata_gene_nt=NULL){
  require(patchwork)
  require(dplyr)
  for(x in list_type_clT){
    ### Tables
    tmp = scRepertoire::clonalProportion(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, exportTable = TRUE)
    tmp = tmp %>% as.data.frame() %>% mutate_all(as.character) %>% as.data.frame(row.names = paste0("Cluster ", row.names(tmp)))
    assign( paste0("clust_clprop_",sub("\\+","_",x)), tmp)
    sobj@misc$scRepertoire$clonalProportion[[paste0("clust_clprop_",sub("\\+","_",x))]]=get(paste0("clust_clprop_",sub("\\+","_",x)))
    ### Plots
    assign(paste0("plot_clust_clprop_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::clonalProportion(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x) + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.5, face="bold")))))
  }
  ### Save
  png(paste0(out.dir,'/clust_clprop', sample.name, '.png'), width =2000, height = 800)
  print(( (plot_clust_clprop_gene_nt | plot_clust_clprop_gene | plot_clust_clprop_nt | plot_clust_clprop_aa ) /
       gridExtra::grid.arrange(gridExtra::tableGrob(clust_clprop_gene_nt, theme = gridExtra::ttheme_default(base_size = 10)) , gridExtra::tableGrob(clust_clprop_gene, theme = gridExtra::ttheme_default(base_size = 10)) , gridExtra::tableGrob(clust_clprop_nt, theme = gridExtra::ttheme_default(base_size = 10)) ,  gridExtra::tableGrob(clust_clprop_aa, theme = gridExtra::ttheme_default(base_size = 10)), nrow=1) ) +
      plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) )
  dev.off()
  return(sobj)
}

## Diversity analysis
Div.c <- function(sobj = NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL, filtred_metadata_aa=NULL, filtred_metadata_nt=NULL, filtred_metadata_gene=NULL, filtred_metadata_gene_nt=NULL){
  require(patchwork)
  require(dplyr)
  for(x in list_type_clT){
    ### Tables
    tmp = scRepertoire::clonalDiversity(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, group = "cluster", exportTable = TRUE)
    tmp[tmp=='NaN']=NA
    assign( paste0("clust_cldiv_",sub("\\+","_",x)), data.frame(lapply(( tmp %>% mutate_all(function(x){ as.numeric(as.character(x)) }) %>% round(2) %>% select("Shannon","Inv.Simpson","Chao","ACE") ), as.character), row.names = paste0("Cluster ", tmp$cluster)))
    # sobj@misc$scRepertoire$clonalDiversity[[paste0("clust_cldiv_",sub("\\+","_",x))]]=get(paste0("clust_cldiv_",sub("\\+","_",x)))
    ### Plots
    assign(paste0("plot_clust_cldiv_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::clonalDiversity(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, group = "cluster") + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.5, face="bold")))))
  }
  ### Save
  png(paste0(out.dir,'/clust_cldiv', sample.name, '.png'), width =2000, height = 800)
  print(( (plot_clust_cldiv_gene_nt | plot_clust_cldiv_gene | plot_clust_cldiv_nt | plot_clust_cldiv_aa ) /
            gridExtra::grid.arrange(gridExtra::tableGrob(clust_cldiv_gene_nt, theme = gridExtra::ttheme_default(base_size = 10)) , gridExtra::tableGrob(clust_cldiv_gene, theme = gridExtra::ttheme_default(base_size = 10)) , gridExtra::tableGrob(clust_cldiv_nt, theme = gridExtra::ttheme_default(base_size = 10)) ,  gridExtra::tableGrob(clust_cldiv_aa, theme = gridExtra::ttheme_default(base_size = 10)), nrow=1) ) +
          plot_annotation(title = sample.name, caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) )
  dev.off()
  return(sobj)
}

## Frequency analysis
Freq.c <- function(sobj = NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL, ident.name=NULL, reduction=NULL, freq_col="Frequency", filtred_metadata_aa=NULL, filtred_metadata_nt=NULL, filtred_metadata_gene=NULL, filtred_metadata_gene_nt=NULL){
  require(patchwork)
  require(dplyr)
  require(ggplot2)
  ### Histogramm of frequencies
  png(paste0(out.dir,'/Frequency_hist', sample.name, '.png'), width = 800, height = 1000)
  ggplot2::ggplot(sobj@meta.data, aes(x = get(ident.name), y = .data[[freq_col]])) +
    ggplot2::geom_boxplot(outlier.alpha = 0, aes(fill = get(ident.name))) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::theme_classic() +
    ggplot2::xlab("Clusters")
  dev.off()

  ### Calculation of clonotypes
  #Top 10 frequencies by clusters
  for(cluster in levels(Seurat::Idents(sobj))){
    if(!any(!is.na(subset(sobj, idents = cluster)@meta.data["highlight_aa_all"]))) next # if no clonotype for this cluster -> next
    #For all clonotype of the cluster
    sobj@meta.data[[paste0("highlight_aa_clust",cluster)]] <- ifelse(sobj@meta.data[[ident.name]]==cluster,sobj@meta.data$highlight_aa_all, NA)
    tmp = sobj@meta.data %>% add_count(.data[[paste0("highlight_aa_clust",cluster)]], name = paste0("Frequency_aa_clust",cluster))
    sobj@meta.data[[paste0("Frequency_aa_clust",cluster)]]=tmp[[paste0("Frequency_aa_clust",cluster)]] #NB: on ne peut pas rajouter la colone directement dans le meta.data sinon ça fait planter les plots pour une raison inconnue.
    sobj@meta.data[[paste0("Frequency_aa_clust",cluster)]] <- ifelse(is.na(sobj@meta.data[[paste0("highlight_aa_clust",cluster)]]), NA, sobj@meta.data[[paste0("Frequency_aa_clust",cluster)]])
    #For top 10 clonotypes of the cluster
    top10_freq = sobj@meta.data[,c(.data[[freq_col]], paste0("Frequency_aa_clust",cluster),"CTaa",paste0("highlight_aa_clust",cluster))] %>% distinct() %>% na.omit() %>% arrange(desc(.data[[paste0("Frequency_aa_clust",cluster)]])) %>% top_n(n = 10, wt = .data[[paste0("Frequency_aa_clust",cluster)]])
    if(dim(top10_freq)[1]>10) top10_freq=top10_freq[1:10,]
    rownames(top10_freq)=top10_freq[[paste0("highlight_aa_clust",cluster)]]
    sobj@meta.data$highlight_aa_top10_freq_cluster <- ifelse(sobj@meta.data[[paste0("highlight_aa_clust",cluster)]] %in% top10_freq[[paste0("highlight_aa_clust",cluster)]],sobj@meta.data[[paste0("highlight_aa_clust",cluster)]], NA)
    #UMAP of top 10 frequencies
    png(paste0(out.dir,'/Frequency_top_10_clust',cluster,'_umap', sample.name, '.png'), width = 600, height = 750)
    print(patchwork::wrap_elements( (Seurat::DimPlot(sobj, reduction = reduction, group.by = "highlight_aa_top10_freq_cluster") + Seurat::DarkTheme()) / gridExtra::tableGrob(top10_freq[,c(.data[[freq_col]], paste0("Frequency_aa_clust",cluster),"CTaa")], theme = gridExtra::ttheme_default(base_size = 10)) +
                                      plot_annotation(title = paste0("Top ", dim(top10_freq)[1], " Clonotypes (by frequencies)"), theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0.5, face="bold")))))
    dev.off()
  }
}

## Clonal Overlap analysis
Overlap.c <- function(sobj = NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL, filtred_metadata_aa=NULL, filtred_metadata_nt=NULL, filtred_metadata_gene=NULL, filtred_metadata_gene_nt=NULL){
  require(patchwork)
  require(dplyr)
  require(ggplot2)
  ### Plots
  for(x in list_type_clT) assign(paste0("plot_clust_clOverlap_",sub("\\+","_",x)),patchwork::wrap_elements(scRepertoire::clonalOverlap(get(paste0("filtred_metadata_", sub("\\+","_",x))), cloneCall = x, method="overlap") + plot_annotation(title = x, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=10, hjust=0.5, face="bold")))))
  ### Save
  png(paste0(out.dir,'/clust_clOverlap_',sample.name,'.png'), width = 2000, height = 600)
  print( (plot_clust_clOverlap_gene_nt | plot_clust_clOverlap_gene | plot_clust_clOverlap_nt | plot_clust_clOverlap_aa ) +
           plot_annotation(title = sample.name,subtitle ="It is looking at the overlap of clonotypes scaled to the smaller number of unique clonotypes in the compared clusters", caption = caption, theme = ggplot2::theme(plot.caption = ggplot2::element_text(hjust=0), plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) )
  dev.off()
  ### Complete table
  for(x in list_type_clT) {
    ### Clusters names identification
    name_clust = names(get(paste0("filtred_metadata_",sub("\\+","_",x))))
    ### Create dataframe for results
    overlap_res=data.frame(tested_clust=character(),
                           nb_unique_seq_tested_clust=double(),
                           ref_clust=character(),
                           nb_unique_seq_ref_clust=double(),
                           nb_overlap=double(),
                           list_sequences=character(),
                           stringsAsFactors=FALSE)
    ### Analyse
    if(x=="gene+nt") cloneCall="CTstrict" else cloneCall=paste0("CT",x)
    df = get(paste0("filtred_metadata_",sub("\\+","_",x)))
    for (i in name_clust){
      df.i <- df[[i]]
      df.i <- df.i[,c("barcode",cloneCall)]
      df.i_unique <- df.i[!duplicated(df.i$barcode),]
      nb_i_unique=dim(df.i_unique)[1]
      for (j in name_clust){
        if(i != j){
          df.j <- df[[j]]
          df.j <- df.j[,c("barcode",cloneCall)]
          df.j_unique <- df.j[!duplicated(df.j$barcode),]
          nb_j_unique=dim(df.j_unique)[1]
          overlap <- intersect(df.i_unique[,cloneCall], df.j_unique[,cloneCall])
          nb_overlap = length(overlap)
          if(nb_overlap==0) overlap='NA'
          overlap_res=rbind(overlap_res, data.frame(tested_clust=i,
                                                    nb_unique_seq_tested_clust=nb_i_unique,
                                                    tested_clust_list_sequences=paste(df.i_unique[,cloneCall], collapse=";"),
                                                    ref_clust=j,
                                                    nb_unique_seq_ref_clust=nb_j_unique,
                                                    ref_clust_list_sequences=paste(df.j_unique[,cloneCall], collapse=";"),
                                                    nb_overlap=nb_overlap,
                                                    overlap_list_sequences=paste(overlap, collapse=";"),
                                                    stringsAsFactors=FALSE))
        }
      }
    }
    ### Save
    write.table(overlap_res, file = paste0(out.dir,"/overlap_",sub("\\+","_",x),"_",sample.name,".txt"), quote = FALSE, sep = "\t", na = "NA", row.names = FALSE, col.names = TRUE, append=FALSE, dec=",")
  }
}

### Get cdr3 aa
Physicochemical_properties.c <- function(sobj = NULL, list_type_clT = c("gene+nt","gene","nt","aa"), out.dir = NULL, caption=NULL, sample.name=NULL, ident.name=NULL, filtred_metadata_aa=NULL, filtred_metadata_nt=NULL, filtred_metadata_gene=NULL, filtred_metadata_gene_nt=NULL, type='TCR'){
  require(patchwork)
  require(ggplot2)
  ### Get cdr3 aa
  if(type == 'TCR'){
    CTaa_split <- as.data.frame(t(sapply(sobj@meta.data$CTaa, function(x) {
      outvec <- rep(NA, 4)
      if (!is.na(x)) {
        tmp <- unlist(strsplit(x = x, split = "\\_"))
        tmp1 <- unlist(strsplit(x = tmp[1], split = "(\\;)"))
        tmp2 <- unlist(strsplit(x = tmp[2], split = "(\\;)"))
        outvec[1:2]=tmp1[1:2]
        outvec[3:4]=tmp2[1:2]
      }
      return(outvec)
    }, USE.NAMES = FALSE)))
    colnames(CTaa_split)=c("TRA_cdr3_1", "TRA_cdr3_2", "TRB_cdr3_1", "TRB_cdr3_2")
  }else if(type =='BCR'){
    CTaa_split <- as.data.frame(t(sapply(sobj@meta.data$CTaa, function(x) {
      outvec <- rep(NA, 2)
      if (!is.na(x)) outvec <- unlist(strsplit(x = x, split = "\\_"))
      return(outvec)
    }, USE.NAMES = FALSE)))
    colnames(CTaa_split)=c("IGH_cdr3", "IGL_cdr3")
  }
  CTaa_split[CTaa_split == "NA"] <- NA # Correcting NAs
  rownames(CTaa_split)=rownames(sobj@meta.data)
  CTaa_split_clust=cbind(sobj[[ident.name]],CTaa_split)
  ### Calculation of aa properties for each TRA and TRB
  for (i in colnames(CTaa_split_clust)[-1]){
    CTaa_split_clust_withoutNA = CTaa_split_clust[(!is.na(CTaa_split_clust[i])),] #suppr NA values
    if(dim(CTaa_split_clust_withoutNA)[1]!=0){
      aaProPerties <- alakazam::aminoAcidProperties(CTaa_split_clust_withoutNA, seq=i, nt=F, trim=F, label=i)
      ### Merging results
      CTaa_split_clust = if(type=='TCR') merge(CTaa_split_clust, aaProPerties[,6:dim(aaProPerties)[2]], by="row.names", all.x=TRUE, all.y=TRUE) else merge(CTaa_split_clust, aaProPerties[,4:dim(aaProPerties)[2]], by="row.names", all.x=TRUE, all.y=TRUE)
      rownames(CTaa_split_clust) <- CTaa_split_clust$Row.names
      CTaa_split_clust$Row.names <- NULL
      ### Plots
      plot_clust_aaProperties_length=patchwork::wrap_elements(ggplot(CTaa_split_clust, aes(x=.data[[ident.name]], y=.data[[paste0(i,"_aa_length")]])) +
                                                                geom_boxplot(aes(fill=.data[[ident.name]])) +
                                                                ggtitle("length") + xlab("") + ylab("") + theme(legend.position='none'))
      names_properties = stringr::str_replace(grep("_aa_",colnames(aaProPerties), value=TRUE), ".+_aa_", "")
      for (j in names_properties){
        assign(paste0("plot_clust_aaProperties_",j),patchwork::wrap_elements(ggplot(CTaa_split_clust, aes(x=.data[[ident.name]], y=.data[[paste0(i,"_aa_", j)]])) +
                                                                               geom_boxplot(aes(fill=.data[[ident.name]])) +
                                                                               ggtitle(paste0(ifelse(j=="gravy", "hydrophobicity", j))) + xlab("") + ylab("") + theme(legend.position='none')))
      }
      ### Formatting Plot
      assign(paste0("plot_clust_aaProperties_",i),patchwork::wrap_elements((plot_clust_aaProperties_length | plot_clust_aaProperties_gravy | plot_clust_aaProperties_bulk | plot_clust_aaProperties_aliphatic | plot_clust_aaProperties_polarity | plot_clust_aaProperties_charge | plot_clust_aaProperties_basic | plot_clust_aaProperties_acidic | plot_clust_aaProperties_aromatic) +
                                                                             patchwork::plot_annotation(title =  paste(unlist(strsplit(x = i, split = "_cdr3_")), collapse=" ")) ) )
    }else{
      assign(paste0("plot_clust_aaProperties_",i),patchwork::wrap_elements((patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer() | patchwork::plot_spacer()) +
                                                                             patchwork::plot_annotation(title =  paste(unlist(strsplit(x = i, split = "_cdr3_")), collapse=" ")) ) )
    }
  }
  ### Save
  png(paste0(out.dir,'/aaProperties_', sample.name, '.png'), width =2500, height = 1000)
  if(type == 'TCR'){
    print( (plot_clust_aaProperties_TRA_cdr3_1 / plot_clust_aaProperties_TRA_cdr3_2 / plot_clust_aaProperties_TRB_cdr3_1 / plot_clust_aaProperties_TRB_cdr3_2 ) +
             plot_annotation(title = sample.name, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) )
  }else if(type =='BCR'){
    print( (plot_clust_aaProperties_IGH_cdr3 / plot_clust_aaProperties_IGL_cdr3 ) +
             plot_annotation(title = sample.name, theme = ggplot2::theme(plot.title = ggplot2::element_text(size=20, hjust=0, face="bold"))) )
  }
  dev.off()
}

### Get references pof packages used to compute results thancks to Materials and Methods slot
find_ref <- function(MandM = NULL, pipeline.path=NULL){
  if(is.null(pipeline.path)) stop("pipeline.path not set!")
  df_ref = readxl::read_excel(paste0(pipeline.path, "/resources/PACKAGES_REFRENCES/references.xlsx"))
  df_res = data.frame(package_name=character(), reference=character())
  for (i in 1:length(df_ref$package_name)){
    for (j in names(MandM)){
      if ((grepl(df_ref$package_name[i], MandM[j], ignore.case = TRUE)) && (! df_ref$package_name[i] %in% df_res$package_name)) df_res <- rbind(df_res, df_ref[i,])
    }
  }
  res=paste0(df_res$package_name, " : ", df_res$reference)
  return(res)
}

#Write Material and Method into texte file
write_MandM <- function(sobj=NULL, output.dir=NULL){
  MandM = c("Materials and Methods","")
  for (pipeline_part in names(sobj@misc$parameters$Materials_and_Methods)){
    if(pipeline_part != "References_packages"){
      if(pipeline_part == "part0_Alignment") pipeline_name <- c("","QC reads, Pseudo-mapping and quantification")
      if(pipeline_part == "part1_Droplets_QC") pipeline_name <- c("","QC data on each sample")
      if(pipeline_part == "part2_Filtering") pipeline_name <- NULL
      if(pipeline_part == "part3_Norm_DimRed_Eval") pipeline_name <- c("","Individual analysis")
      if(pipeline_part == "part4_Clust_Markers_Annot") pipeline_name <- NULL
      if(pipeline_part == "ADT") pipeline_name <- c("","Cell surface proteins (CITE-seq ADT)")
      if(pipeline_part == "TCR") pipeline_name <- c("","Single-cell immune profiling (TCR)")
      if(pipeline_part == "BCR") pipeline_name <- c("","Single-cell immune profiling (Ig)")
      if(pipeline_part == "TCR/BCR") pipeline_name <- c("","Single-cell immune profiling (TCR/Ig)")
      if(pipeline_part == "Cerebro") pipeline_name <- c("","Cerebro")
      
      ### TO DO: Add the same if for integrated and grouped analysis ###
      if(pipeline_part == "Integration_Norm_DimRed_Eval") pipeline_name <- c("","Integration analysis")
      if(pipeline_part == "Integration_Clust_Markers_Annot") pipeline_name <- NULL
      if(pipeline_part == "Grouped_analysis_Norm_DimRed_Eval") pipeline_name <- c("","Grouped analysis")
      if(pipeline_part == "Grouped_analysis_Clust_Markers_Annot") pipeline_name <- NULL
    
      MandM = c(MandM,pipeline_name,sobj@misc$parameters$Materials_and_Methods[[pipeline_part]])
    }
  }
  #add packages reference publications
  if("References_packages" %in% names(sobj@misc$parameters$Materials_and_Methods)) MandM = c(MandM,"","References of tools",sobj@misc$parameters$Materials_and_Methods$References_packages)
  #save
  file<-file(paste0(output.dir,"/Materials_and_Methods.txt"))
  writeLines(MandM, file)
  close(file)
}

