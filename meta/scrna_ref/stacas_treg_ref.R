## sara Tumor/LN KO and WT samples
# Core
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm10)
library(monocle3)
library(BiocParallel)
library(GenomicFeatures)
library(tidyr)
library(Seurat)
library(dplyr)
library(plyr)
library(reshape2)
library(stringr)
library(SeuratWrappers)
library(STACAS)
# Visualization
library(cowplot)
library(ggrastr)
library(ggplot2)
# Annotation
library(SingleR)
library(ProjecTILs)
# QC 
library(scater)
library(DoubletFinder)

###############
#### Setup ####

visualize <- FALSE
seed <- 1234
set.seed(seed)

PDIR <- '/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/treg_ref'
datadir <- file.path(PDIR, "data")
outdir <- file.path(PDIR, "results")
dir.create(datadir, recursive = T, showWarnings = F)
dir.create(outdir, recursive = T, showWarnings = F)
setwd(PDIR)


# Create ENZ -> SYMBOL mapping key
# genome_gse <- org.Mm.eg.db
# txby <- keys(genome_gse, 'SYMBOL')

# Read in gtf file to map ENS->Biotype
gtf_file <- '/cluster/projects/mcgahalab/ref/genomes/mouse/GRCm38/GTF/genome.gtf'
GTF <- rtracklayer::import(gtf_file)
ens2biotype_ids <- with(GTF, setNames(gene_biotype, gene_id))

###########################
#### Global Parameters ####
doublet_quantile_cutoff <- 0.95
visualize_qc <- FALSE

mt_miqc <- FALSE # whether to use miQC for percent.mt filter, else use isOutlier
call_doublets <- FALSE # whether to call doublets

# Preprocessed datasets:
#     Must contain a '*counts.csv.gz' file and '*metadata.csv.gz' file
#     in the /cluster/projects/mcgahalab/ref/scrna/projectils/gse_dat/GSEXXXX/
#     directory
pp_datasets <- c('GSE161345') #, 'GSE184423')

# Annotation
atlas_refdir <- '/cluster/projects/mcgahalab/ref/scrna/projectils'
scgate_refdir <- '/cluster/projects/mcgahalab/ref/scrna/projectils/scgate_model'
projectils_refdir <- file.path(atlas_refdir, "ProjecTILs")

###################
#### Functions ####
##-- Recurring functions ----
source("~/git/mini_projects/mini_functions/iterateMsigdb.R")
source("~/git/mini_projects/mini_functions/gsea2CytoscapeFormat.R")
source("~/git/mini_projects/mini_functions/geneMap.R")
source("~/git/mini_projects/mini_functions/wgcnaComplexHeatmap.R")
source("~/git/mini_projects/mini_functions/linearModelEigen.R")
source("~/git/mini_projects/mini_functions/singlecell/scrna.R")
source("~/git/mini_projects/mini_functions/singlecell/STACAS_fix.R")
source("/cluster/projects/mcgahalab/data/mcgahalab/st2_il33/scrna_tdln_tumor/scripts/doubletFinder_v4.R")

###########################
#### 0. Global objects ####
# Reads in a metadata file containing the GSE, GSM, and dataset label
# to be used when annotating the datasets downstream
ref_seu_dir <- file.path(atlas_refdir, "gse_dat", "Treg_datasets")
ref_celltypes <- read.csv(file.path(ref_seu_dir, "celltypes.csv"), 
                          check.names = F, stringsAsFactors = F) %>%
  tibble::column_to_rownames("GSM")

#######################################
#### 1.a Make ProjecTILs reference ####
# This section of the code reads in all the 10x cellranger
# files and loads them into seurat objects. It saves an
# RDS object which contains a list of seurat objects, one
# for each dataset
ref_groups <- list.files(ref_seu_dir, pattern = "^GSM")
outrds <- file.path(datadir, "seurat_obj", "ref1_unfiltered.rds")
relabel_id <- gsub("^GSM[0-9]*_", "", ref_groups) %>%
  setNames(., ref_groups)

if(file.exists(outrds)){
  cat("Reading in existing 'seus' object...\n")
  seus <- readRDS(file=outrds)
} else {
  seus <- NULL
}

seus <- lapply(ref_groups, function(grp){
  if(relabel_id[grp] %in% names(seus)){
    cat(paste0("Existing ID: ", grp, "\n"))
    seu <- seus[[relabel_id[grp]]]
  } else {
    gsm_id <- gsub("_.*", "", grp)
    print(paste0(grp, "..."))
    seu <- tryCatch({
      mtx <- Read10X(data.dir = file.path(ref_seu_dir, grp), 
                     strip.suffix=TRUE)
      if(class(mtx)=='list') mtx <- mtx[['Gene Expression']]
      CreateSeuratObject(counts = mtx, project = gsm_id)
    }, error=function(e){cat("Failed\n"); NULL})
  }
  return(seu)
})
names(seus) <- relabel_id
seus <- seus[which(!sapply(seus, is.null))]

# For each sampel, drop all cells where there are fewer than 50 genes detected
sapply(seus,dim) %>% as.data.frame.matrix %>% t
seus <- lapply(seus, function(seu_i){
  subset(seu_i, subset = nFeature_RNA > 50)
})
sapply(seus,dim) %>% as.data.frame.matrix %>% t

saveRDS(seus, file=outrds)

#####################################################################
#### 1.b OUTDATED - Load in preprocessed/packaged seurat objects ####
# This section of the code loads in all the data that was
# released on GEO in a non-standard format (i.e. not features,
# barcodes, matrix.mtx.gz files). Since each dataset 
# had specifics on how they released their data, each 
# dataset was manually hardcoded to address and standardize
# the formatting.
outrds <- file.path(datadir, "seurat_obj", "ref1b_filtered.rds")

if(file.exists(outrds)){
  cat("Reading in existing 'seus' object...\n")
  seus_preprocessed <- readRDS(file=outrds)
} else {
  seus_preprocessed <- NULL
}

.relabel_origident <- function(seu, map, colid='orig.ident'){
  idx <- which(seu@meta.data[,colid] %in% names(map))
  seu$tmp <- NA
  seu$tmp[idx] <- as.character(map[seu@meta.data[idx,colid]])
  seu$orig.ident <- seu$tmp
  seu$tmp <- NULL
  return(seu)
}
.filterGSE168158 <- function(){
  cat("Fixing GSE168158... \n")
  seu_gse216166_map <- c('GSE168158_1'='GSM5130034',
                         'GSE168158_2'='GSM51300342')
  
  seu_ppi <- .relabel_origident(seu_ppi, seu_gse216166_map)
  seu_ppi
}
.filterGSE216166 <- function(seu_ppi){
  cat("Fixing GSE216166...\n")
  seu_gse216166_map <- c('S1163175'='GSM6659730',
                         'S1163177'='GSM6659732',
                         'S1163179'='GSM6659734')
  
  Idents(seu_ppi) <- 'tissue_detail'
  seu_ppi <- subset(seu_ppi, idents='naive')
  seu_ppi$percent.mt <- seu_ppi$percent.mito
  seu_ppi <- .relabel_origident(seu_ppi, seu_gse216166_map)
  seu_ppi
}
.filterGSE117131 <- function(seu_ppi, ident_keep='Healthy'){
  cat("Fixing GSE117131...\n")
  seu_gse117131_map <- c('Mouse_Healthy1_001'='GSM3271853',
                         'Mouse_Healthy1_002'='GSM3271859',
                         'Mouse_Healthy2_001'='GSM3271855',
                         'Mouse_Healthy2_002'='GSM3271861',
                         'Mouse_Healthy3_001'='GSM3271857',
                         'Mouse_Healthy3_002'='GSM3271863',
                         'Mouse_Tumor1_001'='GSM3271852', 
                         'Mouse_Tumor1_002'='GSM3271858',
                         'Mouse_Tumor2_001'='GSM3271854',
                         'Mouse_Tumor2_002'='GSM3271860',
                         'Mouse_Tumor3_001'='GSM3271856',
                         'Mouse_Tumor3_002'='GSM3271862')
  
  seu_ppi$tissue_detail <- gsub("[0-9]_[0-9]*$", "", seu_ppi$Sample.Name) %>%
    gsub("Mouse_", "", .)
  Idents(seu_ppi) <- 'tissue_detail'
  # seu_ppi <- subset(seu_ppi, idents=ident_keep)
  seu_ppi <- PercentageFeatureSet(seu_ppi, pattern = "^mt-", col.name = "percent.mt")
  seu_ppi <- .relabel_origident(seu_ppi, seu_gse117131_map, colid='Sample.Name')
  seu_ppi
}

seus_preprocessed <- lapply(pp_datasets, function(pp_i){
  if(pp_i %in% names(seus_preprocessed)){
    cat(paste0("Using existing preprocessed data for ", pp_i, "\n"))
    seu_ppi <- seus_preprocessed[[pp_i]]
  } else {
    cat(paste0("Loading in count data and preprocessing ", pp_i, "\n"))
    dirx <- file.path(atlas_refdir, "gse_dat", pp_i)
    cts_f <- grep("counts.(txt|tsv|csv).gz$", list.files(dirx), value=T, perl=T)
    meta_f <- grep("(metadata|annotation).(txt|tsv|csv).gz$", list.files(dirx), value=T, perl=T)
    
    if(grepl("txt|tsv", cts_f)){
      cts <- read.table(file.path(dirx, cts_f), sep="\t", header = T)
    } else {
      cts <- read.csv(file.path(dirx, cts_f))
    }
    if(grepl("txt|tsv", meta_f)){
      meta <- read.table(file.path(dirx, meta_f), sep="\t", header = T)
    } else {
      meta <- read.csv(file.path(dirx, meta_f))
    }
    
    if('Barcode' %in% colnames(meta)){
      meta <- meta  %>% 
        tibble::column_to_rownames('Barcode')
    }
    if('X' %in% colnames(cts)){
      cts <- cts  %>% tibble::column_to_rownames('X')
    }
    if('gex_symbol' %in% colnames(cts)){
      cts <- cts  %>% tibble::column_to_rownames('gex_symbol')
    }
    
    seu_ppi <-   CreateSeuratObject(
      counts = cts, project = pp_i, 
      assay = "RNA", names.field = 1,
      names.delim = "_", meta.data = meta,
    )
    seu_ppi <- switch(pp_i,
                      GSE216166=.filterGSE216166(seu_ppi),
                      GSE117131=.filterGSE117131(seu_ppi),
                      GSE168158=.filterGSE168158(seu_ppi),
                      seu_ppi)
    
    seu_ppi <- preprocessSeu(
      seu_ppi, ncount_min=1000, ncount_max=Inf, 
      nfeature_min=0, nfeature_max=Inf, 
      mt_max=Inf, org='mouse', numpcs=30, getPCs=FALSE
    )
  }
  return(seu_ppi)
})
names(seus_preprocessed) <- pp_datasets

saveRDS(seus_preprocessed, file=outrds)

###########################################
#### 2. Preprocess to remove doublets  ####
dir.create(file.path(outdir, "gex", "qc"), showWarnings = F, recursive = T)
dir.create(file.path(outdir, "seurat_obj", "per_sample"), showWarnings = F,)
if(!exists("seus")) seus <- readRDS(file.path(datadir, "seurat_obj", "ref1_unfiltered.rds"))

outrds_mnn <- file.path(datadir, "seurat_obj", "ref2_dubfiltered.mnn.rds")
outrds_stacas <- file.path(datadir, "seurat_obj", "ref2_dubfiltered.stacas.rds")
outrds_int <- file.path(datadir, "seurat_obj", "ref2_dubfiltered.intermediate.rds")

### Preprocess and Doublet Calling
names(seus) <- as.character(sapply(seus, function(i) unique(i$orig.ident)))
seus <- lapply(seus, function(seu_i){
  PercentageFeatureSet(seu_i, pattern = "^mt-", col.name = "percent.mt")
})
# seu <- merge(seus[[1]], y = seus[-1], add.cell.ids = names(seus), project = "myeloid_reference")

#.... 1) Identify mitochondrial/nCount/nFeature outliers: ----
outliers <- lapply(seus, function(seu_i){
  # 3-MAD/Mindiff from Median as outliers
  # mitochondrial outlier: scater::isOutlier
  mt_outlier_scuttle <- isOutlier(seu_i$percent.mt, nmads=3, 
                                  min.diff=(10-median(seu_i$percent.mt)), 
                                  batch=seu_i$orig.ident, share.mads=T, 
                                  share.medians=T, type="higher")
  # nCount outlier: scater::isOutlier
  ncnt_outlier_scuttle <- isOutlier(seu_i$nCount_RNA, nmads=4, 
                                    # min.diff=median(seu_i$nCount_RNA), 
                                    batch=seu_i$orig.ident, share.mads=F, 
                                    share.medians=F, type="higher")
  # nFeature outlier: scater::isOutlier
  nft_outlier_scuttle <- isOutlier(seu_i$nFeature_RNA, nmads=3, 
                                   # min.diff=median(seu$nFeature_RNA), 
                                   batch=seu_i$orig.ident, share.mads=F, 
                                   share.medians=F, type="higher")
  
  nmt.max <- min(seu_i$percent.mt[mt_outlier_scuttle]) # 6530
  nfeature.max <- min(seu_i$nFeature_RNA[nft_outlier_scuttle]) # 6530
  ncnt.max <- min(seu_i$nCount_RNA[ncnt_outlier_scuttle]) # 37252
  c('mt'=max(c(nmt.max, 10)), 
    'nft'=max(c(nfeature.max, 5000)), 
    'ncnt'=max(c(ncnt.max, 17000)))
})

pdf(file.path(outdir, "gex", "qc", "ref_qc_featureScatters.pdf"))
# pdf(file.path("~/xfer", "ref_qc_featureScatters.raw.pdf"))
fsl <- lapply(seus, function(seu_i){
  fs1 <- FeatureScatter(seu_i, 'nCount_RNA', 'nFeature_RNA', shuffle=TRUE, raster=TRUE)
  fs2 <- FeatureScatter(seu_i, 'nCount_RNA', 'percent.mt', shuffle=TRUE, raster=TRUE)
  cowplot::plot_grid(fs1, fs2, ncol=2)
})
print(fsl)
dev.off()

#.... 2) Flag Potential Doublets ----
if(file.exists(outrds_int)) {
  seus_post <- readRDS(outrds_int)
  # names(seus_post) <- names(seus)
} else{
  seus_post <- NULL
}

#a) Preprocess for Doublet Calling
s.features <- stringr::str_to_title(cc.genes$s.genes)
g2m.features <- stringr::str_to_title(cc.genes$g2m.genes)

seu_ids <- setdiff(names(seus), names(seus_post))
seus_pre <- lapply(seu_ids, function(seu_id){
  if(seu_id %in% names(seus_post)){
    cat(paste0("Loading in existing preprocessed data: ", seu_id, "\n"))
    seu_i <- seus_post[[seu_id]]
  } else {
    cat(paste0("Preprocessing and adding on ", seu_id, "\n"))
    seu_i <- seus[[seu_id]]
    
    gsm_id <- as.character(unique(seu_i$orig.ident))
    outlier_i <- outliers[[gsm_id]]
    seu_i <- preprocessSeu(seu_i, ncount_min=1000, ncount_max=Inf, 
                           nfeature_min=50, nfeature_max=outlier_i['nft'], 
                           mt_max=outlier_i['mt'], org='mouse', numpcs=30, getPCs=FALSE)
  }
  return(seu_i)
})
names(seus_pre) <- seu_ids

#b) Doublet Calling
if(call_doublets){
  cat("Calling doublets...\n")
  seus_pre2 <- lapply(seus_pre, function(seu_i){
    ## DoubletDetection
    sce_cxds <- runCxds(seu_i, n=25)  # Looks for genes that are mutually exclusive, scores their expression
    seu_i$cxds_z <- sce_cxds$cxds$z_cxds
    seu_i$cxds_doublet <- (seu_i$cxds_z>3)
    
    sce_bcds <- runBcds(seu_i, n=25) # Creates artificial doublets, trains a classifier to predicts them
    seu_i$bcds_score <- sce_bcds$bcds[,1]
    seu_i$bcds_doublet <- (seu_i$bcds_score>0.95)
    
    seu_i
  })
  seus_post <- c(seus_post, seus_pre2)
} else {
  cat("Skipping doublet calling...\n")
  new_ids <- c(names(seus_post), seu_ids)
  seus_post <- c(seus_post, seus_pre)
}
names(seus_post) <- new_ids
saveRDS(seus_post, file=outrds_int)

######################################
#### 3. Identify Treg populations ####
seu <- seus_post[[1]]
recode_map <- function(x, grp=NULL, anno=TRUE){
  x %>% 
    recode("0"='Central Tregs',
           "1"='Central Tregs',
           "2"='Effector Tregs',
           "3"='Effector Tregs',
           "4"='Central Tregs',
           "5"='Effector Tregs',
           "6"='Central Tregs',
           "7"='NLT_like_Tregs',
           "8"='Effector Tregs',
           "9"='NLT_like_Tregs',
           "10"='Effector Tregs',
           "11"='NLT_like_Tregs',
           "12"='STAT1_Tregs',
           "13"='NA',
           "14"='NA')
}

# Features and heatmap selected to annotate based on Fig. 5c of 
# https://doi.org/10.1038/s41423-020-00599-z
feats <- c('Ifit1', 'Ifit3', 'Rsad2',
           'Sell', 'Ccr7', 'Bcl2', 
           'Tnfrsf9', 'Tnfrsf18', 'Izumo1r',
           'S100a4', 'S100a6', 'Ccr2')
pdf("~/xfer/identify_clusters.pdf")
DimPlot(seu, group.by='seurat_clusters', raster=T, reduction='umap', label=T)
DoHeatmap(seu, features=feats, group.by = 'seurat_clusters')
VlnPlot(seu, features=feats, group.by='seurat_clusters')
dev.off()

seu$anno <- seu$seurat_clusters %>% recode_map
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 1500, verbose = FALSE)
ndim = 20
seed = 1234
seu <- ScaleData(seu, verbose = TRUE)

seu <- RunPCA(seu, features = VariableFeatures(seu), ndims.print = 1:5,
                 nfeatures.print = 5, npcs = ndim)
seu <- RunUMAP(seu, reduction = "pca", dims = 1:ndim, seed.use = seed)


ref.treg <- ProjecTILs::make.reference(ref = seu, ndim = ndim, seed = seed, 
                                       recalculate.umap = TRUE, annotation.column = "anno")
pdf("~/xfer/identify_clusters.pdf")
DimPlot(ref.treg, group.by='functional.cluster', raster=T, reduction='umap', label=T)
DimPlot(ref.treg, group.by='seurat_clusters', raster=T, reduction='umap', label=T)
dev.off()

saveRDS(ref.treg, file=file.path(projectils_refdir, "seu.treg.rds"))
